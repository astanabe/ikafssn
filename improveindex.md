# ikafssnindex 高速化・縮重塩基対応計画

## 1. 目的

`ikafssnindex` のk-merスキャンを高速化するとともに、k-mer内に縮重塩基を1文字だけ含む場合はその縮重塩基を展開して複数のk-merとしてインデックスに登録する機能を追加する。

### 高速化

現在の実装では、BLAST DBから `GetAmbigSeq()` でNCBI NA8形式（1バイト/塩基）を取得し、`ncbi_na8_to_char()` でASCII文字に変換して `std::string` を構築し、`KmerScanner::scan()` 内で再度 `encode_base()` により2ビットコードに変換している。

```
BLAST DB (ncbi2na packed, 4塩基/バイト)
  ↓ GetAmbigSeq()     -- Toolkit内部でアンパック＋ambiguity適用
NCBI NA8 (1バイト/塩基)
  ↓ ncbi_na8_to_char() -- switch文による変換
ASCII文字列 std::string
  ↓ encode_base()      -- 256要素LUTによる変換
2ビットコード → k-merスライディングウィンドウ
```

これを以下に置き換える：

```
BLAST DB
  ↓ GetSequence()      -- mmapポインタ直接取得（コピーなし）
ncbi2na packed (4塩基/バイト) + ambiguityデータ
  ↓ 2ビットずつ直接シフト抽出 + ambiguityチェック
k-mer uint16/uint32（変換操作なし）
```

主な改善点：
- `std::string` アロケーション排除（配列あたり数KB〜数MBの動的メモリ確保がゼロに）
- NA8→char→2ビットの二重変換排除
- ncbi2naはikafssn内部の2ビットコードと同一エンコーディング（A=00, C=01, G=10, T=11）なので変換テーブル不要

### 縮重塩基対応

現在はNを含むk-merを全てスキップしている。改善後は、k-mer内に縮重塩基が1箇所だけ含まれる場合、その縮重塩基が表す全ての塩基に展開して複数のk-merを生成し、全てインデックスに登録する。N以外の全てのIUPAC縮重コードにも対応する。

k-mer内に2箇所以上の縮重塩基が含まれる場合は、展開数が最大16（2塩基縮重×2塩基縮重）以上になりコストが高いため、従来通りスキップする。

## 2. BLAST DB内部のncbi2naフォーマット

### 2.1 ncbi2naエンコーディング

| 2ビット値 | 塩基 |
|:---:|:---:|
| 00 | A |
| 01 | C |
| 10 | G |
| 11 | T |

ikafssn内部の `kmer_encoding.hpp` と完全に一致する。変換テーブル不要。

### 2.2 バイトパッキング

1バイトに4塩基を格納。MSB側から順に格納：

```
ビット:  7 6 | 5 4 | 3 2 | 1 0
塩基:  base0 |base1|base2|base3
```

例: バイト `0xE4` = `11 10 01 00` = T G C A

### 2.3 GetSequence() の戻り値

`CSeqDB::GetSequence(int oid, const char** buffer)` は、`.nsq` ファイルのmmap領域への直接ポインタを返す。戻り値はバイト長。

末尾バイトの下位2ビットは、そのバイト内の有効塩基数（remainder）を符号化している：

```
remainder = last_byte & 0x03;
total_bases = (num_bytes - 1) * 4 + remainder;
```

ただし `remainder == 0` は「追加塩基なし（配列長が4の倍数）」を意味する。この場合、末尾バイトは塩基データを持たないパディングバイトである。

**実装上は `seq_length(oid)` で塩基数を別途取得して使い、末尾バイトの処理は塩基数ベースで行う方が安全。**

### 2.4 ncbi2na内の縮重領域

ncbi2naは2ビット/塩基のため縮重塩基を直接表現できない。縮重位置にはランダムに選ばれた非縮重塩基（A/C/G/Tのいずれか）が格納されている。真の縮重情報は別途ambiguityデータとして付加される。

## 3. Ambiguityデータのフォーマット

### 3.1 取得方法

**方法A: `CSeqDBExpert::GetRawSeqAndAmbig()`**

```cpp
#include <objtools/blast/seqdb_reader/seqdbexpert.hpp>

CSeqDBExpert db(db_path, CSeqDB::eNucleotide);
const char* buffer = nullptr;
int seq_length = 0;
int ambig_length = 0;
db.GetRawSeqAndAmbig(oid, &buffer, &seq_length, &ambig_length);
// buffer[0..seq_length-1]       = ncbi2naパックデータ
// buffer[seq_length..seq_length+ambig_length-1] = ambiguityデータ
db.RetSequence(&buffer);
```

ncbi2naデータとambiguityデータの両方を一度に取得できる。ただし `CSeqDBExpert`（`CSeqDB`のサブクラス）を使用する必要がある。

**方法B: `CSeqDB::GetSequence()` + `x_GetAmbChar()` 相当の自前パース**

`GetSequence()` はncbi2naデータのみを返す。ambiguityデータは `.nsq` ファイル内でncbi2naデータの直後に格納されているため、`.nin` ファイルのオフセット情報を使えば自前でアクセスできるが、APIが提供されていないため実装が複雑になる。

**→ 方法Aを採用する。** `CSeqDBExpert` は `CSeqDB` のサブクラスであり、既存の `CSeqDB` 利用箇所を `CSeqDBExpert` に変更するだけで使用可能。

### 3.2 Ambiguityデータのバイナリレイアウト

Ambiguityデータは、ビッグエンディアンの `Int4`（32ビット整数）の配列として格納される。ホストバイトオーダーへの変換が必要。

#### ヘッダワード（最初の4バイト）

```
ビット31 (0x80000000): フォーマットフラグ
  0 = 旧フォーマット（短い配列用、4バイト/エントリ）
  1 = 新フォーマット（長い配列用、8バイト/エントリ）

ビット30-0: エントリに関する情報
  旧フォーマット: エントリ数
  新フォーマット: エントリ数 × 2（Int4ワード数）
```

ambiguityデータが空（長さ0バイト）の場合、その配列には縮重塩基がない。

#### フォーマット選択基準（BLAST DB作成時）

- **旧フォーマット**: 配列長 ≤ 16,777,215塩基 かつ 全ての縮重ラン長 ≤ 15塩基
- **新フォーマット**: 配列長 > 16,777,215塩基 または いずれかの縮重ラン長 > 15塩基

#### 旧フォーマット（4バイト/エントリ）

```
ビット31-28 (4ビット): ncbi4na値（0-15）
ビット27-24 (4ビット): ラン長-1（0-15、実際のラン長 = 値+1、最大16）
ビット23-0  (24ビット): 位置（0-based塩基オフセット）
```

```
┌────┬────┬────────────────────────┐
│VVVV│LLLL│PPPPPPPP PPPPPPPP PPPPPPPP│
│4bit│4bit│        24 bit           │
└────┴────┴────────────────────────┘
  V = ncbi4na値    L = ラン長-1    P = 位置
```

パース：
```cpp
uint32_t word = ntohl(*(uint32_t*)(ambig_ptr + i*4));  // ビッグエンディアン→ホスト
uint8_t  residue  = (word >> 28) & 0xF;       // ncbi4na値
uint32_t run_len  = ((word >> 24) & 0xF) + 1; // 実ラン長
uint32_t position = word & 0x00FFFFFF;         // 塩基位置
```

#### 新フォーマット（8バイト/エントリ）

ワード0:
```
ビット31-28 (4ビット):  ncbi4na値（0-15）
ビット27-16 (12ビット): ラン長-1（0-4095、実際のラン長 = 値+1、最大4096）
ビット15-0  (16ビット): 未使用（ゼロ）
```

ワード1:
```
ビット31-0 (32ビット): 位置（0-based塩基オフセット、最大4GB）
```

```
ワード0: ┌────┬────────────┬────────────────┐
         │VVVV│LLLLLLLLLLLL│0000000000000000│
         │4bit│   12 bit   │    16 bit      │
         └────┴────────────┴────────────────┘
ワード1: ┌────────────────────────────────┐
         │PPPPPPPP PPPPPPPP PPPPPPPP PPPPPPPP│
         │             32 bit               │
         └────────────────────────────────┘
```

パース：
```cpp
uint32_t w0 = ntohl(*(uint32_t*)(ambig_ptr + offset));
uint32_t w1 = ntohl(*(uint32_t*)(ambig_ptr + offset + 4));
uint8_t  residue  = (w0 >> 28) & 0xF;
uint32_t run_len  = ((w0 >> 16) & 0xFFF) + 1;
uint32_t position = w1;
```

## 4. IUPAC縮重コードとncbi4naの対応

### 4.1 ncbi4naエンコーディング

ncbi4naは4ビットのビットフィールドで、各ビットが1つの塩基を表す：

```
ビット0 (0x1) = A
ビット1 (0x2) = C
ビット2 (0x4) = G
ビット3 (0x8) = T
```

### 4.2 全16コード一覧

| ncbi4na値 | 2進 | IUPAC | 表す塩基 | 展開数 | ncbi2na展開先 |
|:---:|:---:|:---:|:---|:---:|:---|
| 0 | 0000 | - (gap) | ギャップ | 0 | （スキップ） |
| 1 | 0001 | A | A | 1 | {0} |
| 2 | 0010 | C | C | 1 | {1} |
| 3 | 0011 | M | A, C | 2 | {0, 1} |
| 4 | 0100 | G | G | 1 | {2} |
| 5 | 0101 | R | A, G | 2 | {0, 2} |
| 6 | 0110 | S | C, G | 2 | {1, 2} |
| 7 | 0111 | V | A, C, G | 3 | {0, 1, 2} |
| 8 | 1000 | T | T | 1 | {3} |
| 9 | 1001 | W | A, T | 2 | {0, 3} |
| 10 | 1010 | Y | C, T | 2 | {1, 3} |
| 11 | 1011 | H | A, C, T | 3 | {0, 1, 3} |
| 12 | 1100 | K | G, T | 2 | {2, 3} |
| 13 | 1101 | D | A, G, T | 3 | {0, 2, 3} |
| 14 | 1110 | B | C, G, T | 3 | {1, 2, 3} |
| 15 | 1111 | N | A, C, G, T | 4 | {0, 1, 2, 3} |

### 4.3 展開アルゴリズム

ncbi4na値 `v` からncbi2na値の集合を生成：

```cpp
for (uint8_t base2 = 0; base2 < 4; base2++) {
    if (v & (1u << base2)) {
        // base2 は展開先の1つ
    }
}
```

ncbi4naのビット配置がncbi2na値と直接対応するため、変換テーブル不要。

## 5. インデックス仕様との互換性

### 5.1 結論：既存インデックスフォーマットは変更不要

縮重塩基を展開すると、同一OID・同一位置から複数の **異なる** k-mer値が生成される。各k-mer値は異なるポスティングリストに格納されるため、個々のポスティングリスト内で同一OID・同一位置が重複登録されることはない。

例：k=4、位置100に縮重塩基R (A or G) を含むk-mer `AC[R]T` の場合：
- `ACAT` (k-mer値 = 0x0023) → kmer=0x0023 のポスティングリストに (OID, pos=100)
- `ACGT` (k-mer値 = 0x002B) → kmer=0x002B のポスティングリストに (OID, pos=100)

この2つは別々のポスティングリストなので干渉しない。

### 5.2 デルタ圧縮への影響

IDデルタ圧縮（同一OIDはdelta=0で「同一シーケンスの続き」を意味）の動作に影響なし。各k-mer値のポスティングリストは独立しており、リスト内のエントリ順序（OID昇順、同一OID内はpos昇順）も保たれる。

### 5.3 counts・total_postingsへの影響

展開により総ポスティング数が増加するが、Phase 1のカウントパスで正しくカウントすれば、`counts[]` と `total_postings` は自動的に正しい値になる。

### 5.4 検索側への影響

検索側（`ikafssnsearch`）はクエリ配列のk-merスキャンに従来の `KmerScanner` を使用し、クエリにNが含まれるk-merはスキップする。インデックス側で縮重展開して登録したk-merは、クエリ側の通常k-merでヒットするため、検索側の変更は不要。

## 6. 実装計画

### Phase 1: BlastDbReaderの変更

**対象ファイル:**
- `src/io/blastdb_reader.hpp`
- `src/io/blastdb_reader.cpp`

**変更内容:**

1. 内部の `CSeqDB` を `CSeqDBExpert` に変更
2. 新メソッド `get_raw_sequence()` を追加：

```cpp
struct RawSequence {
    const char* ncbi2na_data;  // ncbi2naパックデータへのポインタ（mmap領域）
    int ncbi2na_bytes;         // ncbi2naデータのバイト長
    const char* ambig_data;    // ambiguityデータへのポインタ（ncbi2na_data + ncbi2na_bytes）
    int ambig_bytes;           // ambiguityデータのバイト長
    uint32_t seq_length;       // 塩基数
};

RawSequence get_raw_sequence(uint32_t oid) const;
void ret_raw_sequence(const RawSequence& raw) const;  // RetSequence()ラッパー
```

3. `get_sequence()` の内部実装を `get_raw_sequence()` + ncbi2naデコードに書き換える（後述Phase 7参照）
4. `#include <objtools/blast/seqdb_reader/seqdbexpert.hpp>` を追加

### Phase 2: Ambiguityパーサーの実装

**新規ファイル:**
- `src/core/ambiguity_parser.hpp`

**内容:**

```cpp
struct AmbiguityEntry {
    uint32_t position;   // 0-based塩基位置
    uint32_t run_length; // 連続する縮重塩基の数
    uint8_t  ncbi4na;    // ncbi4na値 (0-15)
};

class AmbiguityParser {
public:
    // ambiguityデータをパースしてエントリ一覧を返す
    // ambig_data: ambiguityデータへのポインタ
    // ambig_bytes: ambiguityデータのバイト長
    static std::vector<AmbiguityEntry> parse(const char* ambig_data, int ambig_bytes);
};
```

パース処理：
1. `ambig_bytes == 0` ならば空ベクタを返す（縮重塩基なし）
2. 先頭ワード（ビッグエンディアン）を読み、ビット31でフォーマット判定
3. 旧フォーマット: 4バイト/エントリでパース
4. 新フォーマット: 8バイト/エントリでパース
5. エントリを位置順にソートして返す（通常は既にソート済みだが保証のため）

### Phase 3: ncbi2na直接k-merスキャナーの実装

**新規ファイル:**
- `src/core/packed_kmer_scanner.hpp`

**設計:**

```cpp
template <typename KmerInt>
class PackedKmerScanner {
public:
    PackedKmerScanner(int k);

    // ncbi2naパックデータから直接k-merをスキャンする。
    // ambig_entries: Phase 2でパースしたambiguityエントリ（位置順）
    // seq_length: 塩基数
    //
    // callback(uint32_t pos, KmerInt kmer) は通常のk-mer（縮重塩基なし）に対して呼ばれる。
    // ambig_callback(uint32_t pos, KmerInt base_kmer, uint8_t ncbi4na, int bit_offset)
    //   は縮重塩基1つを含むk-merに対して呼ばれる。
    //   呼び出し元で展開処理を行う。
    //
    // ※ 2箇所以上の縮重塩基を含むk-merはスキップ。
    template <typename Callback, typename AmbigCallback>
    void scan(const char* ncbi2na_data, uint32_t seq_length,
              const std::vector<AmbiguityEntry>& ambig_entries,
              Callback&& callback,
              AmbigCallback&& ambig_callback) const;
};
```

**スキャンアルゴリズム:**

1. ncbi2naデータから2ビットずつ抽出してスライディングウィンドウを維持：
   ```cpp
   // i番目の塩基の2ビットコードを取得
   uint8_t byte = ncbi2na_data[i >> 2];
   uint8_t code = (byte >> (6 - 2 * (i & 3))) & 0x03;
   kmer = ((kmer << 2) | code) & mask;
   ```

2. ambiguityエントリリストを線形スキャンし、現在のk-merウィンドウ `[pos, pos+k)` 内にある縮重塩基の数を追跡。ambiguityエントリのランが k-merウィンドウと重なる縮重塩基位置をカウント。

3. 縮重塩基数に応じた処理：
   - **0個**: 通常の `callback(pos, kmer)` を呼ぶ
   - **1個**: `ambig_callback` を呼ぶ。呼び出し元（index_builder）で展開
   - **2個以上**: スキップ

4. ambiguityエントリリストの走査は、k-merウィンドウの移動に合わせてイテレータを進めるだけなので O(1) per base（全体で O(n + m)、n=塩基数、m=ambiguityエントリ数）。

**縮重塩基追跡の詳細設計:**

k-merウィンドウ内の縮重塩基を効率的に追跡するため、長さ k のリングバッファ（各位置が縮重かどうかのフラグ）を使用する：

```cpp
// ambiguity_bitmap: 配列全体の各塩基が縮重かどうか
// → 事前にambiguityエントリからビットマップを構築するのはメモリ過大
//
// 代替案: ambiguityエントリリストをソート済みとし、
// 「現在ウィンドウ内にある縮重塩基数」をインクリメンタルに更新
//
// ウィンドウが1塩基右にスライドするとき:
//   - 左端から出る塩基が縮重だったら count--
//   - 右端に入る塩基が縮重かチェックし、そうなら count++
```

具体的には、ambiguityエントリを展開して個々の縮重塩基位置のリスト（ソート済み）を生成し、2つのポインタ（ウィンドウ左端用、右端用）で追跡する。ただし、展開するとメモリ使用量が大きくなる可能性がある。

**最終的な設計方針:** ambiguityエントリは（位置, ラン長）のペアのリストとしてそのまま保持する。ウィンドウ内の縮重塩基数カウントは、エントリリスト上の2つのイテレータ（enter側、leave側）で管理する。各イテレータはランの中のオフセットも追跡する。

### Phase 4: index_builderの変更

**対象ファイル:**
- `src/index/index_builder.cpp`

**変更内容:**

Phase 1（カウントパス）とPhase 2-3（ポスティング書き込み）の両方で、スキャンループを以下に変更：

```cpp
// 旧:
std::string seq = db.get_sequence(oid);
scanner.scan(seq.data(), seq.size(), callback);

// 新:
auto raw = db.get_raw_sequence(oid);
auto ambig = AmbiguityParser::parse(raw.ambig_data, raw.ambig_bytes);
packed_scanner.scan(
    raw.ncbi2na_data, raw.seq_length, ambig,
    // 通常k-merコールバック
    [&](uint32_t pos, KmerInt kmer) {
        // 従来と同じ処理
    },
    // 縮重1塩基k-merコールバック
    [&](uint32_t pos, KmerInt base_kmer, uint8_t ncbi4na, int bit_offset) {
        // base_kmer 内の bit_offset 位置にある2ビットを
        // ncbi4naが表す各塩基に置換して展開
        KmerInt clear_mask = ~(KmerInt(0x03) << bit_offset);
        KmerInt cleared = base_kmer & clear_mask;
        for (uint8_t b = 0; b < 4; b++) {
            if (ncbi4na & (1u << b)) {
                KmerInt expanded = cleared | (KmerInt(b) << bit_offset);
                // 登録処理（countインクリメントまたはbufferにpush_back）
            }
        }
    }
);
db.ret_raw_sequence(raw);
```

`bit_offset` は、k-mer整数値内での縮重塩基の2ビット位置。k-merの最右端（最後に追加された塩基）が bit_offset=0、最左端（最初の塩基）が bit_offset=2*(k-1)。

### Phase 5: CMakeLists.txt の変更

`CSeqDBExpert` を使用するため、必要であればリンク対象ライブラリに `seqdb` が含まれていることを確認する（現在既に含まれている可能性が高い）。`seqdbexpert.hpp` のインクルードパスが通っていることを確認する。

### Phase 6: テスト

**テストデータ:**

`test/testdata/` に、全てのIUPAC縮重コードを含むFASTAファイルを追加：

```fasta
>test_ambig_all_iupac
ACGTACGTMRWSYKVHDBN
ACGTACGTACGTACGTACGT
>test_ambig_single_in_kmer
ACGTACGNACGTACGT
>test_ambig_double_in_kmer
ACNNACGTACGT
>test_long_ambig_run
ACGTNNNNNNNNNNNNNNNNACGT
```

**テストケース（`test/test_packed_kmer_scanner.cpp`）:**

1. ncbi2naバイト列からの正しいk-mer抽出（従来のKmerScannerと結果一致）
2. 縮重塩基なしの配列で従来と同一結果
3. 縮重塩基1個含有k-merの正しい展開（全IUPAC縮重コード）
   - 2塩基縮重（M, R, W, S, Y, K）→ 2展開
   - 3塩基縮重（V, H, D, B）→ 3展開
   - 4塩基縮重（N）→ 4展開
4. 縮重塩基2個以上含有k-merのスキップ確認
5. ラン長 > 1 の縮重領域の正しい処理
6. 配列末尾（端数バイト）の正しい処理
7. AmbiguityParserの旧フォーマット・新フォーマット両対応テスト
8. **統合テスト**: 縮重塩基を含むBLAST DBでインデックスを構築し、展開されたk-merが検索でヒットすることを確認

**ベンチマーク:**

テスト用BLAST DB（ntの一部ボリュームなど）で、旧実装と新実装のインデックス構築時間を比較。

### Phase 7: 既存コードのクリーンアップと削除対象

今回の変更により不要となるコードを以下に整理する。

#### 7a. 削除する関数

| 関数 | ファイル | 削除理由 |
|---|---|---|
| `ncbi_na8_to_char()` | `src/io/blastdb_reader.cpp` | `get_sequence()` の内部実装がncbi2naデコードに変わるため、NA8→char変換は不要になる |

#### 7b. 削除する呼び出し・依存

| 呼び出し | ファイル | 削除理由 |
|---|---|---|
| `GetAmbigSeq()` 呼び出し | `src/io/blastdb_reader.cpp` (`get_sequence()` 内) | `get_raw_sequence()` + ncbi2naデコードに置き換え |
| `RetAmbigSeq()` 呼び出し | `src/io/blastdb_reader.cpp` (`get_sequence()` 内) | 上記に伴い不要 |
| `#include <objtools/blast/seqdb_reader/seqdbcommon.hpp>` | `src/io/blastdb_reader.cpp` | `kSeqDBNuclNcbiNA8` 定数の提供元。`GetAmbigSeq()` 廃止により不要 |

#### 7c. `index_builder.cpp` から除去されるインクルードと使用

| 項目 | 詳細 |
|---|---|
| `#include "core/kmer_encoding.hpp"` | `KmerScanner` を使わなくなるため `index_builder.cpp` からは除去。ヘッダ自体は残る |
| `KmerScanner` のインスタンス化 | `index_builder.cpp` 内の3箇所（Phase 1並列/逐次、Phase 2-3）を `PackedKmerScanner` に置き換え |

#### 7d. 残すもの（他で使用中）

| 項目 | 残存利用元 |
|---|---|
| `BlastDbReader::get_sequence()` | `src/ikafssnretrieve/local_retriever.cpp` で部分配列抽出に使用。内部実装を `get_raw_sequence()` + ncbi2naデコードに書き換えるが、メソッド自体とシグネチャは維持 |
| `KmerScanner` クラス | `src/search/volume_searcher.cpp`（クエリk-mer抽出）、`src/ikafssnserver/request_processor.cpp`、`src/ikafssnsearch/main.cpp`、テスト多数で使用 |
| `encode_base()` / `base_encode_table()` / `BASE_ENCODE_INVALID` | `KmerScanner` 内部で使用。`kmer_encoding.hpp` にそのまま残る |
| `kmer_revcomp()` | 検索側で使用。変更なし |

#### 7e. `get_sequence()` の内部書き換え

`get_sequence()` は `ikafssnretrieve` で引き続き使用されるため、メソッドは維持するが、内部実装を以下のように書き換える：

```cpp
// 旧実装（削除）:
//   GetAmbigSeq() → ncbi_na8_to_char() でASCII文字列構築

// 新実装:
//   get_raw_sequence() でncbi2naパック + ambiguityデータ取得
//   ncbi2naバイト列を4塩基ずつ "ACGT" 文字にデコード
//   AmbiguityParser でambiguityエントリ取得
//   該当位置にIUPAC文字を上書き
```

これにより `ncbi_na8_to_char()` と `GetAmbigSeq()` / `kSeqDBNuclNcbiNA8` への依存を完全に排除する。

ncbi2naからASCII文字へのデコードには以下の4要素テーブルを使用する：

```cpp
static constexpr char ncbi2na_to_char[4] = {'A', 'C', 'G', 'T'};
```

ncbi4naからIUPAC文字へのデコードには以下の16要素テーブルを使用する：

```cpp
static constexpr char ncbi4na_to_iupac[16] = {
    '-',  // 0: gap
    'A',  // 1
    'C',  // 2
    'M',  // 3: A,C
    'G',  // 4
    'R',  // 5: A,G
    'S',  // 6: C,G
    'V',  // 7: A,C,G
    'T',  // 8
    'W',  // 9: A,T
    'Y',  // 10: C,T
    'H',  // 11: A,C,T
    'K',  // 12: G,T
    'D',  // 13: A,G,T
    'B',  // 14: C,G,T
    'N',  // 15: A,C,G,T
};
```

## 7. 注意事項

### 7.1 GetSequence() の戻りポインタのライフタイム

`GetSequence()` / `GetRawSeqAndAmbig()` が返すポインタはmmap領域を指す。`RetSequence()` を呼ぶまで有効。TBBの `parallel_for` 内で複数スレッドが同時に異なるOIDの `GetRawSeqAndAmbig()` を呼ぶ場合、NCBI Toolkit内部のスレッドセーフティを確認する必要がある。

現在の `GetAmbigSeq()` はバッファを `new` で確保して返すため自然にスレッドセーフだが、`GetSequence()` / `GetRawSeqAndAmbig()` はmmap領域を共有するため、内部のアトラス（メモリマッピング管理）のロック機構に依存する。現在のTBB並列処理で `get_sequence()` が動作しているのと同様、`GetRawSeqAndAmbig()` も動作するはずだが、テストで確認する。

### 7.2 エンディアン

ambiguityデータはビッグエンディアンで格納されている。x86_64（リトルエンディアン）では `ntohl()` / `__builtin_bswap32()` でバイトスワップが必要。

### 7.3 max_freq_build との相互作用

縮重展開によりポスティング数が増加するが、展開されたk-mer値が高頻度k-merになる可能性は低い（縮重塩基は稀）。`max_freq_build` 閾値は展開後のカウントに基づいて適用されるため、特別な対処は不要。
