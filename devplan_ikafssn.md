# ikafssn 開発計画書

**Independent programs of K-mer-based Alignment-Free Similarity Search for Nucleotide sequences**

バージョン: 0.1.0 (Draft)
最終更新: 2026-02-20

---

## 目次

1. [プロジェクト概要](#1-プロジェクト概要)
2. [用語定義](#2-用語定義)
3. [設計方針と制約](#3-設計方針と制約)
4. [k-mer エンコーディング仕様](#4-k-mer-エンコーディング仕様)
5. [ファイルフォーマット仕様](#5-ファイルフォーマット仕様)
6. [インデックス構築](#6-インデックス構築)
7. [検索アルゴリズム](#7-検索アルゴリズム)
8. [マルチボリューム対応](#8-マルチボリューム対応)
9. [並列化戦略](#9-並列化戦略)
10. [コマンド体系](#10-コマンド体系)
11. [サーバアーキテクチャ](#11-サーバアーキテクチャ)
12. [クライアント-サーバ間プロトコル](#12-クライアントサーバ間プロトコル)
13. [出力フォーマットとパイプライン](#13-出力フォーマットとパイプライン)
14. [ソースツリー](#14-ソースツリー)
15. [ビルドシステム](#15-ビルドシステム)
16. [サイズ見積もり](#16-サイズ見積もり)
17. [テスト計画](#17-テスト計画)
18. [開発フェーズ](#18-開発フェーズ)
19. [将来の拡張](#19-将来の拡張)

---

## 1. プロジェクト概要

### 1.1 目的

NCBI BLAST DB に格納されたヌクレオチド配列を対象に、完全転置索引 (complete inverted index) を構築し、配列アライメントを行わずに k-mer 類似性のみで高速検索を行うプログラム群を開発する。

### 1.2 主要機能

- BLAST DB から全 k-mer を抽出し、直接アドレス方式の転置索引を独自バイナリ形式で構築する。
- クエリ配列に含まれる k-mer を索引から検索し、対象配列上でクエリの k-mer 並び順に沿って共線的 (collinear) にマッチする領域を特定する。
- マッチ領域内の k-mer マッチ数をスコアとして算出し、上位候補を報告する。
- 巨大データベース (NCBI nt 規模) にも、ボリューム単位の分割インデックスと並列検索で対応する。
- ローカル直接検索に加え、クライアント-サーバアーキテクチャによるデーモン運用と HTTP REST API を提供する。

### 1.3 非目標

- 配列アライメント (Smith-Waterman, BLAST 等) の実行は行わない。
- タンパク質配列への対応は行わない。
- 縮重塩基 (R, Y, M, K, W, S, B, D, H, V) への対応は行わない。N を含む k-mer はスキップする。

---

## 2. 用語定義

| 用語 | 定義 |
|---|---|
| k-mer | 長さ k の連続塩基部分文字列 |
| ボリューム | BLAST DB の分割単位ファイル。接尾辞 `.00`, `.01`, ... で識別 |
| posting list | ある k-mer の全出現箇所 (seq_id, pos) のリスト |
| ID posting | posting list から seq_id のみを抽出した列 |
| pos posting | posting list から pos のみを抽出した列 |
| チェイン / チェイニング | クエリと対象配列の間で、順序を保つ k-mer マッチの部分集合を求める操作 |
| 対角線 (diagonal) | ヒット (q_pos, s_pos) に対して d = s_pos − q_pos で定義される値。同一の挿入欠失なしマッチ領域上のヒットは同一の対角線値を持つ |
| Stage 1 | 二段階検索の第一段。ID posting のみを用いた候補配列の高速絞り込み |
| Stage 2 | 二段階検索の第二段。候補配列に限定した位置情報の取得とチェイニング |
| seqidlist | 検索対象を限定または除外するためのアクセッション番号リスト。テキスト形式 (1行1アクセッション) またはバイナリ形式 (`blastdb_aliastool` で生成) |
| OID フィルタ | seqidlist を OID のビットセットに変換したもの。Stage 1 の posting デコード中に適用する |
| グレースフルシャットダウン | サーバ終了時に実行中リクエストの完了を待ってから終了する方式 |

---

## 3. 設計方針と制約

### 3.1 k-mer 長の制限

- サポート範囲: k = 5 〜 13
- k ≤ 8: k-mer 値を `uint16_t` で表現 (最大 16 bits)
- k = 9 〜 13: k-mer 値を `uint32_t` で表現 (最大 26 bits)
- k ≥ 14 はサポートしない。
- k ≤ 10 は小〜中規模 DB 向けである。NCBI nt 規模 (数百 GB) では k ≥ 11 を推奨する (counts テーブルの uint32 オーバーフロー回避のため)。

### 3.2 インデックス方式

全 k-mer 値を配列インデックスとして直接参照する直接アドレステーブル方式に統一する。k ≤ 13 のとき 4^k ≤ 67,108,864 (64M) であり、テーブル全体 (offsets + counts) は最大でも 768 MB に収まる。ソート済み配列方式やハッシュマップ方式は不要である。

### 3.3 ファイル構成

1つの BLAST DB ボリュームにつき3つのインデックスファイルを生成する。

| ファイル | 拡張子 | 内容 |
|---|---|---|
| メインインデックス | `.kix` | ヘッダ、直接アドレステーブル (offsets, counts)、ID posting (デルタ圧縮) |
| ポジションデータ | `.kpx` | ポジション用オフセットテーブル、pos posting (デルタ圧縮) |
| 配列メタデータ | `.ksx` | 配列長テーブル、accession 文字列テーブル |

ID posting と pos posting をファイルレベルで分離する理由は、二段階検索の Stage 1 で kpx にアクセスする必要をなくし、OS のページキャッシュ効率を最大化するためである。

### 3.4 逆相補鎖の扱い

インデックスには forward strand の k-mer のみを格納する。検索時にクエリ配列の逆相補配列を生成し、forward / reverse complement の両方でインデックスを検索する。結果には strand 情報 (`+` / `-`) を付与する。

### 3.5 縮重塩基と N の扱い

A, C, G, T 以外の塩基 (N を含む) が k-mer ウィンドウ内に存在する場合、その k-mer はスキップする。N の出現後 k 塩基分のウィンドウが無効になることを N カウンタで管理する。

### 3.6 外部依存

| ライブラリ | 用途 | 使用コマンド |
|---|---|---|
| NCBI C++ Toolkit | BLAST DB 読み出し (`CSeqDB`) | ikafssnindex, ikafssnsearch, ikafssnretrieve |
| Intel TBB (oneTBB) | スレッドプール、並列アルゴリズム | ikafssnindex, ikafssnsearch, ikafssnserver |
| Drogon | HTTP REST API フレームワーク | ikafssnhttpd |
| libcurl | HTTP クライアント | ikafssnclient (HTTP モード), ikafssnretrieve (リモート取得モード) |
| (標準ライブラリ) | ファイル I/O、mmap、ソケット通信等 | 全コマンド |

各コマンドは必要な依存のみをリンクする。特に `ikafssnclient` は NCBI C++ Toolkit や TBB に依存せず、軽量バイナリとして配布可能である。`ikafssnretrieve` の libcurl 依存はリモート取得機能を有効にした場合のみ必要であり、ローカル BLAST DB のみを使用する場合は不要である。

### 3.7 対象プラットフォーム

Linux (x86_64, aarch64) を主対象とする。バイトオーダーはリトルエンディアン固定とする。ビルド時に `static_assert` でエンディアンを検証する。

---

## 4. k-mer エンコーディング仕様

### 4.1 塩基の 2-bit 符号化

```
A = 0b00 (0)
C = 0b01 (1)
G = 0b10 (2)
T = 0b11 (3)
```

256 要素のルックアップテーブルにより、`char` → 2-bit 値の変換を行う。A/C/G/T 以外 (N 等) の入力に対してはセンチネル値 (0xFF 等) を返し、呼び出し側で N カウンタの制御に用いる。

### 4.2 k-mer の整数表現

k-mer の先頭塩基を上位ビット、末尾塩基を下位ビットとして、2k ビットの整数にパッキングする。

```
k-mer "ACGT" (k=4):
  A=00, C=01, G=10, T=11
  → 0b00011011 = 0x1B
```

有効ビット数は 2k であり、上位ビットは常にゼロである。マスクは `(1 << 2k) - 1` で算出する。

### 4.3 スライディングウィンドウによる逐次計算

```
next_kmer = ((prev_kmer << 2) | encode(base)) & kmer_mask
```

この操作により、配列走査時に各位置の k-mer を O(1) で算出する。

### 4.4 逆相補変換

2-bit 表現では complement(b) = 3 − b = ~b & 0x3 である。逆相補 k-mer は、各 2-bit 要素を complement した上で、2-bit 単位で配列を反転する。

```
手順:
  1. 全ビットを反転 (~kmer)
  2. 2-bit ペア単位でバイトリバーサル
  3. 上位の未使用ビットをシフトで除去 (>> (W - 2k))
     W = 16 (uint16_t) or 32 (uint32_t)
```

### 4.5 テンプレートによる型の切り替え

`KmerInt` 型パラメータにより `uint16_t` (k ≤ 8) と `uint32_t` (k = 9〜13) を切り替える。エンコーディング、逆相補変換、インデックスアクセスの全てをテンプレートで統一する。

ランタイムでの型ディスパッチは、インデックスファイルのヘッダから `kmer_type` を読み取り、テンプレートインスタンスを選択する。検索のホットループ内では仮想関数呼び出しを行わず、テンプレートインスタンスを直接使用する。

---

## 5. ファイルフォーマット仕様

### 5.1 共通規約

- 全ての多バイト整数はリトルエンディアンで格納する。
- セクション境界は 8 バイト境界にアラインする。アライン用パディングは 0x00 で埋める。
- 可変長整数 (varint) は LEB128 形式 (各バイトの最上位ビットが継続フラグ、下位 7 ビットがデータ) を使用する。

### 5.2 .kix (メインインデックス)

#### ヘッダ (64 バイト固定長)

| オフセット | フィールド | 型 | 説明 |
|---|---|---|---|
| 0x00 | `magic` | char[4] | `"KMIX"` |
| 0x04 | `format_version` | uint16 | 現行 = 1 |
| 0x06 | `k` | uint8 | k-mer 長 (5〜13) |
| 0x07 | `kmer_type` | uint8 | 0 = uint16 (k ≤ 8), 1 = uint32 (k ≥ 9) |
| 0x08 | `num_sequences` | uint32 | このボリューム内の配列数 |
| 0x0C | `total_postings` | uint64 | ID posting エントリ総数 |
| 0x14 | `flags` | uint32 | ビットフラグ (後述) |
| 0x18 | `volume_index` | uint16 | ボリューム番号 (0 起算) |
| 0x1A | `total_volumes` | uint16 | 総ボリューム数 |
| 0x1C | `db_name_len` | uint16 | `db_name` 内の有効文字列長 (バイト) |
| 0x1E | `reserved` | uint8[2] | 予約 (0 埋め) |
| 0x20 | `db_name` | char[32] | 元 BLAST DB 名 (固定 32 バイト、未使用部は 0x00 パディング) |

#### flags ビット定義

| ビット | 名前 | 説明 |
|---|---|---|
| 0 | `SEQ_ID_WIDTH` | 0 = seq_id は uint32, 1 = uint64 (将来用) |
| 1 | `HAS_KSX` | 0 = .ksx なし, 1 = .ksx あり |
| 2〜31 | 予約 | 0 埋め |

#### テーブルセクション

| データ | 型 | 要素数 | 説明 |
|---|---|---|---|
| `offsets` | uint64[] | 4^k | 各 k-mer の ID posting セクション内バイトオフセット |
| `counts` | uint32[] | 4^k | 各 k-mer の出現回数 |

`offsets` は ID posting セクション先頭を基準とするバイトオフセットである。空の posting list を持つ k-mer では、次の有効な k-mer のオフセットと同値になる (長さゼロの領域を表現する)。

`counts` は `offsets` の隣接差分から計算可能であるが、以下の理由で独立して保持する:

- 検索時に posting list を読む前に高頻度 k-mer をスキップ判定できる。
- `counts` 配列と `offsets` 配列のメモリアクセスパターンが異なるため、キャッシュ効率上分離した方が有利。

#### ID posting セクション

各 k-mer の posting list に含まれる seq_id のデルタ符号化列を、k-mer 値の昇順に隙間なく連結して格納する。

デルタ符号化規則:
- 各 k-mer の posting list 内で、seq_id は昇順にソートされている。
- 先頭の seq_id は raw varint で格納する。
- 2番目以降は直前の seq_id との差分 (非負整数) を varint で格納する。
- 同一 seq_id が連続する場合、差分は 0 となる。

### 5.3 .kpx (ポジションデータ)

#### ヘッダ (32 バイト固定長)

| オフセット | フィールド | 型 | 説明 |
|---|---|---|---|
| 0x00 | `magic` | char[4] | `"KMPX"` |
| 0x04 | `format_version` | uint16 | 現行 = 1 |
| 0x06 | `k` | uint8 | k-mer 長 (.kix と一致すること) |
| 0x07 | `reserved1` | uint8 | 予約 |
| 0x08 | `total_postings` | uint64 | (.kix と同数) |
| 0x10 | `reserved2` | uint8[16] | 予約 |

#### ポジション用オフセットテーブル

| データ | 型 | 要素数 | 説明 |
|---|---|---|---|
| `pos_offsets` | uint64[] | 4^k | 各 k-mer の pos posting セクション内バイトオフセット |

ID posting と pos posting はデルタ圧縮後のバイト列長が異なるため、独立したオフセットテーブルが必要である。

#### Position posting セクション

各 k-mer の posting list に含まれる pos のデルタ符号化列を、k-mer 値の昇順に隙間なく連結して格納する。kix の ID posting と 1:1 に対応する。

デルタ符号化規則:
- 同一 seq_id 内の pos は昇順にソートされている。
- seq_id が切り替わる位置 (ID posting 側のデルタ ≠ 0) では、pos を raw varint で格納する (デルタをリセット)。
- 同一 seq_id 内 (ID posting 側のデルタ = 0) では、直前の pos との差分を varint で格納する。

この設計により、ID posting と pos posting を同時にデコードする際、ID posting 側のデルタ値が 0 か否かで pos のデルタリセットを判定できる。

### 5.4 .ksx (配列メタデータ)

#### ヘッダ (32 バイト固定長)

| オフセット | フィールド | 型 | 説明 |
|---|---|---|---|
| 0x00 | `magic` | char[4] | `"KMSX"` |
| 0x04 | `format_version` | uint16 | 現行 = 1 |
| 0x06 | `reserved1` | uint16 | 予約 |
| 0x08 | `num_sequences` | uint32 | 配列数 (.kix と一致すること) |
| 0x0C | `reserved2` | uint8[20] | 予約 |

#### 配列長テーブル

| データ | 型 | 要素数 |
|---|---|---|
| `seq_lengths` | uint32[] | num_sequences |

各配列の長さ (塩基数) を OID 順に格納する。チェイニング時の境界チェック、結果表示時の情報付与等に使用する。

#### Accession テーブル

| データ | 型 | 要素数 |
|---|---|---|
| `accession_offsets` | uint32[] | num_sequences + 1 |
| `accession_strings` | char[] | 可変長 |

`accession_offsets[i]` は `accession_strings` 内の配列 i の accession 文字列開始位置を示す。文字列長は `accession_offsets[i+1] - accession_offsets[i]` で算出する。NUL 終端は含まない。

この設計により、検索時に元の BLAST DB を開かずとも accession を取得でき、インデックスファイルのみで完結した検索結果表示が可能になる。

---

## 6. インデックス構築

### 6.1 全体フロー

```
入力: BLAST DB (1ボリュームまたは複数ボリューム)
出力: 各ボリュームにつき .kix, .kpx, .ksx の3ファイル

各ボリュームに対して:
  Phase 0: メタデータ収集 → .ksx 書き出し
  Phase 1: 計数パス      → counts[4^k] 集計
  Phase 2: オフセット計算 → prefix sum
  Phase 3: Posting 書き出し (パーティション × バッファリング)
  Phase 4: ファイナライズ → ヘッダ書き込み、ファイルクローズ
```

### 6.2 Phase 0: メタデータ収集

```
CSeqDB でボリュームを開く
num_sequences を取得
for each OID i in [0, num_sequences):
    seq_lengths[i] = CSeqDB::GetSeqLength(i)
    accession[i] = CSeqDB::GetSeqIDs(i) から主 accession を取得
.ksx に書き出し
```

### 6.3 Phase 1: 計数パス

```
counts[4^k] をゼロ初期化

for each OID i in [0, num_sequences):
    seq = CSeqDB::GetSequence(i)  // 生の塩基配列を取得
    n_count = k - 1  // 最初の k-mer が完成するまでスキップ
    kmer = 0

    for each base in seq:
        if base ∉ {A, C, G, T}:
            n_count = k - 1
            continue

        kmer = ((kmer << 2) | encode(base)) & kmer_mask

        if n_count > 0:
            n_count--
            continue

        counts[kmer]++
```

この段階で `counts` を調査し、構築時高頻度 k-mer 除外オプション (`-max_freq_build`) が指定されている場合は、閾値を超える k-mer の `counts` をゼロにリセットする。これにより低複雑度配列由来の超高頻度 k-mer を索引から除外できる。

なお、Phase 1 の計数ループでは `uint64` アキュムレータを使用し、計数完了後に `uint32` (`counts` テーブルの型) への変換時にオーバーフローを検出する。オーバーフローが発生した場合はエラーを報告し、より大きな k 値の使用を推奨する。

### 6.4 Phase 2: オフセット計算

ID posting と pos posting のオフセットは、デルタ圧縮後のバイト列長に基づいて計算する必要がある。しかし圧縮後のサイズは実際にデルタ符号化しなければ正確にはわからない。

実用的な方法として、2通りのアプローチがある。

**方式 A: 2パス計数 (正確なサイズ計算)**

Phase 1 で counts のみを集計した後、Phase 3 の書き出し時に各パーティション内で一旦メモリ上にデルタ圧縮列を構築し、そのバイト長を確定してからファイルに書き出す。オフセットはパーティションの書き出し進行に伴い逐次確定させる。

**方式 B: 事後オフセット構築**

Phase 3 で各 k-mer の圧縮データを書き出す際に、書き出し先のバイト位置を記録してオフセットテーブルを構築する。全 posting 書き出し完了後にオフセットテーブルをファイル先頭付近に書き戻す (seek + write)。

方式 B の方が実装が単純であり、推奨する。ファイルのテーブルセクション領域をあらかじめ確保 (ゼロ埋め) しておき、posting 書き出し後にテーブルを上書きする。

### 6.5 Phase 3: Posting 書き出し

パーティション方式とバッファリング方式を組み合わせる。

#### パーティション分割

k-mer 値の上位ビット (またはハッシュの上位ビット) により、全 k-mer 空間を P 個のパーティションに分割する。

```
partition_bits = ceil(log2(P))
partition_of(kmer) = (kmer >> (2k - partition_bits)) & (P - 1)

あるいはハッシュベース (k-mer 頻度分布の偏りを軽減):
partition_of(kmer) = (hash(kmer) >> (32 - partition_bits)) & (P - 1)
```

ハッシュベースを採用する場合、構築と検索の両方で同一のハッシュ関数を使用する必要がある。上位ビット方式はハッシュ関数が不要で単純であり、k-mer 空間が小さい (4^k ≤ 64M) ため偏りの実害も限定的である。初期実装では上位ビット方式を推奨する。

#### 各パーティション内のバッファリング

```
buffer_entries = buffer_size / sizeof(TempEntry)
  TempEntry = { kmer_value: KmerInt, seq_id: uint32, pos: uint32 }
  sizeof(TempEntry) = 12 (全 k 値共通。KmerInt = uint16_t 時もアラインメントにより 12 バイト)

パーティション p を処理:
  buffer を確保 (最大 buffer_entries 個)
  
  for each OID i in [0, num_sequences):
      seq = get_sequence(i)
      for each (kmer, pos) in sliding_window(seq, k):
          if partition_of(kmer) ≠ p: continue
          if counts[kmer] == 0: continue  // 構築時除外済み
          buffer に (kmer, seq_id=i, pos) を追加
          if buffer が満杯:
              flush_buffer(buffer)

  残りの buffer を flush

flush_buffer:
  buffer を kmer 値でソート (同一 kmer 内は seq_id, pos 昇順)
  各 kmer の posting をデルタ圧縮し、kix (seq_id) と kpx (pos) に追記
  書き出しバイト位置を記録 → オフセットテーブル更新
```

#### マージ書き出し

同一 k-mer に対して複数回の flush が発生する場合、各 flush の posting list はそれぞれソート済みだが、flush 間のマージが必要になる。

2つの戦略がある:

**戦略 1: パーティション内全データをメモリに蓄積し一括書き出し**

パーティション内の全 posting がバッファに収まる場合 (小規模 DB やパーティション数が多い場合)、ソート後に一括してデルタ圧縮・書き出しを行う。マージは不要。

**戦略 2: 一時ファイルへの中間フラッシュとマージ**

バッファが満杯になるたびに一時ファイルにソート済みブロックを書き出し、パーティション完了後にマルチウェイマージを行って最終的なデルタ圧縮 posting を生成する。

初期実装では戦略 1 を前提とし、パーティション数を十分に増やすことでパーティション内データがバッファに収まるようにする。巨大 DB でバッファに収まらない場合のために、戦略 2 を将来の最適化として残す。

**推奨されるデフォルト動作:**

```
パーティション内の見積もり posting 数 = total_postings / P
必要メモリ ≈ (total_postings / P) × sizeof(TempEntry)

buffer_size ≥ 必要メモリ なら戦略 1 (一括書き出し)
buffer_size < 必要メモリ なら:
  ユーザーに P を増やすか buffer_size を増やすよう警告
  (将来的には戦略 2 にフォールバック)
```

### 6.6 Phase 4: ファイナライズ

- kix のテーブルセクション (offsets, counts) を書き戻す。
- kpx のオフセットテーブル (pos_offsets) を書き戻す。
- 各ファイルのヘッダを書き込む。
- ファイルをクローズし、整合性チェック (total_postings の一致等) を行う。
- 構築中断時 (ディスク容量不足、プロセス強制終了等) に不完全なファイルが残ることを防ぐため、書き出し中はファイル名に `.tmp` 接尾辞を付与し、ファイナライズ完了後に最終名へ `rename` する。

### 6.7 N カウンタの詳細動作

配列走査中に N (および A/C/G/T 以外の塩基) が出現した場合の制御を以下に示す。

```
n_count: 0 になるまでの残りスキップ数。初期値 = k - 1 (k 塩基読み込み後に最初の k-mer が成立)

base を読むたびに:
  if base が無効:
      n_count = k - 1   ← 次の有効な k-mer まで k 塩基必要
  else:
      kmer を更新
      if n_count > 0:
          n_count--     ← まだウィンドウ内に N が残っている
      else:
          この kmer は有効 → counts++ or posting 書き出し
```

このロジックはインデックス構築と検索 (クエリ走査) の両方で共通して使用する。

---

## 7. 検索アルゴリズム

### 7.1 二段階検索の概要

```
入力: クエリ配列 Q, 逆相補 Q', インデックス (kix, kpx, ksx)

各ストランド (Q, Q') について:
  Stage 1: ID posting のみを用い、各 seq_id の k-mer ヒット数を集計
            → 上位候補 seq_id を選出
  Stage 2: 候補 seq_id に限定して pos posting を取得し、チェイニング
            → 最終スコアを算出

両ストランドの結果を統合し、スコア順にソートして出力
```

### 7.2 Stage 1: 候補配列の高速絞り込み

```
score_per_seq[num_sequences] を 0 初期化 (uint32)
oid_filter = OidFilter(seqidlist, negative)  // seqidlist 未指定時は全通過

for each (q_pos, kmer) in query_kmers(Q):
    if counts[kmer] > max_freq: skip   // 高頻度 k-mer を除外
    
    decoder = SeqIdDecoder(kix_id_posting[offsets[kmer]])
    for i in [0, counts[kmer]):
        seq_id = decoder.next()
        if not oid_filter.pass(seq_id): continue   // OID フィルタ
        score_per_seq[seq_id]++

候補選出:
  score_per_seq を走査し、min_stage1_score 以上の seq_id を収集
  スコア上位 stage1_topn 件に絞り込む → 候補集合 C
```

**OID フィルタの実装:**

`-seqidlist` または `-negative_seqidlist` が指定された場合、検索開始前に以下の手順で OID ビットセットを構築する。

```
1. seqidlist ファイルを読み込み (テキストまたはバイナリ形式を自動判別)
2. 各ボリュームの .ksx から accession → OID の逆引きマップを構築
3. seqidlist 内のアクセッションを OID に解決
4. OID ビットセット (std::vector<bool> または bitset) を構築
   - -seqidlist: ビットが立っている OID のみ通過
   - -negative_seqidlist: ビットが立っている OID を除外
5. 解決できなかったアクセッション (ボリュームに存在しない等) は stderr に警告
```

OID ビットセットはボリュームあたり num_sequences ビット (数百万配列で数百 KB) であり、L2 キャッシュに収まる。Stage 1 のデコードループ内でのビット参照は分岐予測が効きやすく、フィルタ適用のオーバーヘッドは小さい。

フィルタが未指定の場合、`oid_filter.pass()` は常に `true` を返すインライン関数とし、コンパイラによる分岐除去を期待する。

**Stage 1 の計算コスト:**

kpx を一切読まない。kix の ID posting はデルタ圧縮で小さく (平均 ~1.2 bytes/entry)、シーケンシャルデコードのみで構成される。score_per_seq 配列は uint32 × num_sequences で、1ボリュームあたり数 MB〜十数 MB に収まる。

### 7.3 Stage 2: チェイニング

#### 7.3.1 候補配列の posting 取得

```
candidate_set = Stage 1 で得た候補 seq_id のハッシュセット

hits_per_seq: map<seq_id, vector<Hit>>
  Hit = { q_pos: uint32, s_pos: uint32 }

for each (q_pos, kmer) in query_kmers(Q):
    if counts[kmer] > max_freq: skip

    id_decoder = SeqIdDecoder(kix_id_posting + offsets[kmer])
    pos_decoder = PosDecoder(kpx_pos_posting + pos_offsets[kmer])
    
    for i in [0, counts[kmer]):
        seq_id = id_decoder.next()
        pos = pos_decoder.next(id_decoder.was_new_seq())  // seq境界でデルタリセット
        
        if seq_id ∈ candidate_set:
            hits_per_seq[seq_id].append(Hit{q_pos, pos})
```

注: `id_decoder.was_new_seq()` は直前の `next()` でデルタが非ゼロだったか (= seq_id が切り替わったか) を返す。pos_decoder はこの情報を用いてデルタリセットを判定する。

#### 7.3.2 対角線フィルタ

```
for each (seq_id, hits) in hits_per_seq:
    diag_counts: map<int32, uint32>
    
    for each hit in hits:
        diag = hit.s_pos - hit.q_pos   // int32 (負値あり)
        diag_counts[diag]++
    
    filtered_hits = []
    for each hit in hits:
        diag = hit.s_pos - hit.q_pos
        if diag_counts[diag] >= min_diag_hits:
            filtered_hits.append(hit)
    
    hits = filtered_hits
```

対角線フィルタにより、散発的なランダムヒットを除去し、チェイニングの入力サイズを削減する。`min_diag_hits` のデフォルトは 2〜3 程度が妥当である。

#### 7.3.3 チェイニング DP

```
hits を q_pos 昇順にソート (同一 q_pos 内は s_pos 昇順)

n = len(hits)
dp[n]: 各ヒットを終端とする最長チェインスコア
prev[n]: トレースバック用ポインタ

for i in [0, n):
    dp[i] = 1        // 自身のみのチェイン
    prev[i] = -1

    for j in [0, i):  // 制約: q_pos[j] < q_pos[i], s_pos[j] < s_pos[i]
        gap_q = q_pos[i] - q_pos[j]
        gap_s = s_pos[i] - s_pos[j]
        diag_diff = |gap_s - gap_q|    // 対角線のずれ (挿入欠失の許容幅)
        
        if s_pos[j] < s_pos[i] and diag_diff ≤ max_gap:
            if dp[j] + 1 > dp[i]:
                dp[i] = dp[j] + 1
                prev[i] = j

best = argmax(dp)
チェインを prev ポインタでトレースバックして復元
```

**計算量について:**

対角線フィルタ後のヒット数 n は通常小さい (数十〜数百)。O(n²) の DP で十分高速である。万一 n が大きい場合 (数千以上) は、Binary Indexed Tree (BIT) を用いた O(n log n) の最適化が可能だが、初期実装では O(n²) で開始する。

#### 7.3.4 スコア算出と結果構成

```
チェインから:
  chain_score = チェイン内の k-mer マッチ数 (= dp[best])
  q_start     = チェイン先頭ヒットの q_pos
  q_end       = チェイン末尾ヒットの q_pos + k
  s_start     = チェイン先頭ヒットの s_pos
  s_end       = チェイン末尾ヒットの s_pos + k
  strand      = "+" or "-" (Q or Q' のどちらで検索したか)

結果レコード:
  (query_id, seq_id, accession, strand, q_start, q_end, s_start, s_end, chain_score)
```

### 7.4 高頻度 k-mer スキップ閾値の自動調整

`-max_freq` が明示指定されない場合、以下の自動調整を行う。

```
mean_count = total_postings / 4^k
default_max_freq = mean_count × 10

ただし下限を 1000、上限を 100000 とする。
```

この値はボリュームごとに異なりうるため、各ボリュームの kix ヘッダから `total_postings` を読み取って算出する。

---

## 8. マルチボリューム対応

### 8.1 BLAST DB のボリューム構造

大規模な BLAST DB は自動的に複数ボリュームに分割される。各ボリュームは独立したファイルセットを持つ。

```
例: nt データベース
  nt.00.nsq, nt.00.nin, nt.00.nhr, ...
  nt.01.nsq, nt.01.nin, nt.01.nhr, ...
  ...
  nt.nal  (エイリアスファイル: 全ボリュームを束ねる)
```

### 8.2 ボリューム検出

構築時: NCBI C++ Toolkit の `CSeqDB::FindVolumePaths()` を使用して、指定された DB プレフィックスに属する全ボリュームのパスを列挙する。

検索時: 指定されたインデックスディレクトリを走査し、`.kix` ファイルをボリューム番号順に列挙する。ファイル名規約は `<db_prefix>.<volume_index>.<kk>mer.kix` とする (`<kk>` は k 値の 2 桁ゼロパディング、例: `09mer`, `13mer`)。`.kpx`, `.ksx` も同一の命名規則に従う。

### 8.3 ボリューム単位の構築

各ボリュームを独立して構築する。ボリューム間の依存関係はない。複数ボリュームの構築を TBB で並列化することも可能だが、I/O 帯域がボトルネックになりやすいため、初期実装ではボリュームを逐次処理する。

### 8.4 ボリューム間の seq_id

各ボリュームの seq_id (OID) はそのボリューム内でのローカル番号 (0 起算) である。検索結果の出力時には、ボリューム番号と OID の組、または .ksx 内の accession で配列を一意に識別する。

グローバル OID (全ボリューム通算の OID) が必要な場合は、各ボリュームの num_sequences の累積和からオフセットを算出する。

---

## 9. 並列化戦略

### 9.1 TBB の使用方針

Intel oneAPI Threading Building Blocks (oneTBB) をスレッドプールおよび並列アルゴリズムの基盤として使用する。

使用する主要コンポーネント:
- `tbb::task_arena`: スレッド数の制御
- `tbb::parallel_for`: インデックス構築の計数パス並列化
- `tbb::parallel_pipeline`: 検索パイプラインのステージ並列化
- `tbb::concurrent_queue` / `tbb::concurrent_vector`: 結果収集

### 9.2 構築時の並列化

```
Phase 1 (計数):
  配列を TBB のチャンクに分割し、スレッドローカルな counts 配列に集計。
  最後に全スレッドの counts をリダクションで統合する。

  tbb::parallel_for で OID 範囲をチャンク分割。
  tbb::combinable<CountsArray> でスレッドローカル counts を管理。

Phase 3 (Posting 書き出し):
  パーティション間の並列化は行わない (I/O 競合を避けるため)。
  パーティション内の DB 走査は Phase 1 と同様に並列化可能だが、
  バッファへの書き込みが競合するため、初期実装ではシングルスレッドとする。
```

### 9.3 検索時の並列化

```
並列化単位: (query, volume) のペア

queries × volumes のジョブ行列を構成し、TBB で並列実行する。

tbb::task_arena arena(num_threads);
arena.execute([&] {
    tbb::parallel_for_each(jobs, [&](const Job& job) {
        // 各ジョブ:
        //   1. volume の kix を mmap (初回のみ、共有)
        //   2. Stage 1 (score_per_seq はスレッドローカル)
        //   3. kpx を mmap (必要に応じて)
        //   4. Stage 2 (チェイニング)
        //   5. 結果をスレッドセーフなキューに投入
    });
});

結果マージ:
  全ジョブ完了後、クエリごとにボリューム横断で結果を統合し、
  スコア順にソートして上位 N 件を出力。
```

### 9.4 mmap の共有

同一ボリュームの kix/kpx は読み取り専用で mmap するため、複数スレッドから安全にアクセスできる。ボリュームごとに1回だけ mmap を行い、全スレッドで共有する。

```
volume_maps: vector<VolumeMap>  // メインスレッドで事前に mmap

struct VolumeMap {
    MmapFile kix;    // .kix 全体
    MmapFile kpx;    // .kpx 全体 (遅延 mmap 可)
    KsxData  ksx;    // .ksx はメモリにロード
};
```

kpx の遅延 mmap: Stage 1 の結果、あるボリュームに候補が存在しない場合、そのボリュームの kpx は mmap しない。これにより不要な仮想アドレス空間の消費を回避する。複数スレッドが同一ボリュームの kpx を同時に要求する場合に備え、`std::once_flag` + `std::call_once` により mmap の実行を排他制御する。

### 9.5 スレッドローカルリソース

Stage 1 の `score_per_seq` 配列はスレッドローカルに確保する。num_sequences × sizeof(uint32) であり、1ボリュームあたり数 MB〜十数 MB 程度のため、スレッド数分確保しても問題ない。

---

## 10. コマンド体系

### 10.1 概要

ikafssn は機能ごとに独立したコマンドとして提供する。サブコマンド方式 (`ikafssn build` 等) は採用しない。理由は以下の通り。

- コマンドごとにリンク依存が大きく異なる。`ikafssnsearch` は NCBI C++ Toolkit と TBB を要するが、`ikafssnclient` はどちらも不要である。単一バイナリに統合するとクライアント専用マシンにも重い依存が強制される。
- デーモン (`ikafssnserver`, `ikafssnhttpd`) はサブコマンドとしての起動が不自然である。
- 各コマンドが独立バイナリであれば、必要なコマンドのみを配布・デプロイできる。

### 10.2 コマンド一覧

| コマンド | 概要 | 主要依存 |
|---|---|---|
| `ikafssnindex` | インデックス構築 | NCBI C++ Toolkit, TBB |
| `ikafssnsearch` | ローカル直接検索 | NCBI C++ Toolkit, TBB |
| `ikafssnretrieve` | マッチ部分配列の抽出 (ローカル BLAST DB または NCBI efetch) | NCBI C++ Toolkit, libcurl (リモート時) |
| `ikafssnserver` | 検索サーバデーモン (UNIX / TCP ソケット) | NCBI C++ Toolkit, TBB |
| `ikafssnhttpd` | HTTP REST API デーモン (`ikafssnserver` のフロントエンド) | Drogon |
| `ikafssnclient` | クライアント (ソケット直接 / HTTP 経由) | libcurl (HTTP モード時) |
| `ikafssninfo` | インデックス / BLAST DB 情報表示 | (NCBI C++ Toolkit: BLAST DB 情報表示時のみ) |

### 10.3 `ikafssnindex`

インデックス構築コマンド。BLAST DB を入力とし、ボリュームごとに `.kix`, `.kpx`, `.ksx` の3ファイルを生成する。

```
ikafssnindex [options]

必須:
  -db <path>              BLAST DB プレフィックス
  -k <int>                k-mer 長 (5〜13)
  -o <dir>                出力ディレクトリ

オプション:
  -buffer_size <size>     バッファサイズ (default: 8G)
                          接尾辞: K, M, G を認識
  -partitions <int>       パーティション数 (default: 4)
                          2の冪乗を推奨 (1, 2, 4, 8, 16, ...)
  -max_freq_build <int>   構築時高頻度 k-mer 除外閾値
                          (default: 0 = 除外なし)
  -threads <int>          計数パスのスレッド数 (default: 1)
  -v, --verbose           詳細ログ出力
```

**使用例:**

```bash
# 小規模 DB, メモリ豊富
ikafssnindex -db mydb -k 11 -o ./index -buffer_size 16G -partitions 1

# 大規模 DB, メモリ限定
ikafssnindex -db nt -k 11 -o ./nt_index -buffer_size 4G -partitions 16

# 高頻度 k-mer を除外して構築
ikafssnindex -db nt -k 11 -o ./nt_index -max_freq_build 50000
```

### 10.4 `ikafssnsearch`

ローカル直接検索コマンド。インデックスファイルを直接 mmap して検索する。サーバプロセスを必要としない独立コマンドである。

```
ikafssnsearch [options]

必須:
  -ix <dir>               インデックスディレクトリ
  -query <path>           クエリ FASTA ファイル (- で標準入力)

オプション:
  -o <path>               出力ファイル (default: 標準出力)
  -threads <int>          並列検索スレッド数 (default: 利用可能な全コア)
  -min_score <int>        最小チェインスコア (default: 3)
  -max_gap <int>          チェイニング対角線ずれ許容幅 (default: 100)
  -max_freq <int>         高頻度 k-mer スキップ閾値
                          (default: 自動計算)
  -min_diag_hits <int>    対角線フィルタ最小ヒット数 (default: 2)
  -stage1_topn <int>      Stage 1 候補数上限 (default: 500)
  -min_stage1_score <int> Stage 1 最小スコア閾値 (default: 2)
  -num_results <int>      最終出力件数 (default: 50)
  -seqidlist <path>       検索対象を指定アクセッションに限定
  -negative_seqidlist <path>  指定アクセッションを検索対象から除外
  -outfmt <tab|json>      出力形式 (default: tab)
  -v, --verbose           詳細ログ出力
```

`-seqidlist` と `-negative_seqidlist` は排他的である (同時指定不可)。ファイル形式はテキスト (1行1アクセッション) と `blastdb_aliastool -seqid_file_in` で生成されるバイナリ形式の両方を受け付ける。ファイル先頭のマジックバイトで自動判別する。

**使用例:**

```bash
# 基本的な検索
ikafssnsearch -ix ./index -query query.fasta -threads 8

# 感度を上げた検索
ikafssnsearch -ix ./index -query query.fasta \
    -min_score 2 -stage1_topn 2000 -max_freq 50000

# seqidlist で検索対象を限定
ikafssnsearch -ix ./index -query query.fasta -seqidlist targets.txt

# negative_seqidlist で特定配列を除外
ikafssnsearch -ix ./index -query query.fasta -negative_seqidlist exclude.txt

# パイプラインで ikafssnretrieve に接続
ikafssnsearch -ix ./index -query query.fasta | ikafssnretrieve -db nt > matches.fasta
```

### 10.5 `ikafssnretrieve`

検索結果に基づき、マッチした部分配列を抽出するコマンド。配列ソースとしてローカル BLAST DB または NCBI Entrez E-utilities (efetch) を選択できる。入力はファイルまたは標準入力で受け付ける。`ikafssnsearch` と `ikafssnclient` の出力形式は同一であるため、どちらの出力にも対応する。

```
ikafssnretrieve [options]

配列ソース (いずれか必須):
  -db <path>              ローカル BLAST DB プレフィックス
  -remote                 NCBI efetch からリモート取得

入力 (いずれか):
  -results <path>         検索結果ファイル (tab 形式)
  (指定なし)              標準入力から読み込み

共通オプション:
  -o <path>               出力 FASTA ファイル (default: 標準出力)
  -context <int>          マッチ領域の前後に付加する塩基数 (default: 0)
  -outfmt <fasta|tab>     出力形式 (default: fasta)
  -v, --verbose           詳細ログ出力

リモート取得オプション (-remote 時):
  -api_key <key>          NCBI API key (環境変数 NCBI_API_KEY でも可)
  -batch_size <int>       バッチあたりの accession 数 (default: 100)
  -retries <int>          リトライ回数 (default: 3)
  -timeout <int>          リクエストタイムアウト秒数 (default: 30)
  -range_threshold <int>  部分取得に切り替える配列長閾値 (default: 100000)
```

**使用例:**

```bash
# ローカル BLAST DB から抽出 (ファイル入力)
ikafssnsearch -ix ./index -query query.fasta -o results.tsv
ikafssnretrieve -db nt -results results.tsv -o matches.fasta

# ローカル BLAST DB から抽出 (パイプ)
ikafssnsearch -ix ./index -query query.fasta | ikafssnretrieve -db nt > matches.fasta

# サーバ経由の検索結果からも同様に利用可能
ikafssnclient -socket /var/run/ikafssn.sock -query query.fasta | ikafssnretrieve -db nt > matches.fasta

# NCBI efetch からリモート取得
ikafssnsearch -ix ./index -query query.fasta | ikafssnretrieve -remote > matches.fasta

# NCBI efetch + API key (高スループット)
ikafssnclient -http http://search.example.com:8080 -query query.fasta \
    | ikafssnretrieve -remote -api_key XXXXXXXX > matches.fasta

# マッチ領域の前後 50bp を含めて抽出
ikafssnretrieve -db nt -results results.tsv -context 50
```

`-context` オプションにより、マッチ領域の前後に指定塩基数を追加して抽出できる。抽出範囲が配列境界を超える場合は配列端で切り詰める。

#### 10.5.1 NCBI efetch リモート取得の仕様

**API エンドポイント:**

```
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi
  ?db=nuccore
  &id=<accession_list>
  &rettype=fasta
  &retmode=text
  &api_key=<key>
```

部分配列取得時は `&seq_start=<start>&seq_stop=<stop>` を付加する (1-based, inclusive)。

**レート制限:**

| 条件 | 制限 | リクエスト間スリープ |
|---|---|---|
| API key なし | 3 リクエスト/秒 | 334 ms |
| API key あり | 10 リクエスト/秒 | 100 ms |

API key が設定されていない場合、起動時に警告メッセージを stderr に出力する。

**バッチ取得戦略:**

検索結果の accession を集約し、efetch のバッチリクエスト (1回の HTTP リクエストで複数 accession を取得) を行う。具体的な手順は以下の通り。

```
1. 検索結果を読み込み、全ヒットの (accession, s_start, s_end) を収集
2. 配列長に基づく取得方式の振り分け:
   a. seq_length ≤ range_threshold (default: 100kb):
      → accession を batch_size ごとにグループ化
      → バッチ efetch で配列全体を取得
      → ローカルで部分配列を切り出し
   b. seq_length > range_threshold:
      → ヒットごとに個別の efetch (seq_start, seq_stop 指定)
      → 染色体レベルの長大配列から小領域を取得する際の帯域節約
3. レスポンス FASTA をパースし、accession でヒットとマッチング
4. 一時的な HTTP エラー (429, 503) は指数バックオフ付きリトライ
```

配列長 (`seq_length`) は検索結果の `s_end` 値から概算するか、`.ksx` の `seq_lengths` テーブルを参照する。`-remote` モードでは `.ksx` にアクセスできない場合があるため、`s_end` 値からの概算をフォールバックとして使用する。

**レスポンスパース:**

バッチ efetch のレスポンスは複数配列の FASTA が連結されて返る。accession の返却順序はリクエスト順と一致しないことがあるため、レスポンス内の defline (`>` 行) から accession を抽出し、リクエスト側のヒットリストとマッチングする。

**エラーハンドリング:**

| エラー | 対処 |
|---|---|
| HTTP 429 (Too Many Requests) | 指数バックオフ付きリトライ (初回 1 秒、最大 `-retries` 回) |
| HTTP 503 (Service Unavailable) | 同上 |
| HTTP 400 (Bad Request) | accession 不正。該当ヒットをスキップし、stderr に警告 |
| HTTP 404 / accession 未発見 | withdrawn 等。該当ヒットをスキップし、stderr に警告 |
| 接続タイムアウト | `-timeout` 秒後にタイムアウト。リトライ |
| レスポンス内の accession 不一致 | リクエストした accession とレスポンスの accession が一致しない場合、stderr に警告 |

### 10.6 `ikafssnserver`

検索サーバデーモン。インデックスをメモリに常駐 (mmap) させ、UNIX ドメインソケットまたは TCP ソケットでクライアントからの検索リクエストを受け付ける。

```
ikafssnserver [options]

必須:
  -ix <dir>               インデックスディレクトリ

リスニング (いずれかまたは両方):
  -socket <path>          UNIX ドメインソケットパス
  -tcp <host>:<port>      TCP リスニングアドレス (例: 0.0.0.0:9100)

オプション:
  -threads <int>          検索ワーカースレッド数 (default: 利用可能な全コア)
  -max_connections <int>  最大同時接続数 (default: 64)
  -shutdown_timeout <int> グレースフルシャットダウンのタイムアウト秒数 (default: 180)
  -pid <path>             PID ファイルパス
  -max_freq <int>         高頻度 k-mer スキップ閾値 (default: 自動計算)
  -min_diag_hits <int>    対角線フィルタ最小ヒット数 (default: 2)
  -stage1_topn <int>      Stage 1 候補数上限 (default: 500)
  -min_stage1_score <int> Stage 1 最小スコア閾値 (default: 2)
  -v, --verbose           詳細ログ出力
```

**使用例:**

```bash
# UNIX ソケットで起動
ikafssnserver -ix ./nt_index -socket /var/run/ikafssn.sock -threads 16

# TCP で起動 (リモートアクセス用)
ikafssnserver -ix ./nt_index -tcp 0.0.0.0:9100 -threads 32

# 両方で同時リスニング
ikafssnserver -ix ./nt_index -socket /var/run/ikafssn.sock -tcp 0.0.0.0:9100

# デーモン運用 (systemd 等から起動)
ikafssnserver -ix ./nt_index -socket /var/run/ikafssn.sock -pid /var/run/ikafssn.pid
```

**運用上の特性:**

- 1プロセスにつき1つの BLAST DB のインデックスのみをサーブする。複数 DB を同時にサーブする場合は、DB ごとに別プロセスを起動する。
- 同一 DB の異なる k-mer サイズのインデックスが `-ix` ディレクトリに存在する場合、全て読み込み、クライアントのリクエストで k を指定できる。
- インデックス更新時はサーバプロセスを再起動する。SIGTERM 受信時はグレースフルシャットダウンを行い、実行中のリクエストの完了を最大 `-shutdown_timeout` 秒 (デフォルト 180 秒) 待つ。待機中は新規リクエストを受け付けない。タイムアウト時は残存リクエストを破棄して終了する。
- リクエストキューはプロセス終了時にフラッシュされるため、再起動後に旧インデックスのリクエストが残存することはない。

### 10.7 `ikafssnhttpd`

HTTP REST API デーモン。`ikafssnserver` に接続し、HTTP REST API として検索サービスを提供する。Drogon フレームワークを使用する。

```
ikafssnhttpd [options]

必須 (いずれか):
  -server_socket <path>           接続先 ikafssnserver の UNIX ソケットパス
  -server_tcp <host>:<port>       接続先 ikafssnserver の TCP アドレス

オプション:
  -listen <host>:<port>   HTTP リスニングアドレス (default: 0.0.0.0:8080)
  -path_prefix <prefix>   API パスプレフィックス (default: なし。例: /nt)
  -threads <int>          Drogon I/O スレッド数 (default: 4)
  -pid <path>             PID ファイルパス
  -v, --verbose           詳細ログ出力
```

**使用例:**

```bash
# ローカルの ikafssnserver に UNIX ソケットで接続
ikafssnhttpd -server_socket /var/run/ikafssn.sock -listen 0.0.0.0:8080

# リモートの ikafssnserver に TCP で接続 (別マシン構成)
ikafssnhttpd -server_tcp 10.0.1.5:9100 -listen 0.0.0.0:8080
```

**設計上の位置づけ:**

- `ikafssnhttpd` は `ikafssnserver` のプロトコル変換プロキシである。HTTP リクエストを受信し、自前バイナリプロトコルに変換して `ikafssnserver` に転送し、レスポンスを HTTP レスポンスに変換して返す。
- `ikafssnhttpd` 自体は検索ロジックを持たない。
- 認証、アクセス制御、TLS 終端は前段のリバースプロキシ (Apache, nginx 等) に委ねる。`ikafssnhttpd` はこれらの機能を実装しない。
- `ikafssnhttpd` と `ikafssnserver` は別マシン上で動作する場合がある。この場合は TCP 接続を使用する。

**REST API エンドポイント (想定):**

| メソッド | パス | 説明 |
|---|---|---|
| POST | `/api/v1/search` | 検索リクエスト (クエリ配列を JSON ボディで送信) |
| GET | `/api/v1/info` | インデックス情報の取得 |
| GET | `/api/v1/health` | ヘルスチェック |

### 10.8 `ikafssnclient`

クライアントコマンド。`ikafssnserver` に直接ソケット接続するか、`ikafssnhttpd` に HTTP 接続して検索結果を取得する。出力形式は `ikafssnsearch` と完全に同一である。

```
ikafssnclient [options]

接続先 (いずれか):
  -socket <path>          ikafssnserver の UNIX ソケットパス
  -tcp <host>:<port>      ikafssnserver の TCP アドレス
  -http <url>             ikafssnhttpd の URL (例: http://example.com:8080)

必須:
  -query <path>           クエリ FASTA ファイル (- で標準入力)

オプション:
  -o <path>               出力ファイル (default: 標準出力)
  -k <int>                使用する k-mer サイズ (サーバが複数保持する場合)
  -min_score <int>        最小チェインスコア (default: サーバ側デフォルト)
  -max_gap <int>          チェイニング対角線ずれ許容幅 (default: サーバ側デフォルト)
  -max_freq <int>         高頻度 k-mer スキップ閾値 (default: サーバ側デフォルト)
  -min_diag_hits <int>    対角線フィルタ最小ヒット数 (default: サーバ側デフォルト)
  -stage1_topn <int>      Stage 1 候補数上限 (default: サーバ側デフォルト)
  -min_stage1_score <int> Stage 1 最小スコア閾値 (default: サーバ側デフォルト)
  -num_results <int>      最終出力件数 (default: 50)
  -seqidlist <path>       検索対象を指定アクセッションに限定
  -negative_seqidlist <path>  指定アクセッションを検索対象から除外
  -outfmt <tab|json>      出力形式 (default: tab)
  -v, --verbose           詳細ログ出力
```

`-seqidlist` / `-negative_seqidlist` 指定時、クライアントはファイルをローカルで読み込み、アクセッションリストをリクエストペイロードに含めてサーバに送信する。サーバ側で OID への解決とフィルタリングが行われる。

**使用例:**

```bash
# UNIX ソケット経由 (ローカル)
ikafssnclient -socket /var/run/ikafssn.sock -query query.fasta

# TCP 直接接続 (ローカルまたはリモート)
ikafssnclient -tcp 10.0.1.5:9100 -query query.fasta

# HTTP 経由 (リモート)
ikafssnclient -http http://search.example.com:8080 -query query.fasta

# seqidlist で検索対象を限定
ikafssnclient -socket /var/run/ikafssn.sock -query query.fasta -seqidlist targets.txt

# パイプラインで ikafssnretrieve に接続
ikafssnclient -socket /var/run/ikafssn.sock -query query.fasta | ikafssnretrieve -db nt > matches.fasta

# 特定の k-mer サイズを指定
ikafssnclient -socket /var/run/ikafssn.sock -query query.fasta -k 9
```

`ikafssnclient` の出力は `ikafssnsearch` と完全に同一であるため、後段の `ikafssnretrieve` は入力元がどちらであるかを意識する必要がない。

### 10.9 `ikafssninfo`

インデックス情報表示コマンド。インデックスファイルの統計情報を表示する。対応する BLAST DB が存在する場合は、DB の情報も併せて表示する。

```
ikafssninfo [options]

必須:
  -ix <dir>               インデックスディレクトリ

オプション:
  -db <path>              BLAST DB プレフィックス (指定時は DB 情報も表示)
  -v, --verbose           詳細ログ出力 (k-mer 頻度分布の詳細等)
```

**出力情報:**

- k-mer 長 (k) および k-mer 整数型 (uint16 / uint32)
- ボリューム数
- 各ボリュームの配列数、総 posting 数、ファイルサイズ (kix, kpx, ksx)
- 全体統計: 総配列数、総 posting 数、k-mer 出現頻度分布 (min, max, mean, median, percentiles)、圧縮率 (実サイズ / 無圧縮理論値)
- `-db` 指定時: BLAST DB のタイトル、配列数、総塩基数、ボリューム構成

---

## 11. サーバアーキテクチャ

### 11.1 全体構成

```
┌─────────────────────────────────────────────────────┐
│ マシン A (検索サーバ)                                 │
│                                                     │
│  ikafssnserver                                      │
│  ├── インデックス mmap (kix, kpx, ksx)              │
│  ├── UNIX ソケット ─────────────────┐               │
│  └── TCP ソケット (:9100) ──────────┤               │
│                                     │               │
│  ikafssnclient (ローカル) ──────────┘               │
│                                                     │
└─────────────────────────────┬───────────────────────┘
                              │ TCP
┌─────────────────────────────┴───────────────────────┐
│ マシン B (HTTP フロントエンド)                        │
│                                                     │
│  nginx (TLS 終端, 認証, レート制限)                  │
│      │                                              │
│      ▼                                              │
│  ikafssnhttpd (:8080)                               │
│      │                                              │
│      └── TCP で ikafssnserver に接続                 │
│                                                     │
│  ikafssnclient -http ─── ikafssnhttpd               │
│                                                     │
└─────────────────────────────────────────────────────┘
```

同一マシン構成の場合は UNIX ドメインソケットで接続し、TCP のオーバーヘッドを回避する。

### 11.2 ikafssnserver の内部構造

```
メインスレッド:
  ├── シグナルハンドラ (SIGTERM, SIGINT → グレースフルシャットダウン)
  ├── ソケットリスナー (accept ループ)
  └── 接続ディスパッチャ → ワーカースレッドプールへ

ワーカースレッドプール (TBB task_arena):
  ├── リクエスト受信・デシリアライズ
  ├── 検索実行 (Stage 1 → Stage 2)
  ├── レスポンスシリアライズ
  └── レスポンス送信

共有リソース (読み取り専用):
  ├── VolumeMap[] (kix, kpx, ksx の mmap)
  └── 設定パラメータ

スレッドローカルリソース:
  ├── score_per_seq[] (Stage 1 用)
  └── ヒットバッファ (Stage 2 用)
```

### 11.3 複数 k-mer サイズのサポート

`ikafssnserver` は起動時に `-ix` ディレクトリ内の全インデックスファイルを走査し、k-mer サイズごとにグループ化して読み込む。

```
nt_index/
  nt.00.09mer.kix, nt.00.09mer.kpx, nt.00.09mer.ksx   ← k=9 のインデックス
  nt.01.09mer.kix, nt.01.09mer.kpx, nt.01.09mer.ksx
  nt.00.11mer.kix, nt.00.11mer.kpx, nt.00.11mer.ksx   ← k=11 のインデックス
  nt.01.11mer.kix, nt.01.11mer.kpx, nt.01.11mer.ksx
```

クライアントはリクエスト内で k を指定する。指定がない場合はサーバのデフォルト k (起動時に設定、または利用可能な最大 k) を使用する。

### 11.4 グレースフルシャットダウン

```
SIGTERM / SIGINT 受信
  │
  ├── accepting_new = false  (新規接続の受付を停止)
  ├── キューイング済み未着手リクエストを破棄 (キューをフラッシュ)
  │
  ├── 実行中リクエスト完了を待機 (最大 shutdown_timeout 秒)
  │     └── 全完了 → 正常終了 (exit 0)
  │
  └── タイムアウト → 残存リクエストを強制終了
        └── exit 1
```

再起動時に旧プロセスのリクエスト状態を引き継ぐ必要はない。systemd の `Restart=on-failure` 等と組み合わせて運用する。

### 11.5 1プロセス1データベースの原則

`ikafssnserver` は1プロセスにつき1つの BLAST DB のインデックスのみを扱う。複数 DB を同時に提供する場合は、DB ごとに別プロセス (別ソケットパス / 別ポート) を起動し、`ikafssnhttpd` 側またはリバースプロキシ側でルーティングする。

```
# 例: nt と refseq を別プロセスで提供
ikafssnserver -ix ./nt_index    -socket /var/run/ikafssn_nt.sock
ikafssnserver -ix ./rs_index    -socket /var/run/ikafssn_rs.sock

# ikafssnhttpd でルーティング (パスベース)
ikafssnhttpd -server_socket /var/run/ikafssn_nt.sock -listen :8080 -path_prefix /nt
ikafssnhttpd -server_socket /var/run/ikafssn_rs.sock -listen :8081 -path_prefix /rs
# nginx 側で /nt → :8080, /rs → :8081 にルーティング
```

---

## 12. クライアント-サーバ間プロトコル

### 12.1 概要

`ikafssnserver` と `ikafssnclient` / `ikafssnhttpd` の間は自前バイナリプロトコルで通信する。このプロトコルは `ikafssnclient` がソケット経由で直接使用し、`ikafssnhttpd` が HTTP リクエストからの変換に使用する。

設計方針:
- シンプルな要求-応答型 (1リクエスト1レスポンス)
- フレーム長先行 (length-prefixed) で境界を明確化
- リトルエンディアン固定
- 拡張性のためにメッセージタイプとバージョンフィールドを含む

### 12.2 フレーム構造

```
┌──────────────────────────────────────────────────────┐
│ Frame Header (12 bytes)                              │
│ ├── magic:          uint32  "IKSV" (0x5653_4B49)    │
│ ├── payload_length: uint32  ペイロード長 (bytes)      │
│ ├── msg_type:       uint8   メッセージタイプ          │
│ ├── msg_version:    uint8   メッセージバージョン      │
│ └── reserved:       uint16  予約 (0)                 │
├──────────────────────────────────────────────────────┤
│ Payload (payload_length bytes)                       │
│ └── メッセージタイプ固有の内容                         │
└──────────────────────────────────────────────────────┘
```

### 12.3 メッセージタイプ

| msg_type | 方向 | 名前 | 説明 |
|---|---|---|---|
| 0x01 | C→S | SEARCH_REQUEST | 検索リクエスト |
| 0x81 | S→C | SEARCH_RESPONSE | 検索レスポンス |
| 0x02 | C→S | INFO_REQUEST | インデックス情報リクエスト |
| 0x82 | S→C | INFO_RESPONSE | インデックス情報レスポンス |
| 0x03 | C→S | HEALTH_REQUEST | ヘルスチェック |
| 0x83 | S→C | HEALTH_RESPONSE | ヘルスチェックレスポンス |
| 0xFF | S→C | ERROR_RESPONSE | エラーレスポンス |

クライアント→サーバは `0x01`〜`0x7F`、サーバ→クライアントは `0x80`〜`0xFF` の範囲を使用する。

### 12.4 SEARCH_REQUEST ペイロード

```
k:              uint8     使用する k-mer サイズ (0 = サーバデフォルト)
min_score:      uint16    最小チェインスコア (0 = サーバデフォルト)
max_gap:        uint16    チェイニング対角線ずれ許容幅 (0 = サーバデフォルト)
max_freq:       uint32    高頻度 k-mer スキップ閾値 (0 = サーバデフォルト)
min_diag_hits:  uint8     対角線フィルタ最小ヒット数 (0 = サーバデフォルト)
stage1_topn:    uint16    Stage 1 候補数上限 (0 = サーバデフォルト)
min_stage1_score: uint16  Stage 1 最小スコア閾値 (0 = サーバデフォルト)
num_results:    uint16    最終出力件数 (0 = サーバデフォルト)
seqidlist_mode: uint8     0 = フィルタなし, 1 = seqidlist (限定), 2 = negative_seqidlist (除外)
reserved:       uint8[2]  予約
num_seqids:     uint32    seqidlist 内のアクセッション数 (seqidlist_mode=0 の場合は 0)
seqids[]:       各アクセッション (num_seqids > 0 の場合):
  accession_len: uint16   アクセッション文字列長
  accession:     char[]   アクセッション文字列
num_queries:    uint16    クエリ配列数
queries[]:      各クエリ:
  query_id_len: uint16    クエリ ID の長さ (bytes)
  query_id:     char[]    クエリ ID 文字列
  seq_len:      uint32    塩基配列長
  sequence:     char[]    塩基配列 (ASCII, A/C/G/T/N)
```

seqidlist はリクエスト単位で指定する。クライアント側でファイルを読み込み、アクセッション文字列のリストとしてペイロードに含める。サーバ側で各ボリュームの `.ksx` を参照して OID に解決し、OID ビットセットを構築してフィルタリングに使用する。数万件程度の seqidlist であればリクエストサイズは数百 KB に収まり、実用上問題ない。

### 12.5 SEARCH_RESPONSE ペイロード

```
status:         uint8     0 = success, 非0 = error
k:              uint8     使用された k-mer サイズ
num_queries:    uint16    結果のクエリ数
results[]:      各クエリの結果:
  query_id_len: uint16    クエリ ID の長さ
  query_id:     char[]    クエリ ID 文字列
  num_hits:     uint16    ヒット数
  hits[]:       各ヒット:
    accession_len: uint16
    accession:     char[]
    strand:        uint8   0 = '+', 1 = '-'
    q_start:       uint32
    q_end:         uint32
    s_start:       uint32
    s_end:         uint32
    score:         uint16  チェインスコア
    volume:        uint16  ボリューム番号
```

### 12.6 ERROR_RESPONSE ペイロード

```
error_code:     uint32    エラーコード
message_len:    uint16    エラーメッセージ長
message:        char[]    エラーメッセージ (UTF-8)
```

---

## 13. 出力フォーマットとパイプライン

### 13.1 出力形式の統一

`ikafssnsearch` と `ikafssnclient` は完全に同一の出力フォーマットを生成する。これにより、後段のツール (`ikafssnretrieve` 等) は入力元を区別する必要がない。

### 13.2 タブ区切り形式 (デフォルト)

```
# query_id	accession	strand	q_start	q_end	s_start	s_end	score	volume
query_001	NM_001234	+	10	450	1020	1460	42	0
query_001	XM_005678	-	15	430	8050	8465	38	2
...
```

先頭に `#` で始まるヘッダ行を出力する。フィールド区切りはタブ文字 (`\t`) である。

### 13.3 JSON 形式

```json
{
  "results": [
    {
      "query_id": "query_001",
      "hits": [
        {
          "accession": "NM_001234",
          "strand": "+",
          "q_start": 10,
          "q_end": 450,
          "s_start": 1020,
          "s_end": 1460,
          "score": 42,
          "volume": 0,
          "volume": 0
        }
      ]
    }
  ]
}
```

### 13.4 パイプライン構成例

```bash
# 基本: ローカル検索 → ローカル BLAST DB から抽出
ikafssnsearch -ix ./index -query query.fasta | ikafssnretrieve -db nt > matches.fasta

# サーバ経由: クライアント → ローカル BLAST DB から抽出
ikafssnclient -socket /var/run/ikafssn.sock -query query.fasta | ikafssnretrieve -db nt > matches.fasta

# ファイル経由: 検索 → 保存 → 抽出
ikafssnsearch -ix ./index -query query.fasta -o results.tsv
ikafssnretrieve -db nt -results results.tsv -o matches.fasta

# HTTP 経由: リモート検索 → ローカル BLAST DB から抽出
ikafssnclient -http http://search.example.com:8080 -query query.fasta | ikafssnretrieve -db nt > matches.fasta

# NCBI efetch: ローカル BLAST DB なしでリモート取得
ikafssnclient -http http://search.example.com:8080 -query query.fasta | ikafssnretrieve -remote > matches.fasta

# NCBI efetch + API key: 高スループットリモート取得
ikafssnsearch -ix ./index -query query.fasta \
    | ikafssnretrieve -remote -api_key XXXXXXXX > matches.fasta
```

出力形式が統一されているため、検索の実行方法 (ローカル / ソケット / HTTP) と配列取得方法 (ローカル BLAST DB / NCBI efetch) は独立に選択でき、全ての組み合わせが動作する。

---

## 14. ソースツリー

```
ikafssn/
├── CMakeLists.txt
├── README.md
├── LICENSE
├── doc/
│   ├── ikafssn.en.md                      使用法ドキュメント (英語)
│   └── ikafssn.ja.md                      使用法ドキュメント (日本語)
├── src/
│   ├── ikafssnindex/
│   │   ├── CMakeLists.txt
│   │   └── main.cpp                        ikafssnindex エントリポイント
│   ├── ikafssnsearch/
│   │   ├── CMakeLists.txt
│   │   └── main.cpp                        ikafssnsearch エントリポイント
│   ├── ikafssnretrieve/
│   │   ├── CMakeLists.txt
│   │   ├── main.cpp                        ikafssnretrieve エントリポイント
│   │   ├── local_retriever.hpp             ローカル BLAST DB からの抽出
│   │   ├── local_retriever.cpp
│   │   ├── efetch_retriever.hpp            NCBI efetch リモート取得
│   │   └── efetch_retriever.cpp
│   ├── ikafssnserver/
│   │   ├── CMakeLists.txt
│   │   ├── main.cpp                        ikafssnserver エントリポイント
│   │   ├── server.hpp                      サーバメインループ
│   │   ├── server.cpp
│   │   ├── connection_handler.hpp          接続受付・ディスパッチ
│   │   ├── connection_handler.cpp
│   │   ├── request_processor.hpp           リクエスト処理 (検索実行)
│   │   └── request_processor.cpp
│   ├── ikafssnhttpd/
│   │   ├── CMakeLists.txt
│   │   ├── main.cpp                        ikafssnhttpd エントリポイント
│   │   ├── http_controller.hpp             Drogon コントローラ
│   │   ├── http_controller.cpp
│   │   ├── backend_client.hpp              ikafssnserver 接続クライアント
│   │   └── backend_client.cpp
│   ├── ikafssnclient/
│   │   ├── CMakeLists.txt
│   │   ├── main.cpp                        ikafssnclient エントリポイント
│   │   ├── socket_client.hpp               ソケット直接接続
│   │   ├── socket_client.cpp
│   │   ├── http_client.hpp                 HTTP 接続 (libcurl)
│   │   └── http_client.cpp
│   ├── ikafssninfo/
│   │   ├── CMakeLists.txt
│   │   └── main.cpp                        ikafssninfo エントリポイント
│   ├── core/
│   │   ├── types.hpp                       基本型定義 (Hit, ChainResult 等)
│   │   ├── config.hpp                      定数, デフォルト値
│   │   ├── kmer_encoding.hpp               2-bit encode/decode, revcomp (テンプレート)
│   │   └── varint.hpp                      LEB128 varint encode/decode
│   ├── index/
│   │   ├── kix_format.hpp                  .kix ヘッダ定義, 定数
│   │   ├── kix_writer.hpp                  .kix 書き出し
│   │   ├── kix_writer.cpp
│   │   ├── kix_reader.hpp                  .kix mmap 読み込み
│   │   ├── kix_reader.cpp
│   │   ├── kpx_format.hpp                  .kpx ヘッダ定義, 定数
│   │   ├── kpx_writer.hpp                  .kpx 書き出し
│   │   ├── kpx_writer.cpp
│   │   ├── kpx_reader.hpp                  .kpx mmap 読み込み
│   │   ├── kpx_reader.cpp
│   │   ├── ksx_format.hpp                  .ksx ヘッダ定義, 定数
│   │   ├── ksx_writer.hpp                  .ksx 書き出し
│   │   ├── ksx_writer.cpp
│   │   ├── ksx_reader.hpp                  .ksx 読み込み, accession→OID 逆引きマップ構築
│   │   ├── ksx_reader.cpp
│   │   ├── index_builder.hpp               構築ロジック統括 (テンプレート)
│   │   ├── index_builder.cpp
│   │   ├── partition_strategy.hpp          パーティション分割
│   │   └── buffer_manager.hpp              バッファ管理・フラッシュ
│   ├── search/
│   │   ├── seq_id_decoder.hpp              Stage 1 用 ID posting デコーダ
│   │   ├── posting_decoder.hpp             Stage 2 用 ID+pos 同時デコーダ
│   │   ├── oid_filter.hpp                  OID ビットセットフィルタ (seqidlist 用)
│   │   ├── oid_filter.cpp
│   │   ├── stage1_filter.hpp               Stage 1 候補絞り込み
│   │   ├── stage1_filter.cpp
│   │   ├── stage2_chaining.hpp             Stage 2 チェイニング DP
│   │   ├── stage2_chaining.cpp
│   │   ├── diagonal_filter.hpp             対角線フィルタ
│   │   ├── diagonal_filter.cpp
│   │   ├── volume_searcher.hpp             単一ボリューム検索
│   │   ├── volume_searcher.cpp
│   │   ├── multi_volume_search.hpp         マルチボリューム並列検索統括
│   │   └── multi_volume_search.cpp
│   ├── protocol/
│   │   ├── frame.hpp                       フレームヘッダ定義
│   │   ├── frame.cpp                       フレーム読み書き
│   │   ├── messages.hpp                    メッセージ型定義 (SearchRequest 等)
│   │   ├── serializer.hpp                  メッセージシリアライズ
│   │   └── serializer.cpp
│   ├── io/
│   │   ├── blastdb_reader.hpp              CSeqDB ラッパー
│   │   ├── blastdb_reader.cpp
│   │   ├── fasta_reader.hpp                クエリ FASTA 読み込み
│   │   ├── fasta_reader.cpp
│   │   ├── seqidlist_reader.hpp            seqidlist 読み込み (テキスト / バイナリ自動判別)
│   │   ├── seqidlist_reader.cpp
│   │   ├── mmap_file.hpp                   mmap RAII ラッパー
│   │   ├── mmap_file.cpp
│   │   ├── result_writer.hpp               出力フォーマッタ (tab/json)
│   │   ├── result_writer.cpp
│   │   ├── result_reader.hpp               検索結果パーサ (ikafssnretrieve 用)
│   │   └── result_reader.cpp
│   └── util/
│       ├── cli_parser.hpp                  コマンドライン解析
│       ├── cli_parser.cpp
│       ├── size_parser.hpp                 "8G", "512M" 等のサイズ文字列解析
│       ├── socket_utils.hpp                UNIX / TCP ソケットユーティリティ
│       ├── socket_utils.cpp
│       ├── progress.hpp                    進捗表示 (stderr)
│       └── logger.hpp                      ログ出力
└── test/
    ├── CMakeLists.txt
    ├── test_kmer_encoding.cpp              2-bit encode/decode, revcomp テスト
    ├── test_varint.cpp                     varint encode/decode テスト
    ├── test_kix_io.cpp                     .kix 読み書きラウンドトリップテスト
    ├── test_kpx_io.cpp                     .kpx 読み書きラウンドトリップテスト
    ├── test_ksx_io.cpp                     .ksx 読み書きラウンドトリップテスト
    ├── test_builder.cpp                    インデックス構築の結合テスト
    ├── test_stage1.cpp                     Stage 1 フィルタの単体テスト
    ├── test_diagonal_filter.cpp            対角線フィルタの単体テスト
    ├── test_chaining.cpp                   チェイニング DP の単体テスト
    ├── test_search_integration.cpp         検索の結合テスト
    ├── test_seqidlist_reader.cpp           seqidlist テキスト/バイナリ読み込みテスト
    ├── test_oid_filter.cpp                 OID フィルタ (限定/除外) のテスト
    ├── test_protocol.cpp                   プロトコルのシリアライズ/デシリアライズテスト
    ├── test_server_client.cpp              サーバ-クライアント結合テスト
    ├── test_result_reader.cpp              検索結果パーサのテスト
    ├── test_efetch_retriever.cpp           efetch URL 生成、レスポンスパース、バッチ戦略のテスト
    ├── testdata/
    │   ├── small.fasta                     テスト用小規模配列データ
    │   └── queries.fasta                   テスト用クエリ配列
    └── scripts/
        └── create_test_blastdb.sh          テスト用 BLAST DB 作成スクリプト
```

---

## 15. ビルドシステム

### 15.1 CMake 構成

```cmake
cmake_minimum_required(VERSION 3.16)
project(ikafssn VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
```

トップレベル CMakeLists.txt で共通ライブラリ (core, index, search, protocol, io, util) をビルドし、各コマンドの CMakeLists.txt で必要なライブラリのみをリンクする。

### 15.2 依存関係

| 依存 | 検出方法 | 使用コマンド | 備考 |
|---|---|---|---|
| NCBI C++ Toolkit | `find_package` またはパス指定 | index, search, retrieve, server, info | `libseqdb`, `libxobjutil` 等 |
| oneTBB | `find_package(TBB REQUIRED)` | index, search, server | `TBB::tbb` |
| Drogon | `find_package(Drogon REQUIRED)` | httpd | Drogon のビルドに jsoncpp, zlib 等が必要 |
| libcurl | `find_package(CURL REQUIRED)` | client (HTTP モード), retrieve (リモート取得) | |

NCBI C++ Toolkit の検出は環境によって煩雑になりうる。CMake 変数 `NCBI_TOOLKIT_DIR` でインストールパスを明示指定できるようにする。

### 15.3 コマンド別リンク依存

```cmake
# ikafssnindex: NCBI Toolkit + TBB
target_link_libraries(ikafssnindex PRIVATE core index io util
    ${NCBI_LIBS} TBB::tbb)

# ikafssnsearch: NCBI Toolkit + TBB
target_link_libraries(ikafssnsearch PRIVATE core index search io util
    ${NCBI_LIBS} TBB::tbb)

# ikafssnretrieve: NCBI Toolkit + libcurl (リモート取得時)
target_link_libraries(ikafssnretrieve PRIVATE core io util ${NCBI_LIBS})
if(ENABLE_REMOTE_RETRIEVE)
  target_link_libraries(ikafssnretrieve PRIVATE CURL::libcurl)
  target_compile_definitions(ikafssnretrieve PRIVATE IKAFSSN_ENABLE_REMOTE)
endif()

# ikafssnserver: NCBI Toolkit + TBB + protocol
target_link_libraries(ikafssnserver PRIVATE core index search protocol io util
    ${NCBI_LIBS} TBB::tbb)

# ikafssnhttpd: Drogon + protocol (NCBI Toolkit 不要)
target_link_libraries(ikafssnhttpd PRIVATE protocol util Drogon::Drogon)

# ikafssnclient: protocol (NCBI Toolkit, TBB 不要。HTTP モード時のみ libcurl)
target_link_libraries(ikafssnclient PRIVATE core protocol io util)
if(ENABLE_CLIENT_HTTP)
  target_link_libraries(ikafssnclient PRIVATE CURL::libcurl)
  target_compile_definitions(ikafssnclient PRIVATE IKAFSSN_ENABLE_HTTP)
endif()

# ikafssninfo: NCBI Toolkit (任意)
target_link_libraries(ikafssninfo PRIVATE core index io util ${NCBI_LIBS})
```

### 15.4 選択的ビルド

全コマンドの依存を揃えることが困難な環境 (例: Drogon が利用できないサーバ) のために、CMake オプションで個別コマンドのビルドを制御する。

```cmake
option(BUILD_HTTPD "Build ikafssnhttpd (requires Drogon)" ON)
option(BUILD_CLIENT "Build ikafssnclient" ON)
option(ENABLE_CLIENT_HTTP "Enable HTTP mode in ikafssnclient (requires libcurl)" ON)
option(ENABLE_REMOTE_RETRIEVE "Enable NCBI efetch in ikafssnretrieve (requires libcurl)" ON)
```

### 15.5 ビルド手順 (想定)

```bash
mkdir build && cd build
cmake .. -DNCBI_TOOLKIT_DIR=/path/to/ncbi-toolkit \
         -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
make test    # 単体テスト実行
```

---

## 16. サイズ見積もり

### 16.1 テーブル部分 (DB サイズに依存しない)

| k | 4^k | offsets (uint64) | counts (uint32) | 合計 |
|---|---|---|---|---|
| 8 | 65,536 | 512 KB | 256 KB | 768 KB |
| 9 | 262,144 | 2 MB | 1 MB | 3 MB |
| 10 | 1,048,576 | 8 MB | 4 MB | 12 MB |
| 11 | 4,194,304 | 32 MB | 16 MB | 48 MB |
| 12 | 16,777,216 | 128 MB | 64 MB | 192 MB |
| 13 | 67,108,864 | 512 MB | 256 MB | 768 MB |

kpx の `pos_offsets` テーブルも同サイズのため、kix + kpx のテーブル合計は上記の約 5/3 倍 (≈1.67 倍) となる (counts は kix のみ)。

### 16.2 posting 部分 (DB サイズに比例)

1つのボリューム内の有効塩基数を B とすると、k-mer 出現総数は約 B × 0.85 (N 除去後) である。

デルタ圧縮後の1エントリあたり平均バイト数の見積もり:
- ID posting (seq_id delta): 約 1.2 bytes
- Pos posting (pos delta/raw): 約 2.0 bytes
- 合計: 約 3.2 bytes/entry

| DB 有効塩基数 | k-mer 出現数 | 無圧縮 (8B/entry) | 圧縮後 (~3.2B/entry) |
|---|---|---|---|
| 1 GB | ~850M | ~6.8 GB | ~2.7 GB |
| 8 GB (nt 1vol) | ~6.8G | ~54 GB | ~22 GB |
| 400 GB (nt 全体) | ~340G | ~2.7 TB | ~1.1 TB |

### 16.3 .ksx

配列数 × (4 + 4 + 平均accession長) 程度。通常は全体の 0.1% 未満であり無視できる。

### 16.4 構築時のメモリ使用量

| コンポーネント | サイズ |
|---|---|
| counts[4^k] | 最大 256 MB (k=13) |
| バッファ | ユーザー指定 (default: 8 GB) |
| CSeqDB 内部バッファ | ~数百 MB |
| その他 | ~数十 MB |
| 合計 (default) | ~9 GB 程度 |

### 16.5 検索時のメモリ使用量

| コンポーネント | サイズ |
|---|---|
| kix mmap (per volume) | テーブル: ~768 MB + ID posting: 可変 (mmap なので常駐は一部) |
| kpx mmap (per volume) | テーブル: ~512 MB + pos posting: 可変 (遅延 mmap) |
| ksx (per volume) | ~数 MB |
| score_per_seq (per thread per volume) | num_sequences × 4 bytes |
| ヒットバッファ (per thread) | 動的 (通常 ~数 MB) |

mmap 領域は物理メモリに全てが常駐するわけではなく、OS のページキャッシュに管理される。物理メモリが十分あれば高速、不足時は SSD からのページインが発生する。

`ikafssnserver` での常駐運用時は、OS のページキャッシュにインデックスの大部分が載った状態が維持されるため、ローカルの `ikafssnsearch` よりも2回目以降のクエリレイテンシが改善されうる。

---

## 17. テスト計画

### 17.1 単体テスト

| テスト対象 | 主な検証項目 |
|---|---|
| `kmer_encoding` | 既知配列の 2-bit 変換が正しいこと。逆相補変換が `revcomp(revcomp(x)) == x` を満たすこと。N カウンタの動作。境界値 (k=5, k=8, k=9, k=13)。 |
| `varint` | varint encode → decode のラウンドトリップ。0, 1, 127, 128, 16383, 16384, UINT32_MAX の各境界値。 |
| `kix_io` | ヘッダの読み書きラウンドトリップ。テーブルの読み書き。デルタ圧縮 ID posting の encode → decode。空の posting list の扱い。 |
| `kpx_io` | pos posting のデルタ圧縮・復号。seq_id 境界でのリセット動作。 |
| `ksx_io` | 配列長テーブル、accession テーブルのラウンドトリップ。 |
| `stage1_filter` | 既知の小規模データで、正しい候補 seq_id が選出されること。max_freq によるスキップ動作。OID フィルタ適用時に対象外の seq_id がスコア集計から除外されること。 |
| `seqidlist_reader` | テキスト形式 (1行1アクセッション) の読み込み。先頭 `>` のトリム。空行・空白行のスキップ。バイナリ形式 (`blastdb_aliastool` 生成) の読み込み。テキスト/バイナリの自動判別。 |
| `oid_filter` | accession → OID 逆引きマップの構築が正しいこと。seqidlist モード (限定) で指定 OID のみ `pass()` が true を返すこと。negative_seqidlist モード (除外) で指定 OID のみ `pass()` が false を返すこと。フィルタ未指定時に全 OID が通過すること。存在しないアクセッションが警告されること。 |
| `diagonal_filter` | 対角線値の計算が正しいこと。min_diag_hits による除去が正しいこと。 |
| `chaining` | 完全一致ケース (チェインスコア = クエリ k-mer 数)。ギャップを含むケース。複数チェインの正しいスコア算出。max_gap による制限。 |
| `protocol` | フレームヘッダの読み書きラウンドトリップ。SearchRequest / SearchResponse のシリアライズ・デシリアライズ。不正フレーム (magic 不一致、長さ超過) の検出。 |
| `result_reader` | タブ区切り出力のパース。ヘッダ行のスキップ。不正行のハンドリング。 |
| `efetch_retriever` | バッチ accession リストの構築。efetch URL の生成 (全体取得 / 部分取得)。レスポンス FASTA のパースと accession マッチング。レート制限スリープ間隔の計算。HTTP エラーコードに応じたリトライ/スキップ判定。 |

### 17.2 結合テスト

| テスト | 内容 |
|---|---|
| Build → Search ラウンドトリップ | 小規模 FASTA から BLAST DB 作成 → インデックス構築 → 既知クエリで検索 → 期待される accession とスコアが得られることを確認 |
| 逆相補検索 | クエリの reverse complement が対象配列に含まれるケースで、正しく `-` strand として検出されること |
| マルチボリューム | 2 ボリュームに分割された BLAST DB に対して構築・検索が正しく動作すること |
| サーバ-クライアント | `ikafssnserver` を起動し、`ikafssnclient` からソケット経由で検索を実行して、`ikafssnsearch` と同一の結果が得られることを確認 |
| 出力一致性 | `ikafssnsearch` と `ikafssnclient` の出力が、同一クエリ・同一インデックスに対してバイト単位で一致すること |
| seqidlist フィルタ | `-seqidlist` 指定時に指定アクセッションのみが結果に含まれること。`-negative_seqidlist` 指定時に指定アクセッションが結果から除外されること。`ikafssnsearch` と `ikafssnclient` で同一の seqidlist を使用した場合に結果が一致すること。 |
| retrieve パイプライン | `ikafssnsearch` の出力を `ikafssnretrieve` にパイプし、正しい部分配列が FASTA として抽出されること |
| retrieve リモート取得 | (CI 外、手動) 既知の accession に対して `-remote` で NCBI efetch から配列を取得し、ローカル BLAST DB からの抽出結果と一致することを確認。バッチ取得と個別取得 (長大配列) の両方を検証。 |
| 大規模テスト | (CI 外、手動) 実際の BLAST DB (RefSeq select 等) に対して構築と検索を実行し、処理時間とメモリ使用量を計測 |

### 17.3 テストデータ

- `testdata/small.fasta`: 10〜20 本の短い配列 (100〜1000 bp)。手動で設計し、期待されるヒットパターンを制御する。
- `testdata/queries.fasta`: small.fasta の一部区間、逆相補、変異 (数塩基置換) を含むクエリ群。
- `testdata/seqidlist.txt`: small.fasta 内の一部アクセッション (5〜10 件) を含むテキスト形式の seqidlist。
- `scripts/create_test_blastdb.sh`: `makeblastdb` を呼び出してテスト用 BLAST DB を生成する。seqidlist のバイナリ版も `blastdb_aliastool` で作成する。

---

## 18. 開発フェーズ

### Phase 1: 基盤 (推定 2〜3 週間)

- コアライブラリ: `kmer_encoding`, `varint`, 基本型定義
- ファイル I/O: `.kix`, `.kpx`, `.ksx` の reader/writer
- mmap ラッパー
- 単体テスト: encoding, varint, ファイル I/O のラウンドトリップ
- CMake ビルドシステム (NCBI Toolkit, TBB の検出)

### Phase 2: インデックス構築 (推定 2〜3 週間)

- `blastdb_reader` (CSeqDB ラッパー)
- Phase 0〜4 の実装 (単一パーティション、バッファリングなし = 最単純ケース)
- 小規模テストデータでの構築テスト
- パーティション方式の追加
- バッファリング方式の追加
- `ikafssnindex` コマンドの実装

### Phase 3: ローカル検索 (推定 3〜4 週間)

- Stage 1 フィルタ
- Stage 2: posting デコーダ、対角線フィルタ、チェイニング DP
- seqidlist リーダー (テキスト / バイナリ自動判別)
- OID フィルタ (accession → OID 逆引き、ビットセット構築)
- `-seqidlist` / `-negative_seqidlist` オプション対応
- 単一ボリューム・単一クエリでの結合テスト
- FASTA reader
- 結果出力 (tab, json)
- `ikafssnsearch` コマンドの実装

### Phase 4: 並列化とマルチボリューム (推定 2〜3 週間)

- TBB によるマルチボリューム並列検索
- ボリューム自動検出
- 結果マージ
- マルチボリューム結合テスト
- 構築時の計数パス並列化

### Phase 5: サーバ・クライアント (推定 3〜4 週間)

- バイナリプロトコルのフレーム/メッセージ定義と実装 (seqidlist フィールド含む)
- プロトコルのシリアライズ・デシリアライズテスト
- `ikafssnserver` の実装 (ソケットリスナー、接続ハンドラ、ワーカープール)
- サーバ側 seqidlist → OID 解決の実装
- グレースフルシャットダウンの実装
- 複数 k-mer サイズのインデックス同時ロード
- `ikafssnclient` のソケット直接接続モード (seqidlist ファイル読み込み・送信対応)
- サーバ-クライアント結合テスト (出力一致性の検証、seqidlist フィルタの検証)

### Phase 6: HTTP、retrieve、リモート取得 (推定 3〜4 週間)

- `ikafssnhttpd` の実装 (Drogon コントローラ、バックエンド接続)
- REST API エンドポイントの実装
- `ikafssnclient` の HTTP 接続モード (libcurl)
- `ikafssnretrieve` の実装 (結果パーサ、ローカル BLAST DB 部分配列抽出)
- `ikafssnretrieve` の NCBI efetch リモート取得機能:
  - efetch URL 構築 (バッチ / 個別部分取得)
  - レート制限制御 (API key 有無に応じたスリープ)
  - レスポンス FASTA パースと accession マッチング
  - 指数バックオフ付きリトライ
  - 配列長に基づくバッチ取得 / 個別取得の自動切替
- パイプライン結合テスト (ローカル / リモート)

### Phase 7: 補助ツールと最適化 (推定 2〜3 週間)

- `ikafssninfo` コマンドの実装
- プロファイリングとボトルネック解消
- madvise チューニング (MADV_SEQUENTIAL, MADV_RANDOM の適切な使用)
- max_freq 自動調整のキャリブレーション
- 使用法ドキュメントの整備 (`doc/ikafssn.en.md`, `doc/ikafssn.ja.md`)
- systemd ユニットファイルのサンプル

合計推定: 17〜24 週間

---

## 19. 将来の拡張

以下は初期リリースには含めないが、将来的に追加を検討する項目である。

### 19.1 minimizer 対応

全 k-mer の代わりに minimizer (ウィンドウ内の最小 k-mer) のみを索引化することで、インデックスサイズを大幅に削減する。検索感度とのトレードオフがある。

### 19.2 構築の外部ソートモード

バッファサイズとパーティション数の組み合わせでもメモリに収まらない極端なケースで、一時ファイルへの中間フラッシュとマルチウェイマージを行う構築モード。

### 19.3 Spaced seed

k-mer の代わりに spaced seed (特定位置のみマッチさせるパターン) を用いることで、置換に対する頑健性を向上させる。

### 19.4 タンパク質配列対応

アルファベットサイズが 20 になり、エンコーディングとテーブルサイズの設計が大幅に変わる。5-bit 表現で k=6 なら 30 bits (uint32) に収まる。

### 19.5 seq_id の uint64 拡張

ヘッダの `SEQ_ID_WIDTH` フラグにより uint64 の seq_id をサポートする。posting entry のサイズが変わるため、デコーダの型パラメータ化が必要。

### 19.6 分散検索

複数マシンにボリュームを分散配置し、`ikafssnserver` を各マシンで起動する。`ikafssnhttpd` またはその上位のディスパッチャが複数サーバにリクエストを分配し、結果をマージする。現行のマルチボリューム結果マージロジックの延長で実現可能。

### 19.7 インクリメンタル更新

既存インデックスに新規配列を追加する機能。全再構築を回避するため、追加分の posting を別セグメントに格納し、検索時にマージする。

---

*以上*
