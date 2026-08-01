# ikafssn ユーザガイド

**ikafssn** (Independent programs of K-mer-based Alignment-Free Similarity Search for Nucleotide sequences) は、NCBI BLAST DB の塩基配列に対して転置インデックスを構築し、k-mer マッチングとコリニアチェイニングによるアライメントフリー類似配列検索を行うツール群です。

- **プライマリリポジトリ**: <https://github.com/astanabe/ikafssn>
- **セカンダリリポジトリ**: <https://gitlab.com/astanabe/ikafssn>

## 概要

ikafssn は 7 つの独立したコマンドラインプログラムで構成されます。

| コマンド | 概要 |
|---|---|
| `ikafssnindex` | BLAST DB から k-mer 転置インデックスを構築 |
| `ikafssnsearch` | ローカル直接検索 (mmap インデックス) |
| `ikafssnretrieve` | マッチした部分配列の抽出 |
| `ikafssnserver` | 検索デーモン (UNIX/TCP ソケット) |
| `ikafssnhttpd` | ikafssnserver の HTTP REST プロキシ |
| `ikafssnclient` | クライアント (ソケットまたは HTTP) |
| `ikafssninfo` | インデックス / DB / サーバ情報表示 |

各コマンドは必要な依存のみをリンクする独立した実行ファイルです。

## クイックスタート

```bash
# 1. BLAST DB からインデックスを構築
ikafssnindex -db mydb -k 11 -o ./index

# 2. クエリ FASTA で検索
ikafssnsearch -ix ./index/mydb -query query.fasta

# 3. マッチした部分配列を抽出
ikafssnsearch -ix ./index/mydb -query query.fasta | ikafssnretrieve -db mydb > matches.fasta
```

## コマンド

### 共通オプション

すべてのコマンドで以下のオプションが利用可能です。

```
  -h, --help              使用法の表示 (バージョンヘッダ付き)
  -version, --version     バージョン・ビルド情報の表示
  -v, --verbose           詳細出力
```

`-version` の出力形式:

```
<コマンド名>: <バージョン>
 Package: ikafssn <バージョン>, build <ISO 8601 タイムスタンプ>
```

### 圧縮入出力

`ikafssnsearch`、`ikafssnclient`、`ikafssnretrieve` は FASTA クエリ入力および
TSV / JSON / FASTA 出力に対し、透過的なコーデックラッパを通します。

- **入力**は先頭バイトのマジック判定で自動検出されるため、`.gz` / `.bz2`
  / `.xz` / `.zst` の FASTA ファイル（およびそれらをパイプした標準入力）
  はユーザー側でフラグを付けることなく展開されます。連結ストリーム
  (`cat a.fa.gz b.fa.gz > c.fa.gz`) も 4 種すべてのコーデックで正しく展開されます。
- **出力**はパス末尾の拡張子でコーデックを選択します
  (`out.tsv.gz`、`out.json.zst`、`out.fa.bz2` など)。`-o` 未指定または
  `-o -` の場合は他の設定によらず標準出力には無圧縮のまま書き込まれます。
  標準出力を圧縮したい場合は外部エンコーダにパイプしてください。
- `-compression_level <int>` で選択したコーデックの圧縮レベルを指定します。
  各コーデックの受理範囲と既定値: gzip 0..9 (既定 6)、bzip2 1..9 (既定 9)、
  xz 0..9 (既定 6)、zstd は `ZSTD_minCLevel()..ZSTD_maxCLevel()` (既定 3)。
  範囲は I/O オープン前にパース時点で検証されます。
- SAM / BAM 出力は 4 種の圧縮拡張子を **拒否** します — BAM 自体が
  既に BGZF であり、二重ラップを意図的に避けるためです。圧縮された
  結果ファイルが必要な場合は `-output_format tsv` または `-output_format json` を
  指定してください。

### ikafssnindex

BLAST DB から k-mer 転置インデックスを構築します。各ボリュームに対して `.kix` (ID ポスティング)、`.kpx` (位置ポスティング、`-mode 1` の場合は省略)、`.ksx` (配列メタデータ) のファイルを生成します。`-max_freq_build` 使用時は共有 `.khx` (構築時除外ビットセット) も生成されます。`.khx` ファイルは全ボリューム共通 (k 値ごとに 1 つ) です。

**インデックス単位:** 各親 BLAST OID は構築時に 1 つ以上の **fragment**（フラグメント）に分割されます。フラグメントは `.kix` / `.kpx` / `.ksx` で 1 つの内部 SeqId として登録される単位です。同じ親由来の隣接フラグメントは `-overlap_length` 塩基を共有するため、境界をまたぐチェインヒットも、必ずどちらか一方のフラグメント内に完全に収まります。`-min_length_split 0` を指定した場合は、各親に対して親全体をカバーする 1 つのフラグメントだけが登録される退化レイアウトになります。検索側の dedup によりフラグメント相対座標は **親相対座標** に畳み戻されるため、下流ツールには親 OID あたり 1 行の正準ヒットだけが届きます。

```
ikafssnindex [options]

必須:
  -db <path>              BLAST DB プレフィックス
  -k <int>                k-mer 長 (5〜15)
  -o <dir>                出力ディレクトリ

オプション:
  -mode <1|2|3>           インデックスがサポートする検索モード (デフォルト: 1)
                          1 = Stage 1 のみ (.kpx 生成スキップ、ディスク・時間節約。デフォルト)
                          2 = Stage 1+2
                          3 = Stage 1+2+3 (インデックス構築では 2 と同一動作)
  -min_seq_length <int>   インデックスに登録する親 OID の最短長 (デフォルト: 64)。
                          これより短い親はインデックス時に弾かれ、フラグメント
                          としても登録されません。値は .kix / .kpx / .ksx
                          ヘッダに記録され、ikafssnsearch / ikafssnclient は
                          自身の -min_query_length とこの値の整合性を起動時
                          にチェックします（短すぎるクエリが silently に 0 ヒット
                          になるのを防ぐため）。
  -min_length_split <int> フラグメント分割の閾値 (デフォルト: -mode 1 のとき
                          50000、-mode 2/3 のとき 0)。これを超える親は
                          kafssstore の calcsegment2 公式に従って複数の
                          重なり付きフラグメントに分割されます (DNA2 モード:
                          ncbi4na==0xF / N ラン位置で親を有効セグメントに
                          分け、各セグメントを更にこの長さに近づけて分割)。
                          0 を指定すると分割を無効化し、各親はそのまま 1
                          フラグメントとして登録されます。
  -overlap_length <int>   隣接フラグメントの末尾オーバーラップ塩基数
                          (デフォルト: -mode 1 のとき 500、それ以外で 0)。
                          有効セグメント内の隣接フラグメントはこの塩基数を
                          共有します。分割を有効化する場合は
                          0 < overlap_length < min_length_split / 2 を満たす
                          必要があります。検索時にはこの値が許容クエリ長の
                          上限となり (これより長いクエリは
                          kSkipQueryTooLong でスキップされる)、親相対の
                          dedup キーが必ず隣接 2 フラグメント以内に収まる
                          ことを保証します。
  -memory_limit <size>    ボリュームごとのソートバッファ予算 (デフォルト: 物理メモリの半分)
                          接尾辞 K, M, G を認識。ポスティングパスで k-mer
                          エントリをこの予算に収まるようパーティション分割
                          することで、ボリュームあたりのピーク RAM を制限する。
  -max_freq_build <num>   構築時高頻度 k-mer 除外閾値 (ボリューム横断集計)
                          1 または 1.0: 無効 (除外なし、デフォルト)
                          0〜1 未満: 全ボリューム合計 NSEQ に対する割合
                          1 超: 絶対カウント閾値
                          0: エラー (使用不可)
                          カウントは全ボリュームで合算後にフィルタリング
  -freq_threshold_part <int>
                          .kpx の (k-mer, seq_id) 単位パーティション閾値 (デフォルト: 8、上限: 255)
                          1 つの (k-mer, seq_id) クラスタの出現回数がこの値を超えると
                          独立パーティションとして切り出し、それ以下のクラスタは共有
                          short bucket にまとめる (染色体級配列のブロック幅が絶対位置の
                          大きさに引きずられないようにし、`.kpx` の圧縮率を改善)
  -nthread_highfreq_filter <int>
                          ボリューム横断高頻度フィルタリングのスレッド数
                          (デフォルト: min(8, nthread))
  -max_degen_expand <int> 縮重塩基展開の最大数/k-mer (デフォルト: 4、最大: 16、0/1: 無効)
                          IUPAC 縮重塩基を含む k-mer から生成する非縮重 k-mer の最大数を制御。
                          各位置の変異数の積がこの上限以下の場合に展開を実行。
  -t <int>                スペースドシード用テンプレート長 (デフォルト: 0)
                          0: 連続 k-mer (従来方式)
                          13, 15, 18: スペースドシードテンプレート長 (-k 8 または 9 が必要)
                          16, 18, 21: スペースドシードテンプレート長 (-k 11 または 12 が必要)
  -template_type <str>    スペースドシードのテンプレート種別 (-t 指定時は必須)
                          coding: コーディングテンプレートのみ
                          optimal: オプティマルテンプレートのみ
                          both: coding と optimal のインデックスを順次構築
  -nthread <int>          使用スレッド数 (デフォルト: 利用可能な全コア)
                          ボリューム内の計数・パーティションスキャン・
                          ソート、メタデータパスの OID 走査を並列化。
                          ボリュームは 1 つずつ順次処理され、各ボリュームに
                          -nthread 全数が割り当てられる。
  -force_rebuild <0|1>    完全再構築を強制する (デフォルト: 0)
                          0 = ボリューム単位の resume: 既存の .kix /
                              .kpx / .ksx を厳密検査 (ファイルサイズ、
                              magic / format_version、全 posting list
                              ウォーク) し、合格したボリュームは再利用、
                              失敗または欠落したものは再構築する。
                              前回中断時の .ksx.tmp / .kix.tmp /
                              .kpx.tmp も再利用対象。ビルド完了後の
                              バリデーション失敗時は該当パスを表示し
                              て非ゼロ終了する (自動削除はしない)。
                          1 = 該当ビルドパラメータの tmp / final
                              ファイルを全削除し、全ボリュームを再構築。
  -v, --verbose           詳細ログ出力
```

**使用例:**

```bash
# 小規模 DB、メモリ豊富
ikafssnindex -db mydb -k 11 -o ./index -memory_limit 128G

# 大規模 DB、メモリ制限、マルチスレッド
ikafssnindex -db nt -k 11 -o ./nt_index -memory_limit 32G -nthread 32

# 高頻度 k-mer を除外して構築 (絶対値指定)
ikafssnindex -db nt -k 11 -o ./nt_index -max_freq_build 50000

# 全ボリューム合計配列数の 1% を超える k-mer を除外
ikafssnindex -db nt -k 11 -o ./nt_index -max_freq_build 0.01

# mode 1 インデックスを構築 (Stage 1 のみ、.kpx なし)
ikafssnindex -db mydb -k 11 -o ./index -mode 1

# スペースドシードインデックスを構築 (k=11, t=16, coding テンプレート)
ikafssnindex -db mydb -k 11 -t 16 -template_type coding -o ./index

# コーディングテンプレートのみでスペースドシードインデックスを構築
ikafssnindex -db mydb -k 11 -t 18 -template_type coding -o ./index

# coding と optimal の両方のインデックスを順次構築
ikafssnindex -db mydb -k 11 -t 18 -template_type both -o ./index

# k=12, t=21 でスペースドシードインデックスを構築
ikafssnindex -db mydb -k 12 -t 21 -o ./index

# PCR 最適化スペースドシードインデックスを構築 (k=8, t=13)
ikafssnindex -db mydb -k 8 -t 13 -o ./index
```

### ikafssnsearch

ローカル直接検索コマンドです。インデックスファイルを直接 mmap して検索します。サーバプロセスを必要としません。

```
ikafssnsearch [options]

必須:
  -ix <prefix>            インデックスプレフィックス (blastn -db と同様)
  -query <path>           クエリ FASTA ファイル (- で標準入力)

オプション:
  -k <int>                使用する k-mer サイズ (複数の k 値が存在する場合は必須)
  -o <path>               出力ファイル (デフォルト: 標準出力)
  -nthread <int>          検索と出力整形のスレッド数
                          (デフォルト: 利用可能な全コア)
  -memory_limit <size>    検索メモリ予算 (デフォルト: 物理メモリの半分)
                          .khx/.ksx メタデータを常駐させ、残余で Stage 3 の
                          posting heap を制限。接尾辞 K, M, G を認識
  -mode <1|2|3>           検索モード (デフォルト: 1)
                          1=Stage 1 のみ、2=Stage 1+2、3=Stage 1+2+3
  -db <path>              モード 3 用 BLAST DB パス (デフォルト: -ix と同じ)
  -min_query_length <int> クエリ配列の最短長 (デフォルト: 64)。これより短い
                          クエリは kSkipQueryTooShort でスキップされます。
                          値はインデックスの min_seq_length (.kix / .ksx
                          ヘッダに記録) 以上である必要があり、それより小さい
                          値を指定するとインデックスが登録していない短い
                          配列までクエリ対象に含めてしまうため、整合性
                          チェックで起動が拒否されます。
  -stage1_max_nhit_per_subject <int>  parent ごとの Stage 1 候補数上限
                          (デフォルト: 1、0=無制限)。
                          (query, parent, volume, strand) 単位で適用
  -stage1_max_nhit_per_subject_mode <1|2|3|4>  per-subject 選別モード
                          (デフォルト: 3); 1/2=上位N (タイ非包含)、3/4=上位N + タイ包含
  -stage1_max_nhit_per_volume <int>  (query, volume, strand) ごとの Stage 1
                          候補数上限 (デフォルト: 0=無制限)
  -stage1_max_nhit_in_total <int>  全 volume・strand を通したクエリごとの
                          Stage 1 候補数上限 (デフォルト: 0=無制限)。これを
                          設定すると、-stage1_max_nhit_per_volume を明示しない
                          限り同じ値が設定されます (それ以上である必要があります)。
  -stage1_min_score <num> Stage 1 最小スコア閾値 (デフォルト: 0.5)
                          整数 (>= 1): 絶対閾値
                          小数 (0 < P < 1): クエリ k-mer に対する割合
                            クエリごとに ceil(Nqkmer * P) - Nhighfreq に解決
  -stage2_min_score <int> 最小チェインスコア (デフォルト: 0 = 適応的)
                          0 = Stage 1 の解決済み閾値を最小値として使用
                          1 以上: 絶対的な最小チェインスコア
  -stage2_max_gap <int>   チェイニング対角線ずれ許容幅 (デフォルト: 100)
  -stage2_max_lookback <int>  チェイニング DP 探索窓サイズ (デフォルト: 64、0=無制限)
  -stage2_max_nhit_per_subject <int>  サブジェクトあたりの最大チェイン数 (デフォルト: 1、0=無制限)
  -stage2_max_nhit_per_subject_mode <1|2|3|4>  サブジェクト単位の選択モード (デフォルト: 3)
                           1/2=上位N位 (タイ非包含)、3/4=上位N位+スコアタイ包含;
                           1/3=parent ごとに strand 統合、2/4=strand 分離
  -stage2_max_nhit_in_total <int>  全ボリュームを通したクエリあたりの最大チェイン数、
                          chainscore 基準 (デフォルト: 0=無制限)
  -stage2_min_nhit_diag <int>  対角線フィルタ最小ヒット数 (デフォルト: 1)
  -context_extend <value>        モード 3 のコンテクスト拡張 (デフォルト: 2.0)
                          整数: 拡張する塩基数; 小数: クエリ長に対する倍率
  -stage3_traceback <0|1> モード 3 でトレースバックを有効化 (デフォルト: 0)
  -stage3_gapopen <int>   モード 3 のギャップオープンペナルティ (デフォルト: 10)
  -stage3_gapext <int>    モード 3 のギャップ伸長ペナルティ (デフォルト: 1)
  -stage3_min_ppositive <num>  モード 3 の最小正スコア率フィルタ (デフォルト: 0)
  -stage3_min_npositive <int>  モード 3 の最小正スコア塩基数フィルタ (デフォルト: 0)
  -stage3_max_nhit_per_subject <int>  サブジェクトあたりの最大ヒット数、alnscore 基準 (デフォルト: 1、0=無制限)
  -stage3_max_nhit_per_subject_mode <1|2|3|4>  サブジェクト単位の選択モード (デフォルト: 3)
                           1/2=上位N位 (タイ非包含)、3/4=上位N位+スコアタイ包含;
                           1/3=parent ごとに strand 統合、2/4=strand 分離
  -stage3_max_nhit_in_total <int>  全ボリュームを通したクエリあたりの最大ヒット数、
                          alnscore 基準 (デフォルト: 0=無制限)
  -stage3_score_matrix <str>  モード 3 のスコア行列 (デフォルト: degmatch)
                          利用可能: degmatch、dnafull、nuc44
  -seqidlist <path>       検索対象を指定アクセッションに限定
  -negative_seqidlist <path>  指定アクセッションを検索対象から除外
  -strand <-1|1|2>       検索する鎖 (デフォルト: 2)
                          1=プラス鎖のみ、-1=マイナス鎖のみ、2=両鎖
  -accept_qdegen <0|1>    縮重塩基を含むクエリを許可 (デフォルト: 1)
  -max_degen_expand <int> 縮重塩基展開の最大数/k-mer (デフォルト: 16、最大: 256、0/1: 無効)。
                          これはフィルタではなくクエリ用パラメータであり、構築時の
                          max_degen_expand だけが異なる索引変種が複数残った場合の
                          タイブレークにも使われる (構築値がこのクエリ値と一致する変種を
                          優先し、なければ最大の構築値を持つ変種を選択)。
  -t <int>                スペースドシード用テンプレート長 (デフォルト: 0 = 連続)
                          0: 連続 k-mer (従来方式)。これは既定値でありワイルドカードではない
                          13, 15, 18: スペースドシードテンプレート長 (-k 8 または 9 が必要)
                          16, 18, 21: スペースドシードテンプレート長 (-k 11 または 12 が必要)
  -template_type <str>    スペースドシードのテンプレート種別 (デフォルト: both; -t 0 とは併用不可)
                          coding: coding インデックスのみ使用
                          optimal: optimal インデックスのみ使用
                          both: coding と optimal のインデックスを検索時にマージ
  -min_seq_length <int>   この min_seq_length で構築された索引変種を選択 (デフォルト: 任意)
  -min_length_split <int> この min_length_split で構築された索引変種を選択 (デフォルト: 任意)
  -overlap_length <int>   この overlap_length で構築された索引変種を選択 (デフォルト: 任意)
  -max_freq_build <int>   この max_freq_build で構築された索引変種を選択 (デフォルト: 任意)
  -output_format <tsv|json|sam|bam>  出力形式 (デフォルト: tsv)
  -compression_level <int> 出力圧縮レベル (gzip/bzip2/xz/zstd; 既定値: gzip=6, bzip2=9, xz=6, zstd=3)
                          コーデックは -o の拡張子 (.gz/.bz2/.xz/.zst) で選択。SAM/BAM では拒否
  -v, --verbose           詳細ログ出力
```

(注: -query および -primer は gzip/bzip2/xz/zstd 圧縮された FASTA を
 先頭マジックバイトで自動判定するため、フラグ指定は不要です。)

`-ix` オプションには `blastn -db` と同様にインデックスのプレフィックスパス (拡張子なし) を指定します。例えば、マルチボリューム BLAST DB (`nt`、ボリューム `nt.00`、`nt.01`) に対して `ikafssnindex -db nt -k 11 -o /data/index` が以下のファイルを生成した場合:

```
/data/index/nt.00.k11.minlen64.minsplit50000.ovllen500.kix
/data/index/nt.00.k11.minlen64.minsplit50000.ovllen500.ksx
/data/index/nt.01.k11.minlen64.minsplit50000.ovllen500.kix
/data/index/nt.01.k11.minlen64.minsplit50000.ovllen500.ksx
/data/index/nt.k11.minlen64.minsplit50000.ovllen500.kvx
```

`-ix /data/index/nt` と指定します。プレフィックス `/data/index/nt` はディレクトリ `/data/index/` とベース名 `nt` に分割され、`.kvx` マニフェストファイルからボリュームベースネームの一覧を読み取ってボリュームを検出します。集約データベース (例: `combined` が `foo` と `bar` を束ねたもの) の場合、インデックスファイルは `foo.k11.…kix`、`bar.k11.…kix`、マニフェストは `combined.k11.…kvx` となります。

**索引変種の選択。** 索引変種は、ファイル名に符号化された 8 つのパラメータ `k`、`t`、`template_type`、`min_seq_length`、`min_length_split`、`overlap_length`、`max_freq_build`、`max_degen_expand` で識別されます。`ikafssnsearch` はプレフィックス配下の全変種を列挙し、上記オプションで絞り込みます。各オプションは無指定ならワイルドカード (任意の値に一致)、指定すればその値での厳密フィルタとして働きます。フィルタ後、残った変種は最初の 7 パラメータ (`max_degen_expand` を除く) で一意に定まらなければならず、そうでなければ「index filter is ambiguous」エラー (絞り込みに使えるオプションを列挙) で、一件も一致しなければ「no index variant matches」エラーで検索を中止します。したがって変種が 1 つしかなければ全オプションを省略でき、複数あれば一意に定まるだけのオプションを与えれば十分です。

`max_degen_expand` だけは一意でなくても許容されます。`max_degen_expand` のみが異なる変種が複数残った場合は、上記のとおり `-max_degen_expand` でタイブレークされます。構築時の `-max_degen_expand` (`ikafssnindex`) とクエリ時の `-max_degen_expand` は一致しなくても構いません。

template_type の解決は次のとおりです:

| オプション | 結果 |
|---|---|
| `-t 0`(または `-t` 省略) | 連続変種(単一) |
| `-t>0 -template_type coding` | coding 変種(単一) |
| `-t>0 -template_type optimal` | optimal 変種(単一) |
| `-t>0 -template_type both` | coding+optimal ペア(both 必須。片方欠落でエラー) |
| `-t>0`(`-template_type` 無指定) | 両方存在すれば coding+optimal ペア、片方のみ存在すればその単独(cod のみ→coding、opt のみ→optimal) |

both ペアの場合、coding と optimal は `template_type` 以外の全識別パラメータ(**`max_degen_expand` を含む**)で一致しなければなりません(ちゃんぽん不可)。共有する `max_degen_expand` は**両側に存在する値**の中でタイブレーク(query 値優先・なければ共通の最大値)し、両側で共通する `max_degen_expand` が無ければ混在させずにエラーで中止します。

注意すべき非対称が 2 点あります。第一に、`-t` の既定値は `0` (連続インデックス) であり、これはワイルドカードではなく具体値です。`-t` を省略すると連続変種が選択され、スペースドシード変種には一致しません (一方 `ikafssninfo` と `ikafssnserver` では無指定の `-t` をワイルドカードとして扱います)。第二に、`-t 0` と `-template_type` の併用は、連続インデックスにテンプレート種別が存在しないためエラーになります。

`-accept_qdegen` が 0 の場合、IUPAC 縮重塩基 (R, Y, S, W, K, M, B, D, H, V, N) を含むクエリは `degen_rejected` のスキップマーカー (TSV: `*SKIPPED:degen_rejected`、JSON: `"status": "skipped"`、SAM: unmapped レコード + `XR:Z:degen_rejected`) と stderr 警告付きでスキップされ、終了コードは 2 になります。詳しい理由一覧と各フォーマットでの表現は下記「スキップ理由」を参照してください。`-accept_qdegen 1` を指定すると縮重塩基を含むクエリも受け付けます。1 文字の縮重塩基を含む k-mer は全可能バリアントに展開して検索に使用されます (例: R→A,G で 2 k-mer、N→A,C,G,T で 4 k-mer)。位置ごとの展開積が `-max_degen_expand` を超える窓は emit されません。この場合、クエリごとに 1 回、当該クエリ名と当該 k-mer がスキップされた旨の警告が標準エラーに出力されます (注: これは窓単位の unemit であり、クエリ全体のスキップではありません — 残りの窓は検索に使われ、emit されなかった位置はフラクショナル閾値の `Nhighfreq` に反映されます)。サーバモード (`ikafssnserver`) ではこの警告がプロトコル経由で `ikafssnclient` に伝播され、クライアント側でも同じメッセージが表示されます。この処理はインデックス構築時のサブジェクト配列に対する処理と同等です。

`-seqidlist` と `-negative_seqidlist` は排他的 (同時指定不可) です。ファイル形式はテキスト (1 行 1 アクセッション) と `blastdb_aliastool -seqid_file_in` で生成されるバイナリ形式の両方を受け付け、先頭のマジックバイトで自動判別します。

後段ステージでのみ有効なオプションを、それより前の `-mode` とともに指定した場合は、黙って無視せずエラー (非ゼロ終了) になります。Stage 2 オプション (`-stage2_min_score`、`-stage2_max_gap`、`-stage2_max_lookback`、`-stage2_min_nhit_diag`、`-stage2_max_nhit_per_subject`、`-stage2_max_nhit_per_subject_mode`、`-stage2_max_nhit_in_total`) は `-mode 2` 以上を、Stage 3 オプション (`-stage3_traceback`、`-stage3_gapopen`、`-stage3_gapext`、`-stage3_min_ppositive`、`-stage3_min_npositive`、`-stage3_score_matrix`、`-stage3_max_nhit_per_subject`、`-stage3_max_nhit_per_subject_mode`、`-stage3_max_nhit_in_total`、`-context_extend`、`-db`) は `-mode 3` を必要とします。同じ規則は `ikafssnclient` (CLI) と `ikafssnhttpd` (JSON リクエストボディ) でも適用されます。

`ikafssnindex` も矛盾するオプションの組み合わせをエラーにします: `-template_type` は `-t > 0` が必要、`-freq_threshold_part` は `-mode 2` または `3` が必要、`-nthread_highfreq_filter` は `-max_freq_build` の有効化が必要です。`ikafssnretrieve` はローカル `-db` 実行時に efetch 専用オプション (`-api_key`、`-batch_size`、`-max_nretry`、`-timeout`、`-range_threshold`) を拒否します。`ikafssnclient` と `ikafssninfo` では、HTTP 認証オプション (`-user`、`-http_user`、`-http_password`、`-netrc_file`) は `-http` が必要で、3 つの方式は相互排他であり、`-http_password` は `-http_user` を必要とします。

**所要時間の出力.** ステージごとの件数に加えて、実行のたびに 2 行の所要時間が標準エラー出力に書かれます。いずれも `info` レベルなので `-v` なしで出力されます。

```
[INFO] Timing run_search (s): s1_open=0.512 s1_compute=12.345 s1_fold=0.031 s1_intotal=0.002 s2_open=0.220 s2a=45.678 s2b=3.210 s2_free=0.015 dedup=0.012 parent_topn=0.034 s2_intotal=0.005 total=62.064
[INFO] Timing overall (s): open_index=1.234 preprocess=0.456 run_search=62.064 convert=2.345 stage3=0.000 stage3_select=0.001 write=3.456 total=69.556
```

`Timing run_search` は検索オーケストレータの内訳です。インデックスボリュームの open / close (`s1_open`、`s2_open`)、Stage 1 の集計 (`s1_compute`)、mode 1 の結果畳み込み (`s1_fold`)、Stage 1 と Stage 2 の in-total (L) 制限 (`s1_intotal`、`s2_intotal`)、Stage 2A と Stage 2B (`s2a`、`s2b`)、ext_job ごとの Stage 2 状態の解放 (`s2_free`)、フラグメント重複の dedup (`dedup`)、parent 単位の上位 N 選択 (`parent_topn`) を表します。`-mode 1` では Stage 1 のキーのみが出力されます。`Timing overall` はその外側の工程、すなわちインデックスの open、クエリ前処理、`run_search` 本体、出力レコードへの変換、Stage 3 アライメント、Stage 3 の各種制限、出力の書き出しを表します。`total` は各キーの和ではなく実測値なので、和との差はどのステージにも帰属されていない時間を示します。オーケストレータは共通なので、`ikafssnserver` でも検索リクエストごとに `Timing run_search` の行が記録されます。

**使用例:**

```bash
# 基本的な検索
ikafssnsearch -ix ./index/mydb -query query.fasta -nthread 8

# k-mer サイズを指定 (インデックスに複数の k 値が含まれる場合は必須)
ikafssnsearch -ix ./index/mydb -k 11 -query query.fasta

# 感度を上げた検索 (クエリごとの候補数を増やす)
ikafssnsearch -ix ./index/mydb -query query.fasta \
    -stage2_min_score 2 -stage1_max_nhit_per_subject 0 -stage1_max_nhit_in_total 2000

# seqidlist で検索対象を限定
ikafssnsearch -ix ./index/mydb -query query.fasta -seqidlist targets.txt

# negative_seqidlist で特定配列を除外
ikafssnsearch -ix ./index/mydb -query query.fasta -negative_seqidlist exclude.txt

# 割合指定の Stage 1 閾値 (クエリ k-mer の 50%)
ikafssnsearch -ix ./index/mydb -query query.fasta -stage1_min_score 0.5

# モード 3: トレースバック付きペアワイズアライメント
ikafssnsearch -ix ./index/mydb -query query.fasta -mode 3 -stage3_traceback 1 -stage3_max_nhit_in_total 5

# モード 3: SAM 出力
ikafssnsearch -ix ./index/mydb -query query.fasta -mode 3 -stage3_traceback 1 -output_format sam -o result.sam

# モード 3: BAM 出力
ikafssnsearch -ix ./index/mydb -query query.fasta -mode 3 -stage3_traceback 1 -output_format bam -o result.bam

# モード 3: 正スコア率でフィルタ
ikafssnsearch -ix ./index/mydb -query query.fasta -mode 3 -stage3_traceback 1 -stage3_min_ppositive 90

# モード 3: コンテクスト拡張 (前後各50塩基)
ikafssnsearch -ix ./index/mydb -query query.fasta -mode 3 -context_extend 50 -stage3_max_nhit_in_total 5

# モード 3: コンテクスト拡張 (前後各クエリ長の0.1倍)
ikafssnsearch -ix ./index/mydb -query query.fasta -mode 3 -context_extend 0.1 -stage3_max_nhit_in_total 5

# パイプラインで ikafssnretrieve に接続
ikafssnsearch -ix ./index/mydb -query query.fasta | ikafssnretrieve -db nt > matches.fasta
```

#### In-Silico PCR（プライマーモード）

`ikafssnsearch` および `ikafssnclient` は `-primer` オプションでプライマーモードを利用できます。このモードでは、プライマーペアから仮想的な PCR アンプリコンクエリを構成し、インデックスに対して検索を行います。`ikafssnclient` ではクライアント側でプライマーからクエリへの変換を行うため、サーバ側の変更は不要です。

**オプション:**

```
  -primer <path>              プライマー FASTA ファイル (-query と排他)
  -insert_length <int>        予想インサート長 (-primer 時必須)
  -stage1_primer_score <num>  プライマーの Stage 1 閾値 (デフォルト: 0.5)
                              0 < v ≤ 1: クエリ k-mer に対する割合
                              v ≥ 2: 絶対値
  -stage2_primer_score_add <int>  Stage 2 プライマースコア加算値 (デフォルト: 1)
```

`-primer` と `-query` は排他的であり、同時に指定できません。`-primer` 使用時は `-insert_length` が必須です。

**プライマーペアの処理方法:**

プライマー FASTA ファイルには偶数本の配列を含める必要があります。連続する 2 配列がペア (奇数番目 = フォワードプライマー、偶数番目 = リバースプライマー) として解釈されます。

各プライマーペアから以下の構成でクエリ配列が生成されます:

```
fwd + N×insert_length + RC(rev)
```

ここで `fwd` はフォワードプライマー配列、`N×insert_length` は指定長の N 文字列 (未知領域)、`RC(rev)` はリバースプライマーの逆相補配列です。

**閾値解決ロジック:**

プライマーモードでは、クエリ全体の k-mer のうちプライマー由来の k-mer のみが実際にマッチに寄与します (N 領域の k-mer は N を含むためスキップされます)。`-stage1_primer_score` は以下のように解決されます:

- `0 < v ≤ 1` (割合): Stage 1 閾値は `ceil(Nprimer_kmer * v)` に解決されます (Nprimer_kmer はプライマー由来の k-mer 数)
- `v ≥ 2` (絶対値): そのまま絶対閾値として使用されます

Stage 2 では、最小チェインスコアが `max(Lf, Lr) + stage2_primer_score_add` に設定されます。ここで Lf と Lr はそれぞれフォワードプライマーとリバースプライマーの k-mer 位置数です。また、`-stage2_max_gap` はコマンドラインで明示的に指定しない限り、自動的にインサート長に設定されます。

**使用例:**

```bash
# 基本的な In-Silico PCR
ikafssnsearch -ix ./index/mydb -primer primers.fasta -insert_length 500

# Stage 1 閾値をプライマー k-mer の 80% に設定
ikafssnsearch -ix ./index/mydb -primer primers.fasta -insert_length 300 \
    -stage1_primer_score 0.8

# モード 3 でアライメントまで実行
ikafssnsearch -ix ./index/mydb -primer primers.fasta -insert_length 500 \
    -mode 3 -stage3_traceback 1 -stage3_max_nhit_in_total 10
```

### ikafssnretrieve

検索結果に基づきマッチした部分配列を抽出します。配列ソースとしてローカル BLAST DB または NCBI E-utilities (efetch) を選択できます。

**FASTA defline:** 各レコードは
`>parent_accession:start-end query=<qseqid> strand=<+|-> score=<chainscore|alnscore>`
形式で出力されます。`start` / `end` は 1-based 包括の **親相対座標** で、左端のアクセッションは常に親 OID のもの（フラグメント由来の合成名は使われません）。検索出力の `sstart` / `send` とは異なり、この範囲は常に昇順です。`-` 鎖の行でも範囲は正鎖番号で同じ区間を指し、出力される配列がその逆相補になります。ローカルパスでは `BlastDbReader::get_subsequence(parent_oid, start, end)` 経由でサブシーケンスを取得するため、染色体級の親でも要求された区間だけがデコードされ、OID 全体を読み込むことはありません。

```
ikafssnretrieve [options]

配列ソース (いずれか必須):
  -db <path>              ローカル BLAST DB プレフィックス
  -remote                 NCBI efetch からリモート取得

入力:
  -tsv <path>         検索結果ファイル (TSV 形式)
  (指定なし)              標準入力から読み込み

共通オプション:
  -o <path>               出力 FASTA ファイル (デフォルト: 標準出力)
                          拡張子 (.gz/.bz2/.xz/.zst) で出力コーデックを選択
  -context_extend <value>        コンテクスト拡張 (デフォルト: 2.0)
                          整数: マッチ領域の前後に付加する塩基数
                          小数: クエリ長 (qlen) に対する倍率
  -compression_level <int> 出力圧縮レベル (既定値: gzip=6, bzip2=9, xz=6, zstd=3)
  -v, --verbose           詳細ログ出力

(注: -tsv は gzip/bzip2/xz/zstd 圧縮された検索結果ファイルを先頭マジック
 バイトで自動判定するため、フラグ指定は不要です。)

リモート取得オプション (-remote 時):
  -api_key <key>          NCBI API key (環境変数 NCBI_API_KEY でも可)
  -batch_size <int>       バッチあたりのアクセッション数 (デフォルト: 100)
  -max_nretry <int>          リトライ回数 (デフォルト: 3)
  -timeout <int>          リクエストタイムアウト秒数 (デフォルト: 30)
  -range_threshold <int>  部分取得に切り替える配列長閾値 (デフォルト: 100000)
```

**使用例:**

```bash
# ローカル BLAST DB から抽出 (ファイル入力)
ikafssnsearch -ix ./index/mydb -query query.fasta -o results.tsv
ikafssnretrieve -db nt -tsv results.tsv -o matches.fasta

# ローカル BLAST DB から抽出 (パイプ)
ikafssnsearch -ix ./index/mydb -query query.fasta | ikafssnretrieve -db nt > matches.fasta

# サーバ経由の検索結果からも同様に利用可能
ikafssnclient -socket /var/run/ikafssn.sock -ix nt -query query.fasta | ikafssnretrieve -db nt > matches.fasta

# NCBI efetch からリモート取得
ikafssnsearch -ix ./index/mydb -query query.fasta | ikafssnretrieve -remote > matches.fasta

# NCBI efetch + API key (高スループット)
ikafssnclient -http http://search.example.com:8080 -ix nt -query query.fasta \
    | ikafssnretrieve -remote -api_key XXXXXXXX > matches.fasta

# マッチ領域の前後 50bp を含めて抽出
ikafssnretrieve -db nt -tsv results.tsv -context_extend 50

# コンテクストをクエリ長の倍率で指定 (前後各 0.1 倍)
ikafssnretrieve -db nt -tsv results.tsv -context_extend 0.1
```

### ikafssnserver

検索デーモンです。インデックスを mmap でメモリに常駐させ、UNIX ドメインソケットまたは TCP ソケットで検索リクエストを受け付けます。

```
ikafssnserver [options]

必須:
  -ix <prefix>            インデックスプレフィックス (繰り返し指定でマルチ DB 対応)

リスニング (いずれかまたは両方):
  -socket <path>          UNIX ドメインソケットパス
  -tcp <host>:<port>      TCP リスニングアドレス

索引変種のロードフィルタ (無指定の項目はワイルドカード):
  -k <int>                この k-mer 長のみロード (デフォルト: 任意)
  -t <int>                このテンプレート長のみロード: 0=連続, 13/15/16/18/21
                          (デフォルト: 任意; -t 0 は連続インデックスのみロード)
  -template_type <str>    coding, optimal, contiguous, both (デフォルト: 任意; -t 0 とは併用不可)
  -min_seq_length <int>   この min_seq_length の変種のみロード (デフォルト: 任意)
  -min_length_split <int> この min_length_split の変種のみロード (デフォルト: 任意)
  -overlap_length <int>   この overlap_length の変種のみロード (デフォルト: 任意)
  -max_freq_build <int>   この max_freq_build の変種のみロード (デフォルト: 任意)
  -max_degen_expand <int> この max_degen_expand の変種のみロード (デフォルト: 任意)。
                          ロードフィルタ専用。query 側の max_degen_expand はクライアントが
                          指定する(リクエストが省略した場合のサーバ内部フォールバックは 16)。

オプション:
  -nthread <int>          ワーカースレッド数 (デフォルト: 利用可能な全コア)
  -max_queue_size <int>   同時処理クエリ配列数のグローバル上限 (デフォルト: 1024)
  -max_nseq_per_req <int> 1 リクエストあたりの受理配列数上限 (デフォルト: スレッド数)
  -max_concurrent_search <int>  -memory_limit の残余 posting_budget を共有する
                          同時検索数の上限。in-flight な posting ヒープ
                          (Stage 3) を抑制 (デフォルト: 0 = 無制限)
  -pid <path>             PID ファイルパス
  -db <path>              モード 3 用 BLAST DB パス (繰り返し指定、-ix と対応;
                          デフォルト: 対応する -ix プレフィックスと同じ)
  -stage1_max_nhit_per_subject <int>  parent ごとの Stage 1 候補数上限のデフォルト (デフォルト: 1、0=無制限)
  -stage1_max_nhit_per_subject_mode <1|2|3|4>  per-subject 選別モードのデフォルト (デフォルト: 3)
  -stage1_max_nhit_per_volume <int>  (query,volume,strand) ごとの Stage 1 候補数上限のデフォルト (デフォルト: 0)
  -stage1_max_nhit_in_total <int>  全 volume を通したクエリごとの Stage 1 候補数上限のデフォルト (デフォルト: 0)
  -stage1_min_score <num> デフォルト Stage 1 最小スコア閾値 (デフォルト: 0.5)
                          整数 (>= 1) または小数 (0 < P < 1)
  -stage2_min_score <int> デフォルト最小チェインスコア (デフォルト: 0 = 適応的)
  -stage2_max_gap <int>   デフォルトチェイニング対角線ずれ許容幅 (デフォルト: 100)
  -stage2_max_lookback <int>  デフォルトチェイニング DP 探索窓サイズ (デフォルト: 64、0=無制限)
  -stage2_max_nhit_per_subject <int>  デフォルトサブジェクトあたりの最大チェイン数 (デフォルト: 1、0=無制限)
  -stage2_max_nhit_per_subject_mode <1|2|3|4>  デフォルトのサブジェクト単位選択モード (デフォルト: 3)
                           1/2=上位N位 (タイ非包含)、3/4=上位N位+スコアタイ包含;
                           1/3=parent ごとに strand 統合、2/4=strand 分離
  -stage2_max_nhit_in_total <int>  デフォルトの全ボリューム通算 Stage 2 最大チェイン数 (デフォルト: 0)
  -stage2_min_nhit_diag <int> デフォルト対角線フィルタ最小ヒット数 (デフォルト: 1)
  -context_extend <value>        デフォルトコンテクスト拡張 (デフォルト: 2.0)
                          整数: 拡張する塩基数; 小数: クエリ長に対する倍率
  -stage3_traceback <0|1> デフォルトトレースバックモード (デフォルト: 0)
  -stage3_gapopen <int>   デフォルトギャップオープンペナルティ (デフォルト: 10)
  -stage3_gapext <int>    デフォルトギャップ伸長ペナルティ (デフォルト: 1)
  -stage3_min_ppositive <num>  デフォルト最小正スコア率 (デフォルト: 0)
  -stage3_min_npositive <int>  デフォルト最小正スコア塩基数 (デフォルト: 0)
  -stage3_max_nhit_per_subject <int>  デフォルトの Stage 3 サブジェクトあたり最大ヒット数 (デフォルト: 1、0=無制限)
  -stage3_max_nhit_per_subject_mode <1|2|3|4>  デフォルトのサブジェクト単位選択モード (デフォルト: 3)
  -stage3_max_nhit_in_total <int>  デフォルトの全ボリューム通算 Stage 3 最大ヒット数 (デフォルト: 0)
  -stage3_score_matrix <str>  デフォルトスコア行列 (デフォルト: degmatch)
                          利用可能: degmatch、dnafull、nuc44
  -accept_qdegen <0|1>    デフォルト縮重塩基クエリ許可 (デフォルト: 1)
  -memory_limit <size>    検索メモリ予算 (デフォルト: 物理メモリの半分)
                          .khx/.ksx メタデータを常駐させ、残余で Stage 3 の
                          posting heap と同時検索プールを制限。接尾辞 K, M, G を認識
  -shutdown_timeout <int> グレースフルシャットダウンのタイムアウト秒数 (デフォルト: 30)
  -v, --verbose           詳細ログ出力
```

ロードフィルタはサーバがロードする索引変種を制限し、`ikafssninfo` ローカルモードと同じワイルドカード意味論を用います(無指定の項目はワイルドカード、無指定の `-t` は完全ワイルドカード、`-t 0` と `-template_type` の併用はエラー)。**`-max_degen_expand` を含む 8 識別パラメータすべて**がここではロードフィルタです。`-max_degen_expand` はサーバでは**ロードフィルタ専用**で、query 側の `max_degen_expand`(k-mer 展開+タイブレーク)はクライアントが指定し、サーバはリクエストが省略した場合のみ内部フォールバック 16 を適用します。`-max_degen_expand` 無指定なら構築時の全 maxexpand 値をロードし(リクエスト時タイブレークが選択可能)、指定するとロード対象が絞られ、クライアントはロードされた変種からしか選択できなくなります。

サーバはロードした各変種をワイヤ上で別々のグループとして公開します(info 応答は各グループに 8 つの識別パラメータすべてを載せます)。検索リクエストはクライアントが解決した変種識別を載せ、サーバは `t` と `template_type` をリテラル一致(0=連続)、4 つの索引フィールドを厳密一致させ、`max_degen_expand` 違いの変種が複数残った場合はリクエストの `max_degen_expand` に一致するもの(なければ最大)を選択します。`template_type=both` のリクエストでは、`max_degen_expand` を共有する coding+optimal ペアを選択します(混在不可)。`-min_query_length` と長尺クエリ上限は、選択された変種の `min_seq_length` と `overlap_length` からリクエストごとに導出されます。

**使用例:**

```bash
# UNIX ソケットで起動
ikafssnserver -ix ./nt_index -socket /var/run/ikafssn.sock -nthread 16

# TCP で起動 (リモートアクセス用)
ikafssnserver -ix ./nt_index -tcp 0.0.0.0:9100 -nthread 32

# 両方で同時リスニング
ikafssnserver -ix ./nt_index -socket /var/run/ikafssn.sock -tcp 0.0.0.0:9100

# デーモン運用 (systemd 等から起動)
ikafssnserver -ix ./nt_index -socket /var/run/ikafssn.sock -pid /var/run/ikafssn.pid

# モード 3 対応: BLAST DB と Stage 3 デフォルトを指定
ikafssnserver -ix ./nt_index -db nt -socket /var/run/ikafssn.sock \
    -stage3_traceback 1 -context_extend 50

# マルチ DB: 1 プロセスで 2 つのデータベースを同時サーブ
ikafssnserver -ix ./nt_index -db nt -ix ./rs_index -db refseq_genomic \
    -socket /var/run/ikafssn.sock -nthread 32
```

**運用上の特性:**

- 1 プロセスで複数の BLAST DB インデックスを同時にサーブできます。`-ix` (および必要に応じて `-db`) を複数回指定して複数データベースをロードします。各データベースは `-ix` プレフィックスのベースネーム (パスの最終コンポーネント) で識別され、サーバが複数 DB をホストする場合、クライアントは `-db <name>` でターゲット DB を指定する必要があります。
- `-db` を指定する場合、その数は `-ix` の数と一致する必要があります (順番に対応)。`-db` を省略した DB は `-ix` プレフィックスを BLAST DB パスとして使用します。`-db` パスが未指定の DB はモード 1-2 のみ対応 (max_mode=2)、`-db` を指定するとモード 3 も利用可能 (max_mode=3) になります。
- `-ix` プレフィックスに対応する索引変種が複数存在する場合 (8 つの識別パラメータ — k, t, template_type, min_seq_length, min_length_split, overlap_length, max_freq_build, max_degen_expand — のいずれかが異なる)、ロードフィルタに従って全てを読み込み、各クライアントリクエストが解決済みの変種識別を載せて 1 つを選択します。どの変種をロードするか自体はロードフィルタで制限できます。
- SIGTERM/SIGINT 受信時はグレースフルシャットダウンを行います。新規接続の受付を停止し、実行中のリクエストの完了を最大 `-shutdown_timeout` 秒待ちます。
- **配列単位の同時実行制御:** サーバは接続単位ではなく、配列単位で同時実行数を制御します。リクエストが到着すると、有効なクエリ配列ごとにパーミットの取得を試みます。グローバル上限 (`-max_queue_size`) に達した場合、超過分の配列はリトライ用に「拒否」としてクライアントに返されます。`-max_nseq_per_req` は 1 リクエストが取得できるパーミット数の上限を設定し、大量配列を含む単一リクエストによるスロットの独占を防ぎます。
- **リクエスト間の posting budget プール (`-max_concurrent_search`):** デフォルト値 `0` では、同時実行中の各リクエストが独立に残余 `posting_budget` (`-memory_limit` から `.khx` / `.ksx` の `WILLNEED` を差し引いた残り) の全量を使用できます。`N >= 1` を指定すると、リクエスト間でこの予算をリース/リリース方式で共有します。各リクエストは検索全体を通じて 1 個のリースを保持し、同時に実行できる検索は最大 `N` 個に制限されます。これにより接続数によらず in-flight な posting ヒープ (Stage 3 のアライメントバッチャが消費) を `posting_budget` 以下に抑えられますが、`N` を超える追加検索はリースが解放されるまで待機します。1 リクエスト分にメモリ上限を合わせるのは不便だが、プロセス当たりの posting メモリには厳格な上限を設けたい運用に向きます。

### ikafssnhttpd

HTTP REST API デーモンです。1 つ以上の `ikafssnserver` インスタンスに接続し、SQLite ジョブストアとワーカープールを介した非同期ポーリング型 REST API として検索サービスを提供します。Drogon フレームワークを使用します。複数のバックエンドを指定することで、マルチデータベース対応や同一 DB のレプリカ間の負荷分散が可能です。

起動時にすべてのバックエンドに接続してサーバ能力情報をキャッシュします (最大 30 秒間、指数バックオフでリトライ)。同じデータベース名が複数のバックエンドに存在する場合、共通する (db, k) の組み合わせについて合計配列数・合計塩基数の一致をクロスサーバ検証し、不一致があれば起動エラーとなります。バックエンド間で k 値セットが異なることは許容されます (例: サーバ A が k=10、サーバ B が k=10 と k=11 を持つ構成)。統合された能力情報には全バックエンドの k 値グループの和集合が含まれます。検索リクエストはまず統合された能力情報に対して同期的にチェック (バックエンドへの通信なし) し、明らかに無効なリクエストを即座に拒否します。検証通過後はジョブストアに永続化し、ワーカープールが優先度と空きスロット状況に基づいて最適なバックエンドへ非同期にディスパッチします。

**ルーティングとヘルスチェック:**

- **優先度**: バックエンドは CLI 引数の順序で優先度が決まります (先 = 最高優先度)。
- **選択**: 検索リクエストごとに、最高優先度かつ有効な空き容量があるバックエンドを選択します。有効な空き容量は、スロット空き数 (`max_queue_size - queue_depth`) とリクエストあたり上限 (`max_nseq_per_req`) の小さい方で決まります。全バックエンドが満杯の場合は、最高優先度のバックエンドが使用されます。
- **事前チェック**: 検索前に選択されたバックエンドへ最新の info リクエストを送信し、接続を確認します。
- **除外**: バックエンドが応答しない場合 (info/search の接続エラー)、`-exclusion_time` 秒間除外されます。除外されたバックエンドはハートビート時に自動的に再チェックされ、到達可能になれば再有効化されます。
- **ハートビート**: バックグラウンドスレッドが `-heartbeat_interval` 秒ごとにすべてのバックエンドの info を更新します。
- **httpd 側でのリトライ**: ディスパッチ失敗時はワーカープールが `-max_nretry` 回まで再試行し、それでも失敗したジョブのみ `failed` としてマークします。バックエンドから返される個別クエリの拒否もワーカーがリトライします。

```
ikafssnhttpd [options]

バックエンド接続 (1 つ以上必須、指定順 = 優先度):
  -server_socket <path>      ikafssnserver の UNIX ソケットパス
  -server_tcp <host>:<port>  ikafssnserver の TCP アドレス

ジョブストア (非同期 REST):
  -sqlite_db <path>           SQLite ジョブストアのパス
                              (デフォルト: /var/lib/ikafssnhttpd/jobs.db)
  -query_dir <path>           ジョブクエリファイル保存ディレクトリ。queued/running
                              ジョブのリクエストボディを 1 ジョブ 1 ファイル
                              (zstd 圧縮) として保持し、ジョブが終局状態へ
                              遷移した時点で削除します
                              (デフォルト: /var/lib/ikafssnhttpd/queries)
  -query_compression_level <int>  プレーンテキスト JSON
                              (Content-Type: application/json) アップロード時に
                              適用する zstd 圧縮レベル。Content-Type:
                              application/zstd で送られたボディはクライアント
                              選択の圧縮レベルを保ったままそのまま保存され、
                              本オプションの影響を受けません
                              (デフォルト: 3、許容範囲: 1-22)
  -result_dir <path>          ジョブ結果ファイル保存ディレクトリ。完了ジョブごとに
                              1 つの zstd ファイルを書き、sendfile(2) でクライアントへ
                              直接送出します
                              (デフォルト: /var/lib/ikafssnhttpd/results)
  -result_compression_level <int>  結果ファイルの zstd 圧縮レベル
                              (デフォルト: 3、許容範囲: 1-22)
  -retention_time <int>       完了/失敗ジョブの保持秒数 (デフォルト: 86400)
                              これを過ぎるとハウスキーパが該当行と対応する
                              結果ファイルを削除します
  -max_nretry <int>           1 ジョブあたりのディスパッチ再試行上限 (デフォルト: 3)
  -nthread_worker <int>       バックエンド検索ワーカープールのスレッド数 (デフォルト: 4)

オプション:
  -listen <host>:<port>       HTTP リスニングアドレス (デフォルト: 0.0.0.0:8080)
  -path_prefix <prefix>       API パスプレフィックス (例: /nt)
  -nthread <int>              Drogon I/O スレッド数 (デフォルト: 利用可能な全コア)
  -heartbeat_interval <int>   ハートビート間隔 (秒、デフォルト: 3600)
  -exclusion_time <int>       バックエンド除外時間 (秒、デフォルト: 3600)
  -pid <path>                 PID ファイルパス
  -v, --verbose               詳細ログ出力
```

**REST API エンドポイント (非同期ポーリング型):**

| メソッド | パス | 説明 |
|---|---|---|
| POST | `/api/v1/jobs` | 検索ジョブの投入。ボディは JSON ドキュメント。`ikafssnclient` は常に単一 zstd フレームに圧縮し `Content-Type: application/zstd` で送信します。curl やスクリプトなど一時的な HTTP クライアントから直接投入する場合は、事前圧縮を省略するために `Content-Type: application/json` で生 JSON を送ることもできます。それ以外の Content-Type は HTTP 415 で拒否されます。`{job_id, status}` を即時返却 |
| GET  | `/api/v1/jobs/{job_id}` | ジョブステータス (`queued` / `running` / `done` / `failed` / `timeout`) |
| GET  | `/api/v1/jobs/{job_id}/result` | シリアライズ済み SearchResponse blob (status = `done` 時のみ有効)。レスポンスボディは単一 zstd フレームで `Content-Type: application/zstd`。デシリアライズ前に `zstd -d` または `ZSTD_decompress` で展開してください |
| GET  | `/api/v1/info` | 全バックエンドの統合インデックス情報 |
| GET  | `/api/v1/health` | ヘルスチェック (いずれかのバックエンドが到達可能なら OK) |

`/api/v1/info` レスポンスは全 healthy バックエンドのデータベースを統合して返します。複数バックエンドで提供されるデータベースの場合、容量情報は kmer_group 内の `modes` 配列にモードごとに集約され、全提供バックエンドの `max_queue_size`、`queue_depth`、および `max_nseq_per_req` (各バックエンドの `min(空きスロット, per_req)` の合計として算出) が表示されます。トップレベルにも全モードの最小値として `max_nseq_per_req` フィールドが出力されます。

`POST /api/v1/jobs` はリクエストボディを検証した後、まず `<query_dir>/<ab>/<job_id>.json.zst` (ここで `<ab>` は `job_id` 先頭 2 文字。シャーディングによりディレクトリエントリ数を抑える) に zstd ファイルとして書き出し、その後 SQLite 行を INSERT します。`Content-Type: application/zstd` でアップロードされたボディはクライアントが選んだ圧縮レベルを保ったまま検証後そのまま保存され (再圧縮は行いません)、`Content-Type: application/json` のプレーンテキスト JSON は `-query_compression_level` で圧縮してから保存します。その後ワーカープールがバックエンドへ非同期にディスパッチします。ジョブが正常終了するとワーカーはシリアライズ済み `SearchResponse` を `<result_dir>/<ab>/<job_id>.bin.zst` に書き出し、SQLite 行のステータスを `done` に更新します。`GET /api/v1/jobs/{job_id}/result` はそのファイルを `sendfile(2)` でゼロコピー送出 (`Content-Type: application/zstd`) するため、`ikafssnhttpd` は結果データをメモリに載せません。同じタイミングで対応するクエリファイルも削除されるため、`-query_dir` 配下に残るのは queued / running 中のジョブのみとなり、滞留が `jobs.db` の容量へ波及しません (新スキーマで `jobs.db` はメタデータのみ)。クライアントは `GET /api/v1/jobs/{job_id}` を `done` になるまでポーリングし、結果を取得します。完了 (成功・失敗) ジョブは `-retention_time` 秒間保持され、ハウスキーパが期限切れ行を SQLite から削除した後、対応する結果ファイルおよび残存クエリファイル (保険) を順に削除します。起動時には強制終了で書きかけ状態だった `*.bin.zst.tmp` および `*.json.zst.tmp` を一括スイープし、また定期的なオーファンスイープで SQLite 行が消えた残存結果ファイル・クエリファイル (まれ: `mark_done` 失敗後、または `mark_done` とクエリ unlink の間でのワーカークラッシュ後のみ発生) も回収します。

**使用例:**

```bash
# 単一バックエンド: UNIX ソケット
ikafssnhttpd -server_socket /var/run/ikafssn.sock -listen 0.0.0.0:8080

# 単一バックエンド: TCP
ikafssnhttpd -server_tcp 10.0.1.5:9100 -listen 0.0.0.0:8080

# 複数バックエンドで負荷分散 (同じ DB を 2 台のサーバで提供)
ikafssnhttpd -server_tcp server1:9100 -server_tcp server2:9100 -listen :8080

# 異なる DB を持つ複数バックエンド
ikafssnhttpd -server_socket /var/run/nt.sock -server_socket /var/run/rs.sock -listen :8080

# ミックス構成: プライマリ + フェイルオーバー
ikafssnhttpd -server_socket /var/run/primary.sock -server_tcp backup:9100 -listen :8080
```

### ikafssnclient

クライアントコマンドです。`ikafssnserver` に直接ソケット接続するか、`ikafssnhttpd` に HTTP 接続して検索結果を取得します。出力形式は `ikafssnsearch` と同一です。クエリ送信前にサーバの能力情報を取得し、指定されたデータベース名・k-mer サイズ・モードの妥当性を事前検証 (プリフライトチェック) します。続いて `ikafssnsearch` と全く同じ手順で索引変種を解決します — サーバの変種別グループ一覧に同じ `-k` / `-t` / `-template_type` / `-min_seq_length` / `-min_length_split` / `-overlap_length` / `-max_freq_build` フィルタと同じ template_type 解決(単一・both 必須・ワイルドカードで coding/optimal フォールバック。both ペアは `max_degen_expand` を共有し混在しない)を適用し、解決した識別をリクエストに載せて送信します。`ikafssnsearch` と同様、`max_degen_expand` 以外の 7 パラメータで単一の変種(または 1 つの both ペア)に収束しなければならず、`max_degen_expand` 自体はサーバのタイブレークに委ねます。曖昧または該当なしの場合は、クエリデータの送信前にローカルでエラーになります。無効なパラメータが指定された場合は、利用可能なデータベース一覧を含むエラーメッセージが表示されます。クライアントはサーバの `max_nseq_per_req` と空きスロット数に基づいてクエリを適切なサイズのバッチに自動分割し、部分的に拒否されるような過大なリクエストを回避します。各バッチ内でサーバが同時実行制限によりクエリ配列を拒否した場合、拒否された配列を指数バックオフ (30 秒、60 秒、120 秒、120 秒、…) で自動リトライし、全配列の処理が完了するまで繰り返します。

**ソケット/TCP モード (同期) — チェックポインティング:** `-socket` または `-tcp` で接続した場合、クライアントはバッチ処理中の中間結果を一時ディレクトリに自動保存します。プロセスが中断された場合 (例: Ctrl+C、ネットワーク障害)、同じコマンドを再実行すると中断箇所から再開し、処理済みクエリをスキップします。一時ディレクトリの命名は `{出力}.{入力}.{ix名}.{kk}.ikafssn.tmp/` で、正常完了後に自動削除されます。ディレクトリベースのロックにより同一パラメータでの同時実行を防止します。再開時の検証では検索パラメータ、入力ファイルの SHA256、各バッチファイルの整合性をチェックします。

**HTTP モード (非同期ポーリング型):** `-http` で接続した場合、クライアントは入力をバッチ単位の複数ジョブに分割し、`POST /api/v1/jobs` で `ikafssnhttpd` に投入します。リクエストボディは送信前に zstd 圧縮され (`Content-Type: application/zstd`)、巨大なクエリバッチでも転送帯域を抑えられます (サーバ側で透過的に展開)。各ジョブのメタデータは `~/.ikafssnclient/<group_id>/<job_id>.json` に永続化され、defline キャッシュ (`<job_id>.deflines.zst`) も併せて保存されます。クライアントは `GET /api/v1/jobs/{job_id}` を終端ステータスに達するまでポーリングし、結果を取得します。サーバから返るレスポンスボディは既に単一 zstd フレームなので、`<job_id>.result.bin.zst` にはそのままバイト列として書き込み (二重圧縮しない)、続いて展開・デシリアライズしてユーザの出力ファイルへマージします。状態がすべてディスク上に保持されているため、再起動を跨いでも復帰可能です — `-resume <id>` は既存のグループ/ジョブのポーリングループに再入し、`-submit_only` は投入直後に `group_id` を返して終了するため、別プロセス (または後続の `-resume` 呼び出し) で結果を回収できます。注意: 旧クライアントが書いた `*.result.bin.zst` は zstd を二重に掛けた旧フォーマットであり、新コードとは互換性がありません。アップグレード時は未完了のグループは破棄してください。

```
ikafssnclient [options]

接続先 (いずれか):
  -socket <path>           ikafssnserver の UNIX ソケットパス
  -tcp <host>:<port>       ikafssnserver の TCP アドレス
  -http <url>              ikafssnhttpd の URL (例: http://example.com:8080)
                           HTTP モードは非同期 REST ポーリングを使用

非同期 REST ジョブ管理 (HTTP モード専用、検索とは排他):
  -submit_only             バッチを投入し、group_id を表示して終了 (ポーリングしない)
                           回収は後で -resume を使用
  -jobs                    ~/.ikafssnclient/ 配下のローカル管理ジョブグループを一覧
                           (サーバとは通信しない)
  -job_detail <id>          グループ内のジョブ一覧、または単一ジョブの詳細を表示
  -resume <id>             既存グループ/単一ジョブのポーリングを再開

必須 (新規検索時):
  -query <path>            クエリ FASTA ファイル (- で標準入力)
  -ix <name>               サーバ上のターゲットデータベース名

プライマーモード (-query の代替):
  -primer <path>           プライマーペア FASTA ファイル (-query と排他)
  -insert_length <int>     予想インサート長 (-primer 時必須)
  -stage1_primer_score <num>  Stage 1 閾値 (0<v≤1: 割合、v≥2: 絶対値; デフォルト: 0.5)
  -stage2_primer_score_add <int>  Stage 2 閾値加算値: max(Lf,Lr) + N (デフォルト: 1)

オプション:
  -o <path>                出力ファイル (デフォルト: 標準出力)
  -k <int>                 使用する k-mer サイズ (デフォルト: サーバ側デフォルト)
  -mode <1|2|3>            検索モード (デフォルト: サーバ側デフォルト)
  -stage1_max_nhit_per_subject <int>  parent ごとの Stage 1 候補数上限 (デフォルト: サーバ側デフォルト)
  -stage1_max_nhit_per_subject_mode <1|2|3|4>  per-subject 選別モード (デフォルト: サーバ側デフォルト)
  -stage1_max_nhit_per_volume <int>  (query,volume,strand) ごとの Stage 1 候補数上限 (デフォルト: サーバ側デフォルト)
  -stage1_max_nhit_in_total <int>  全 volume を通したクエリごとの Stage 1 候補数上限 (デフォルト: サーバ側デフォルト)
  -stage1_min_score <num>  Stage 1 最小スコア閾値 (デフォルト: サーバ側デフォルト)
                           整数 (>= 1) または小数 (0 < P < 1)
  -stage2_min_score <int>  最小チェインスコア (デフォルト: サーバ側デフォルト)
                           0 = 適応的モードを明示的にリクエスト
  -stage2_max_gap <int>    チェイニング対角線ずれ許容幅 (デフォルト: サーバ側デフォルト)
  -stage2_max_lookback <int>  チェイニング DP 探索窓サイズ (デフォルト: サーバ側デフォルト)
  -stage2_max_nhit_per_subject <int>  サブジェクトあたりの最大チェイン数 (デフォルト: サーバ側デフォルト)
  -stage2_max_nhit_per_subject_mode <1|2|3|4>  サブジェクト単位の選択モード (デフォルト: 3)
                           1/2=上位N位 (タイ非包含)、3/4=上位N位+スコアタイ包含;
                           1/3=parent ごとに strand 統合、2/4=strand 分離
  -stage2_max_nhit_in_total <int>  全ボリュームを通したクエリあたりの最大チェイン数 (-mode 2 以上が必要)
  -stage2_min_nhit_diag <int> 対角線フィルタ最小ヒット数 (デフォルト: サーバ側デフォルト)
  -context_extend <value>         コンテクスト拡張 (デフォルト: 2.0)
                           整数: 拡張する塩基数; 小数: クエリ長に対する倍率
  -stage3_traceback <0|1>  トレースバック有効化 (デフォルト: サーバ側デフォルト)
  -stage3_gapopen <int>    ギャップオープンペナルティ (デフォルト: サーバ側デフォルト)
  -stage3_gapext <int>     ギャップ伸長ペナルティ (デフォルト: サーバ側デフォルト)
  -stage3_min_ppositive <num> 最小正スコア率フィルタ (デフォルト: サーバ側デフォルト)
  -stage3_min_npositive <int> 最小正スコア塩基数フィルタ (デフォルト: サーバ側デフォルト)
  -stage3_max_nhit_per_subject <int>  サブジェクトあたりの最大ヒット数 (-mode 3 が必要)
  -stage3_max_nhit_per_subject_mode <1|2|3|4>  サブジェクト単位の選択モード (デフォルト: 3; -mode 3 が必要)
  -stage3_max_nhit_in_total <int>  全ボリュームを通したクエリあたりの最大ヒット数 (-mode 3 が必要)
  -stage3_score_matrix <str> スコア行列 (デフォルト: サーバ側デフォルト)
                           利用可能: degmatch、dnafull、nuc44
  -seqidlist <path>        検索対象を指定アクセッションに限定
  -negative_seqidlist <path>  指定アクセッションを検索対象から除外
  -strand <-1|1|2>         検索する鎖: 1=プラス、-1=マイナス、2=両鎖 (デフォルト: サーバ側デフォルト)
  -accept_qdegen <0|1>     縮重塩基を含むクエリを許可 (デフォルト: 1)
  -max_degen_expand <int>  縮重塩基展開の最大数 (デフォルト: 16、最大: 256、0/1: 無効)。
                           ikafssnsearch と同一の query パラメータ(サーバに委譲しない)。
                           max_degen_expand 変種が複数残った場合のタイブレークにも使われる。
  -min_query_length <int>  クエリ配列の最短長 (デフォルト: 64)。クライアント
                           側で短いクエリを事前にフィルタするため、サーバ側
                           の同じチェックが発火する前に弾かれます。値は
                           InfoResponse 経由でサーバから取得した
                           min_seq_length 以上である必要があり、それより
                           小さい場合は pre-flight 検証で拒否されます。
  -t <int>                 スペースドシード用テンプレート長 (デフォルト: 0 = 連続、ワイルドカードではない)
                           0: 連続 k-mer; 13, 15, 18 (k=8-9); 16, 18, 21 (k=11-12)
  -template_type <str>     スペースドシードのテンプレート種別 (デフォルト: -t>0 で both; -t 0 とは併用不可)
                           coding: coding インデックスのみ使用
                           optimal: optimal インデックスのみ使用
                           both: coding と optimal のインデックスを検索時にマージ
  -min_seq_length <int>    この min_seq_length の索引変種を選択 (デフォルト: 任意)
  -min_length_split <int>  この min_length_split の索引変種を選択 (デフォルト: 任意)
  -overlap_length <int>    この overlap_length の索引変種を選択 (デフォルト: 任意)
  -max_freq_build <int>    この max_freq_build の索引変種を選択 (デフォルト: 任意)
  -output_format <tsv|json|sam|bam>  出力形式 (デフォルト: tsv)
  -compression_level <int> 出力圧縮レベル (既定値: gzip=6, bzip2=9, xz=6, zstd=3)
                           コーデックは -o の拡張子 (.gz/.bz2/.xz/.zst) で選択。SAM/BAM では拒否
  -v, --verbose            詳細ログ出力

(注: -query および -primer は gzip/bzip2/xz/zstd 圧縮された FASTA を
 先頭マジックバイトで自動判定するため、フラグ指定は不要です。)

HTTP 認証 (HTTP モード専用):
  -user <user:password>   認証情報 (curl 形式)
  -http_user <USER>       ユーザー名 (wget 形式)
  -http_password <PASS>   パスワード (-http_user と併用)
  -netrc_file <path>      .netrc ファイルのパス
```

**使用例:**

```bash
# UNIX ソケット経由 (ローカル)
ikafssnclient -socket /var/run/ikafssn.sock -ix nt -query query.fasta

# TCP 直接接続
ikafssnclient -tcp 10.0.1.5:9100 -ix nt -query query.fasta

# HTTP 経由
ikafssnclient -http http://search.example.com:8080 -ix nt -query query.fasta

# seqidlist で検索対象を限定
ikafssnclient -socket /var/run/ikafssn.sock -ix nt -query query.fasta -seqidlist targets.txt

# パイプラインで ikafssnretrieve に接続
ikafssnclient -socket /var/run/ikafssn.sock -ix nt -query query.fasta | ikafssnretrieve -db nt > matches.fasta

# 特定の k-mer サイズを指定
ikafssnclient -socket /var/run/ikafssn.sock -ix nt -query query.fasta -k 9

# HTTP Basic 認証 (curl 形式)
ikafssnclient -http http://search.example.com:8080 -ix nt -query query.fasta -user admin:secret

# HTTP Basic 認証 (wget 形式)
ikafssnclient -http http://search.example.com:8080 -ix nt -query query.fasta -http_user=admin -http_password=secret

# .netrc ファイルによる認証
ikafssnclient -http http://search.example.com:8080 -ix nt -query query.fasta -netrc_file ~/.netrc

# モード 3: トレースバック付きペアワイズアライメント
ikafssnclient -socket /var/run/ikafssn.sock -ix nt -query query.fasta -mode 3 -stage3_traceback 1

# モード 3: SAM 出力
ikafssnclient -socket /var/run/ikafssn.sock -ix nt -query query.fasta -mode 3 -stage3_traceback 1 -output_format sam -o result.sam

# In-Silico PCR (プライマーモード)
ikafssnclient -socket /var/run/ikafssn.sock -ix nt -primer primers.fasta -insert_length 500

# プライマーモードでカスタム閾値を指定
ikafssnclient -tcp 10.0.1.5:9100 -ix nt -primer primers.fasta -insert_length 300 \
    -stage1_primer_score 0.8

# プライマーモードでモード 3 アライメント
ikafssnclient -socket /var/run/ikafssn.sock -ix nt -primer primers.fasta -insert_length 500 \
    -mode 3 -stage3_traceback 1 -stage3_max_nhit_in_total 10
```

### ikafssninfo

インデックス / データベース情報表示コマンドです。**ローカルモード** (インデックスファイルを直接読み取り) と**リモートモード** (稼働中の `ikafssnserver` または `ikafssnhttpd` に問い合わせ) の 2 つのモードをサポートします。

```
ikafssninfo [options]

必須 (いずれか):
  -ix <prefix>             インデックスプレフィックス [ローカルモード]
  -socket <path>           ikafssnserver の UNIX ソケットパス [リモートモード]
  -tcp <host>:<port>       ikafssnserver の TCP アドレス [リモートモード]
  -http <url>              ikafssnhttpd の URL [リモートモード]

ローカルモードオプション (索引変種フィルタ; 無指定の項目はワイルドカード):
  -db <path>               BLAST DB プレフィックス (デフォルト: -ix から自動検出)
  -k <int>                 k-mer 長 (デフォルト: 任意)
  -t <int>                 テンプレート長: 0=連続, 13/15/16/18/21
                           (デフォルト: 任意; -t 0 は連続インデックスのみ選択)
  -template_type <str>     coding, optimal, contiguous, both (デフォルト: 任意; -t 0 とは併用不可)
  -min_seq_length <int>    min_seq_length でフィルタ (デフォルト: 任意)
  -min_length_split <int>  min_length_split でフィルタ (デフォルト: 任意)
  -overlap_length <int>    overlap_length でフィルタ (デフォルト: 任意)
  -max_freq_build <int>    max_freq_build でフィルタ (デフォルト: 任意)
  -max_degen_expand <int>  max_degen_expand でフィルタ (デフォルト: 任意)
  -stats <0|1>             k-mer 出現頻度分布を計算 (デフォルト: 0; 低速)

リモート HTTP 認証:
  -user <user:password>   認証情報 (curl 形式)
  -http_user <USER>       ユーザー名 (wget 形式)
  -http_password <PASS>   パスワード (-http_user と併用)
  -netrc_file <path>      .netrc ファイルのパス

オプション:
  -v, --verbose            詳細出力
```

`-ix` とリモートオプション (`-socket`、`-tcp`、`-http`) は排他的です。リモートオプションは同時に 1 つのみ指定可能です。

**ローカルモード:** インデックスファイルを直接読み取り、詳細な統計情報を表示します。`-db` 未指定の場合、インデックスプレフィックスパスが有効な BLAST DB に対応するかを確認し、自動検出を試みます。`ikafssnsearch` / `ikafssnclient` と異なり、ローカルモードは単一変種への解決を行いません。フィルタを通過した**すべて**の索引変種を、それぞれ `=== Index variant: … ===` の見出しの下に報告し、統計を変種間で合算することはありません。フィルタオプションはワイルドカード意味論を用います — 無指定の `-t` は任意のテンプレート長に一致し (`ikafssnsearch` では無指定の `-t` は連続を意味する)、`-max_degen_expand` はここでは通常のフィルタ (タイブレークではない) であり、`-t 0` と `-template_type` の併用はエラーになります。

ローカルモードの出力情報:

- k-mer 長 (k) および k-mer 整数型 (uint16/uint32)
- ボリューム数
- 各ボリュームの親 (BLAST OID) 数、フラグメント (内部 SeqId) 数、フラグメント長分布 (min / median / mean / max)、総ポスティング数、ファイルサイズ、除外 k-mer 数 (`.khx` 存在時)
- 全体統計: 総親数、総フラグメント数、集約フラグメント長分布、総ポスティング数、総インデックスサイズ、圧縮率
- ファイル名から parse した変種の 8 識別パラメータ (`k`, `t`, `template_type`, `min_seq_length`, `min_length_split`, `overlap_length`, `max_freq_build`, `max_degen_expand`)
- `-stats 1` 指定時: k-mer 出現頻度分布 (min, max, mean, パーセンタイル)
- `-db` 指定時 (または自動検出時): BLAST DB のタイトル、配列数、総塩基数、ボリューム構成

「Parents」と「Fragments」の値が一致しない場合はインデックスがフラグメント分割を有効にしている (`min_length_split > 0`) ことを意味し、一致する場合は「親 1 つにつきフラグメント 1 つ」の退化レイアウトです。

**リモートモード:** 稼働中のサーバに問い合わせ、サーバの能力情報を表示します。

リモートモードの出力情報:

- キュー深度 / 最大キューサイズ (queue_depth / max_queue_size)
- データベースごとの情報: 名前、デフォルト k、最大モード、および索引変種ごとに 1 グループ。各グループは 8 つの識別パラメータ (`k`, `t`, `template_type`, `min_seq_length`, `min_length_split`, `overlap_length`, `max_freq_build`, `max_degen_expand`) と、ボリューム数・配列数・総塩基数・ポスティング数の統計を表示
- `-v` 指定時: 各変種グループ内のボリュームごとの詳細 (配列数、総塩基数、ポスティング数)

**使用例:**

```bash
# ローカル: 基本的なインデックス情報
ikafssninfo -ix ./index/mydb

# ローカル: BLAST DB 情報も表示
ikafssninfo -ix ./index/mydb -db mydb

# ローカル: 詳細な頻度分布を表示
ikafssninfo -ix ./index/mydb -stats 1

# リモート: UNIX ソケット経由でサーバに問い合わせ
ikafssninfo -socket /var/run/ikafssn.sock

# リモート: TCP 経由でサーバに問い合わせ
ikafssninfo -tcp 10.0.1.5:9100

# リモート: HTTP 経由でサーバに問い合わせ
ikafssninfo -http http://search.example.com:8080

# リモート: 詳細出力 (ボリュームごとの情報を表示)
ikafssninfo -socket /var/run/ikafssn.sock -v

# リモート: HTTP 認証付き
ikafssninfo -http http://search.example.com:8080 -user admin:secret
```

## 検索パイプライン

ikafssn は 3 段階の検索パイプラインを使用します。

デフォルトパラメータはスループットを優先しています。`stage1_min_score=0.5` (割合指定) でクエリ k-mer の 50% 以上のマッチを要求してフィルタリングします。各ステージは同じ候補数制限のファミリ — per-subject (N)、per-volume (M、Stage 1 のみ)、in-total (L) — を備えており、選択した `-mode` の終端ステージの L で最終出力件数を制限できます (mode 1 → Stage 1 L、mode 2 → Stage 2 L、mode 3 → Stage 3 L)。Stage 1 はデフォルトで parent ごとに最大 1 候補のみ残します (`stage1_max_nhit_per_subject=1`)。その per-volume (`stage1_max_nhit_per_volume`) と in-total (`stage1_max_nhit_in_total`) の上限はデフォルトで無制限です。Stage 1・Stage 2 の制限は高コストな Stage 3 アラインメントの前に候補を絞り込むため、Stage 3 コストを抑えるレバーになります。一方 Stage 3 の制限はアラインメント後 (alnscore 基準) に適用されるため、アラインメント処理量自体は減らさず出力を絞り込みます。

1. **Stage 1 (候補選択):** クエリの各 k-mer に対して ID ポスティングをスキャンし、配列ごとに **coverscore** (配列にマッチしたクエリ k-mer 位置の数) を集計します。`stage1_min_score` 以上のスコアを持つ配列を候補として選出します。続いて候補集合を N → M → L の順 (いずれも coverscore 基準・タイ包含) で絞り込みます: `stage1_max_nhit_per_subject` (N) は各 (query, volume, strand) 内で parent ごとに上位 N 件を残し、`stage1_max_nhit_per_volume` (M) は (query, volume, strand) ごとに上位 M 件を残し、`stage1_max_nhit_in_total` (L) は全 volume・strand を通したクエリごとに上位 L 件を残します。

2. **Stage 2 (コリニアチェイニング):** 各候補に対して `.kpx` から位置レベルのヒットを収集し、対角線フィルタを適用した後、チェイニング DP により最良のコリニアチェインを求めます。チェインの長さが **chainscore** として報告されます。`chainscore >= stage2_min_score` のチェインが結果に含まれます。DP の内側ループは `-stage2_max_lookback` (デフォルト: 64) で制限され、各ヒットは直前の B 個のヒットのみを前駆候補として参照します。これにより、単一クエリ×サブジェクト間のヒット数が非常に多い場合の最悪計算量を O(n²) から O(n×B) に削減します。0 を指定すると探索窓が無制限 (O(n²)) になります。`-stage2_max_nhit_per_subject` が 1 より大きい値 (または 0 で無制限) の場合、貪欲な最良チェイン除去により同一サブジェクトから重複のない複数のチェインを抽出します: 最良チェインを見つけてそのヒットを除去し、残りのヒットで DP を再実行する処理を、制限に達するか `min_score` を満たすチェインがなくなるまで繰り返します。`-stage2_max_nhit_per_subject` (N) は 2 段階で適用されます — 抽出時に (query, fragment, strand) 単位で 1 回、overlap 同一領域の重複除外後に (query, parent[, strand]) 単位でもう 1 回 — いずれも chainscore を基準とします。`-stage2_max_nhit_per_subject_mode` (デフォルト: 3) は両段階を制御します: モード 1/2 はタイを扱わず上位 N 件を保持、モード 3/4 は上位 N 件に加え N 位の chainscore に並ぶチェインをすべて保持します; モード 1/3 は 2 本の strand を parent ごとに 1 つのグループへ統合し、モード 2/4 は strand を分離します。センチネル値 0 はモード 3 に解決され、明示値は 1〜4 のみ有効です。本オプションは `-mode 2` 以上が必要です。per-subject 選択の後、`-stage2_max_nhit_in_total` (L) が全 volume・strand を通したクエリごとに chainscore 基準で上位 L 件 (タイ包含、0=無制限) に絞り込み、Stage 3 へ進むチェイン数を制限します。

3. **Stage 3 (ペアワイズアライメント):** Stage 2 の各ヒットに対して、BLAST DB からサブジェクト部分配列を取得し (`-context_extend` による拡張オプション付き)、Parasail ライブラリを使って半大域ペアワイズアライメントを実行します (`-stage3_score_matrix` で指定されたスコア行列を使用、デフォルト: DEGMATCH)。全ヒットに対してアライメントスコア (**alnscore**) が計算されます。`-stage3_traceback 1` を指定すると、CIGAR 文字列、正スコア率、正スコア塩基数、負スコア数、ギャップ付きアライメント配列も計算されます。`-stage3_min_ppositive` と `-stage3_min_npositive` によるフィルタリングが可能です (トレースバックモードのみ)。アライメントと重複除外の後、`-stage3_max_nhit_per_subject` (N、デフォルト 1) が (qseqid, sseqid[, sstrand]) グループごとに alnscore 基準で上位 N 件を残し (`-stage3_max_nhit_per_subject_mode` がタイ/strand の扱いを Stage 1/2 と同様に制御)、続いて `-stage3_max_nhit_in_total` (L) が全 volume を通したクエリごとに alnscore 基準で上位 L 件 (タイ包含、0=無制限) に絞り込みます。これらの制限はアライメント後に適用されるため、アライメント処理量は減らさず結果集合を絞り込みます。サブジェクト部分配列は `-nthread` 全数によるヒット並列でフェッチされ、各 TBB タスク内では (volume, OID) 順に走査して mmap の連続アクセス局所性を保持します。

**コンテクスト拡張 (`-context_extend`):** 整数は塩基数、小数はそのヒットのクエリ長に対する倍率です (`2` は 2 塩基、`2.0` はクエリ長の 2 倍)。拡張量は Stage 2 チェインの**前後それぞれ**に適用され、親の範囲でクランプされます。したがって 100 bp のクエリに `-context_extend 2.0` を指定すると、100 bp のチェインは最大で `200 + 100 + 200` = 500 塩基まで広がります。クランプは前後独立で補償は行いません。親の先頭から始まるチェインでも、後方に伸ばす量は指定どおりのままです。

アライメントは半大域なので、無料の端ギャップでスキップできるのは各端につき 2 配列のうち片方だけです。そのためクエリと subject 窓の短い方が必ず全長消費されます。倍率指定では通常 窓 > クエリ長 になるため、両端がどこにも一致しないクエリでもクエリ全長がアラインされます (`qstart` = 1、`qend` = `qlen`)。subject に存在しない配列 (プライマー、アダプタ等) をクエリが含む可能性があり、報告される境界を一致領域で止めたい場合は `-context_extend 0` か小さい塩基数を指定してください。

**適応的 `-stage2_min_score` (デフォルト):** `-stage2_min_score 0` (デフォルト) の場合、最小チェインスコアはクエリごとに適応的に設定され、解決済みの Stage 1 閾値が使用されます。割合指定の `-stage1_min_score` (例: `0.5`) との組み合わせでは、各クエリの k-mer 構成に基づくクエリごとの適応的閾値が設定されます。絶対値指定の `-stage1_min_score` の場合は、その設定値がそのまま使用されます。固定閾値を使用する場合は `-stage2_min_score` に正の整数を指定してください。

**Mode 1 (Stage 1 のみ):** `-mode 1` を指定すると Stage 2, 3 が省略されます。`.kpx` ファイルへのアクセスが不要となり、I/O とメモリを節約できます。結果には Stage 1 スコアのみが含まれ、位置フィールド (qstart, qend, sstart, send) と chainscore は省略されます。ソート基準は Stage 1 スコアに強制されます。

**Mode 3 (全パイプライン):** `-mode 3` を指定すると全 3 段階が実行されます。BLAST DB が必要です (`-db` で指定、デフォルトはインデックスプレフィックスと同じ)。ソート基準は alnscore に自動設定されます。SAM/BAM 出力には `-mode 3` と `-stage3_traceback 1` の両方が必要です。

デフォルトではクエリのフォワード鎖とリバースコンプリメント鎖の両方を検索します。`-strand 1` でプラス (フォワード) 鎖のみ、`-strand -1` でマイナス (リバースコンプリメント) 鎖のみの検索に制限できます。

### スペースドシードテンプレートマスク

スペースドシードが有効な場合 (`-t > 0`)、discontiguous megablast 方式のビットマスクテンプレートを使って k-mer が抽出されます。各マスクは t 塩基のウィンドウから k 個の位置を選択します。各 (k, t) の組み合わせに対して **coding** (コーディング領域最適化) と **optimal** (非コーディング領域最適化) の 2 種類のテンプレートが利用可能です。インデックス構築時には `coding` または `optimal` を指定して単一テンプレートのインデックスを構築します。検索時には **both** オプションで coding と optimal の個別インデックスをマージして両方の結果を統合できます。

全テンプレートは discontiguous MegaBLAST テンプレートの設計原則に基づいています。**coding** テンプレートは第 2 コドン位置を最大限カバーする周期的「110」構造を持ち、余剰ギャップは第 1 コドン位置に配置されます。**optimal** テンプレートは対応する coding テンプレートとの重複を最小化するよう設計され、非周期的な構造を使用します。ブックエンド構造はテンプレート長に応じて変化します (t=13: 両端「111」、t=15: 前端「111」+後端「11」、t=18: 前端「111」or「11」+後端「11」)。k=11/12 のテンプレートは discontiguous MegaBLAST のネイティブテンプレートであり、k=8/9 のテンプレートは PCR アンプリコン解析に適した短いシード重みのために同じ原則に従って新規設計されたものです。

有効な (k, t) の組み合わせとマスク値:

| k | t | 種別 | マスク (バイナリ、左→右 = 位置 0〜t−1) | マスク (hex) |
|---|---|---|---|---|
| 8 | 13 | coding | `1101100101101` | 0x1B2D |
| 8 | 13 | optimal | `1110011010011` | 0x1CD3 |
| 8 | 15 | coding | `100101100101101` | 0x4B2D |
| 8 | 15 | optimal | `111010010010011` | 0x7493 |
| 8 | 18 | coding | `100100100100101101` | 0x2492D |
| 8 | 18 | optimal | `110010001010010011` | 0x32293 |
| 9 | 13 | coding | `1101101101101` | 0x1B6D |
| 9 | 13 | optimal | `1110111010011` | 0x1DD3 |
| 9 | 15 | coding | `100101101101101` | 0x4B6D |
| 9 | 15 | optimal | `111010011010011` | 0x74D3 |
| 9 | 18 | coding | `100100101100101101` | 0x24B2D |
| 9 | 18 | optimal | `111010001010010011` | 0x3A293 |
| 11 | 16 | coding | `1101101101101101` | 0xDB6D |
| 11 | 16 | optimal | `1110010110110111` | 0xE5B7 |
| 11 | 18 | coding | `101101100101101101` | 0x2D96D |
| 11 | 18 | optimal | `111010010110010111` | 0x3A597 |
| 11 | 21 | coding | `100101100101100101101` | 0x12CB2D |
| 11 | 21 | optimal | `111010010100010010111` | 0x1D2897 |
| 12 | 16 | coding | `1111101101101101` | 0xFB6D |
| 12 | 16 | optimal | `1110110110110111` | 0xEDB7 |
| 12 | 18 | coding | `101101101101101101` | 0x2DB6D |
| 12 | 18 | optimal | `111010110010110111` | 0x3ACB7 |
| 12 | 21 | coding | `100101101101100101101` | 0x12DB2D |
| 12 | 21 | optimal | `111010010110010010111` | 0x1D2C97 |

**KmerInt 型の選択:** k-mer 値に使用する整数型は必要ビット幅 `2*k` で決定されます。ビット幅が 16 を超える場合は `uint32_t`、以下の場合は `uint16_t` が使用されます。したがって k <= 8 では `uint16_t` (16 ビット以下)、k >= 9 では `uint32_t` が使用されます。

### 高頻度 k-mer フィルタリング

高頻度 k-mer フィルタリングはインデックス構築時の `-max_freq_build` のみで行われます。検索時の高頻度フィルタは存在せず、検索ホットパスでクエリごとにボリューム横断の k-mer カウント集計を行うコストは発生しません。

**構築時除外** (`-max_freq_build`): `-max_freq_build` を指定してインデックスを構築すると、高頻度 k-mer がインデックスから完全に除外されます。k-mer カウントは全ボリュームで合算された後に閾値と比較されるため、各ボリュームでは閾値未満だが合計では閾値を超える k-mer も正しく除外されます。除外された k-mer は共有 `.khx` ファイル (k 値ごとに 1 つ、ボリュームごとではない) に記録されます。小数値 (0 < x < 1) を指定した場合、閾値は全ボリューム合計 NSEQ に基づいて解決されます。検索時に割合指定の `-stage1_min_score` を使用する場合、構築時に除外された k-mer が `.khx` ファイルから認識され、閾値計算 (下記の `Nhighfreq`) から差し引かれます。

### 割合指定の Stage 1 閾値

`-stage1_min_score` を小数 (0 < P < 1) で指定すると、閾値はクエリごとに以下の式で解決されます:

```
threshold = ceil(Nqkmer * P) - Nhighfreq
```

各変数の意味:
- **Nqkmer**: 純粋な窓カウント `max(0, seq_len - span + 1)`。spaced seed (`-t > 0`) のときは `span = t`、それ以外は `span = k`。クエリ内容には依存しません。
- **Nhighfreq**: 検索に寄与しない窓位置の総数で、以下 3 ケースの和集合として定義されます:
  1. **構築時除外**: 当該位置で展開された k-mer のうち**いずれか 1 つでも** `.khx` で除外されている。
  2. **縮重展開超過**: 当該位置の窓に含まれる IUPAC 縮重塩基の展開積が `-max_degen_expand` を超えるため emit されない。
  3. **窓内 truly-invalid 文字**: 当該窓に `[ACGT]` ∪ IUPAC 以外の文字が含まれる (注: クエリ全体に invalid 文字があると、閾値計算前に `kSkipInvalidChar` で全体がスキップされます。下記参照)。

  ikafssn の実装は `Nhighfreq = #{ANY 除外された位置} + (Nqkmer − #emit された位置数)` と等価な計算を行います。

検索対象のすべてのストランドで閾値が `< 1` に落ちる場合、サイレントに 0 ヒットを返すのではなく `kSkipThresholdUnreachable` でクエリ単位スキップになります (下記参照)。

### スキップ理由

クエリが Stage 1〜3 を通せない場合、ikafssn は理由を表すマーカーレコードを各出力フォーマットに emit します:

| 理由 (enum / 文字列) | 意味 |
|---|---|
| `kSkipQueryTooShort` / `query_too_short` | `seq_len < span` (k-mer / spaced seed の窓が乗らない)、または `seq_len < -min_query_length`。 |
| `kSkipDegenRejected` / `degen_rejected` | `-accept_qdegen 0` 指定時に IUPAC 縮重塩基を含む。 |
| `kSkipInvalidChar` / `invalid_char` | `[ACGT]` ∪ IUPAC 以外の文字を含む。詳細メッセージに位置を含む。 |
| `kSkipThresholdUnreachable` / `threshold_unreachable` | フラクショナル閾値が検索対象のすべてのストランドで `< 1` になる。 |
| `kSkipQueryTooLong` / `query_too_long` | フラグメント分割インデックス で `seq_len > overlap_length`。親相対 dedup キーは「すべてのチェインヒットが隣接 2 フラグメント内に収まること」を前提とするため、それより長いクエリはフラグメント単位の partial chain を生成する前に弾かれます。`min_length_split == 0` (分割無効) のインデックスは `overlap_length == 0` を報告するためこのチェックが無効化されます。 |

各出力フォーマットでの表現:

- **TSV**: 1 行のスキップマーカーが emit され、`sseqid = "*SKIPPED:<reason>"`、`sstrand = '*'`、数値列はすべて 0、`qlen` は元の長さ。`result_reader` はこのマーカー行を silently に読み飛ばすので既存パイプラインへの影響はありません。
- **JSON**: 該当クエリオブジェクトに `"status": "skipped"`, `"skip_reason": "<reason>"`, `"skip_detail": "<detail>"` が付き、`"hits": []`。通常クエリは `"status": "ok"`。
- **SAM/BAM**: unmapped レコード (`FLAG = 4`, `RNAME = *`, `POS = 0`, `CIGAR = *`) として出力され、aux タグ `XR:Z:<reason>` と `XD:Z:<detail>` が付く。`samtools view -f 4` で抽出可能。
- **stderr**: preprocess 時に `Warning: query 'X' skipped: <reason> (<detail>)` が出力されます。

スキップが 1 件以上発生した検索は終了コード `2` を返します。

### `-template_type both`: クロステンプレート統合

spaced seed で `-template_type both` を指定すると、ikafssn は coding と optimal の両インデックスを同時に検索し、Stage 1 スコアを `[0, Nqkmer]` に収まる単一の統合スコアにマージします。これにより、テンプレートごとに独立してスコアリングする素朴な方式の 2 つの問題を回避します:

1. **クロステンプレート二重カウント**: ある配列位置が coding と optimal の両 mask にマッチした場合、素朴な方式では当該クエリ位置に `+2` が加算され、スコアレンジが `[0, 2 × Nqkmer]` に膨らみます。ikafssn は共有の (sid, q_pos) 単位 dedup バッファを通じて両ストリームを `q_pos` 昇順でインターリーブ処理するため、各クエリ位置は最大 `+1` しか寄与しません。
2. **ハイブリッドクエリの取りこぼし**: クエリの領域 A は coding のみ、領域 B は optimal のみマッチするケースは、加算的な `thr_cod + thr_opt` 閾値では検出を逃します。ikafssn は代わりに `min(thr_cod, thr_opt)` の統合閾値を適用するため、ハイブリッドマッチも保持されます。

これにより、`-template_type both` で報告される Stage 1 スコアの意味は `coding` 単独 / `optimal` 単独と同じです: 構築時 `.khx` 除外を考慮したうえでサブジェクトとマッチしたクエリ窓位置の数で、上限は `Nqkmer`。`Nhighfreq` は各サイドで独立に計算されます。

### スコア種別

ikafssn は 3 種類のスコアを計算します。

| スコア | 説明 | 計算ステージ |
|---|---|---|
| **coverscore** | 参照配列にマッチしたクエリ k-mer 位置の数。各クエリ位置は参照配列あたり最大 1 回カウントされます (その k-mer が参照配列側の複数位置に出現しても重複計上されません)。これが Stage 1 スコアです。 | Stage 1 |
| **chainscore** | チェイニング DP が求めた最良コリニアチェインの長さ (k-mer ヒット数)。`.kpx` の位置データを使用します。 | Stage 2 |
| **alnscore** | Parasail による半大域ペアワイズアライメントスコア (`-stage3_score_matrix` で指定、デフォルト: degmatch)。BLAST DB からのサブジェクト配列取得が必要です。 | Stage 3 |

- Stage 1 スコアは常に coverscore です。
- ソート基準はモードにより自動決定: mode 1 は coverscore、mode 2 は chainscore、mode 3 は alnscore。
- `-mode 1` では Stage 1 スコアのみが利用可能で、chainscore と alnscore は計算されません。

### Stage 3 スコア行列

Stage 3 のペアワイズアライメントで使用するスコア行列は `-stage3_score_matrix` オプションで選択できます。利用可能な行列は以下の 3 種類です。

| 行列名 | 説明 |
|---|---|
| **degmatch** (デフォルト) | 縮重塩基ペアが少なくとも 1 つのヌクレオチドを共有する場合に正のスコアを付与します。プライマーマッチングや縮重塩基を含む検索に適しています。 |
| **dnafull** | EMBOSS の DNAfull 行列 (Todd Lowe 作成)。あいまいなヌクレオチドコードを確率に基づいて最も近い整数に丸めた値を使用します。部分重複する縮重ペアには負または低いスコアが付与されます。 |
| **nuc44** | NCBI BLAST のヌクレオチド行列。dnafull と類似していますが、縮重塩基のスコアリングがやや異なります。 |

**CIGAR `=`/`X` 判定:** 拡張 CIGAR オペレータ `=` (sequence match) と `X` (sequence mismatch) は、厳密な塩基一致ではなくアライメントスコアに基づいて判定されます。スコア行列で正のスコア (score > 0) が付与される位置は `=`、0 以下のスコアが付与される位置は `X` として報告されます。つまり、DEGMATCH 行列では N-A のような縮重塩基ペア (スコア 1) は match (`=`) としてカウントされますが、NUC44/DNAFULL では mismatch (`X`) としてカウントされます。

**カラム名:** 出力カラムは `ppositive` (正スコア率)、`npositive` (正スコア塩基数)、`nnegative` (負スコア数) です。DEGMATCH 行列では、これらのカウントは厳密な一致位置ではなく正スコア位置を表します。対応するフィルタオプションは `-stage3_min_ppositive` / `-stage3_min_npositive` です。

**アラインメント領域のみを報告します:** Stage 3 はコンテキスト拡張したサブジェクト窓に対してクエリをアラインするため、生のアラインメントには窓の残りを埋める末端ギャップが含まれます。報告前にこれを除去するので、`qstart` / `qend` / `sstart` / `send`、`cigar`、`qseq` / `sseq`、および `ppositive` の分母はいずれも末端ギャップを除いたアラインメント領域を表し、同じアラインメントに対して BLAST が報告する値と一致します。したがって `ppositive` は窓の長さではなくアラインメント長に対する `npositive` の割合です。

**DEGMATCH 行列 (16×16):**

行と列は 16 種類のシンボルに対応します: A, T, G, C, S (G/C), W (A/T), R (A/G), Y (T/C), K (G/T), M (A/C), B (G/T/C), V (A/G/C), H (A/T/C), D (A/G/T), N (A/T/G/C), `*` (stop/invalid)。

|   | A | T | G | C | S | W | R | Y | K | M | B | V | H | D | N | * |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **A** | 5 | -4 | -4 | -4 | -4 | 3 | 3 | -4 | -4 | 3 | -4 | 2 | 2 | 2 | 1 | -4 |
| **T** | -4 | 5 | -4 | -4 | -4 | 3 | -4 | 3 | 3 | -4 | 2 | -4 | 2 | 2 | 1 | -4 |
| **G** | -4 | -4 | 5 | -4 | 3 | -4 | 3 | -4 | 3 | -4 | 2 | 2 | -4 | 2 | 1 | -4 |
| **C** | -4 | -4 | -4 | 5 | 3 | -4 | -4 | 3 | -4 | 3 | 2 | 2 | 2 | -4 | 1 | -4 |
| **S** | -4 | -4 | 3 | 3 | 4 | -4 | 1 | 1 | 1 | 1 | 2 | 2 | 1 | 1 | 1 | -4 |
| **W** | 3 | 3 | -4 | -4 | -4 | 4 | 1 | 1 | 1 | 1 | 1 | 1 | 2 | 2 | 1 | -4 |
| **R** | 3 | -4 | 3 | -4 | 1 | 1 | 4 | -4 | 1 | 1 | 1 | 2 | 1 | 2 | 1 | -4 |
| **Y** | -4 | 3 | -4 | 3 | 1 | 1 | -4 | 4 | 1 | 1 | 2 | 1 | 2 | 1 | 1 | -4 |
| **K** | -4 | 3 | 3 | -4 | 1 | 1 | 1 | 1 | 4 | -4 | 2 | 1 | 1 | 2 | 1 | -4 |
| **M** | 3 | -4 | -4 | 3 | 1 | 1 | 1 | 1 | -4 | 4 | 1 | 2 | 2 | 1 | 1 | -4 |
| **B** | -4 | 2 | 2 | 2 | 2 | 1 | 1 | 2 | 2 | 1 | 3 | 1 | 1 | 1 | 1 | -4 |
| **V** | 2 | -4 | 2 | 2 | 2 | 1 | 2 | 1 | 1 | 2 | 1 | 3 | 1 | 1 | 1 | -4 |
| **H** | 2 | 2 | -4 | 2 | 1 | 2 | 1 | 2 | 1 | 2 | 1 | 1 | 3 | 1 | 1 | -4 |
| **D** | 2 | 2 | 2 | -4 | 1 | 2 | 2 | 1 | 2 | 1 | 1 | 1 | 1 | 3 | 1 | -4 |
| **N** | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | -4 |
| **\*** | -4 | -4 | -4 | -4 | -4 | -4 | -4 | -4 | -4 | -4 | -4 | -4 | -4 | -4 | -4 | -4 |

スコアルール:
- 確定塩基同士の一致 (A-A, T-T, G-G, C-C) = 5
- 確定塩基同士の不一致 (例: A-T, A-G) = -4
- 2 重縮重の自己対 (S-S, W-W, R-R, Y-Y, K-K, M-M) = 4
- 2 重縮重 vs 含有される確定塩基 (例: R-A, R-G) = 3
- 2 重縮重 vs 部分重複する別の 2 重縮重 (例: R-K) = 1
- 2 重縮重 vs 重複なし (例: S-W) = -4
- 3 重縮重の自己対 (B-B, V-V, H-H, D-D) = 3
- 3 重縮重 vs 含有される確定塩基/2 重縮重 (例: B-T, B-K) = 2
- 3 重縮重 vs 部分重複 = 1
- N vs 任意 (`*` 以外) = 1, N-N = 1
- `*` vs 任意 = -4

## 出力形式

**座標規約:** TSV / JSON 出力は BLAST `-outfmt 6` に準拠します。

- `qstart` / `qend` / `sstart` / `send` の 4 つはすべて **1-based 包括**です。
- `qstart` / `qend` は **クエリ相対座標** (与えたクエリ配列上の位置) なので、どちらの鎖でも `qstart` < `qend` です。
- `sstart` / `send` は正鎖番号での**親相対座標**ですが、アライメントの進行方向に並ぶため、逆鎖ヒット (`sstrand` が `-`) では **`sstart` > `send`** になります。

`sseqid` は **親 OID** のアクセッション (フラグメント由来の合成名は使われない)、`slen` は親 OID の全長です。Stage 2 / Stage 3 の dedup によりフラグメント単位のチェインが畳み戻されるため、出力は `(qseqid, sseqid, sstrand, send, alnscore)` の組ごとに 1 行の正準ヒットのみを含みます。

SAM / BAM 出力はこの規約に従いません。`POS` は htslib が扱う 0-based の最左座標であり、鎖によって入れ替わることはありません ([SAM/BAM 形式](#sambam-形式) を参照)。

### Tab 形式 (デフォルト)

**Mode 2** (デフォルト):

タブ区切りのカラムです (Stage 1 スコアのカラムは常に `coverscore`)。

```
# qseqid  sseqid  sstrand  qstart  qend  qlen  sstart  send  slen  coverscore  chainscore  volume
```

**Mode 1** (`-mode 1`):

```
# qseqid  sseqid  sstrand  qlen  slen  coverscore  volume
```

**Mode 3, traceback=0** (`-mode 3`):

```
# qseqid  sseqid  sstrand  qend  qlen  send  slen  coverscore  chainscore  alnscore  volume
```

注: トレースバックなしでは正確なアライメント開始位置が得られないため、`qstart` と `sstart` は省略されます。逆鎖ではアライメントの終端が subject 座標の小さい側になるため、この形式の `send` は traceback 形式が報告する 2 つの subject 座標のうち小さい方に対応します。

**Mode 3, traceback=1** (`-mode 3 -stage3_traceback 1`):

```
# qseqid  sseqid  sstrand  qstart  qend  qlen  sstart  send  slen  coverscore  chainscore  alnscore  ppositive  npositive  nnegative  cigar  qseq  sseq  volume
```

### JSON 形式

`results` には `qseqid` ごとに 1 つのオブジェクトが、入力中で最初に現れた順に並びます。各オブジェクトの `hits` は検索が生成した順序のままです。同じ `qseqid` を持つクエリ配列はそれぞれ独立に検索されヒットも別々に得られますが、オブジェクトは名前をキーとするため 1 つに統合されます。

**Mode 2** (デフォルト):

```json
{
  "results": [
    {
      "qseqid": "query1",
      "hits": [
        {
          "sseqid": "NC_001234.5",
          "sstrand": "+",
          "qstart": 1,
          "qend": 150,
          "qlen": 200,
          "sstart": 1001,
          "send": 1150,
          "slen": 5000,
          "coverscore": 8,
          "chainscore": 12,
          "volume": 0
        }
      ]
    }
  ]
}
```

**Mode 1** (`-mode 1`):

```json
{
  "results": [
    {
      "qseqid": "query1",
      "hits": [
        {
          "sseqid": "NC_001234.5",
          "sstrand": "+",
          "qlen": 200,
          "slen": 5000,
          "coverscore": 8,
          "volume": 0
        }
      ]
    }
  ]
}
```

**Mode 3, traceback=0** (`-mode 3`):

```json
{
  "results": [
    {
      "qseqid": "query1",
      "hits": [
        {
          "sseqid": "NC_001234.5",
          "sstrand": "+",
          "qend": 150,
          "qlen": 200,
          "send": 1150,
          "slen": 5000,
          "coverscore": 8,
          "chainscore": 12,
          "alnscore": 240,
          "volume": 0
        }
      ]
    }
  ]
}
```

**Mode 3, traceback=1** (`-mode 3 -stage3_traceback 1`):

```json
{
  "results": [
    {
      "qseqid": "query1",
      "hits": [
        {
          "sseqid": "NC_001234.5",
          "sstrand": "+",
          "qstart": 1,
          "qend": 150,
          "qlen": 200,
          "sstart": 1001,
          "send": 1150,
          "slen": 5000,
          "coverscore": 8,
          "chainscore": 12,
          "alnscore": 240,
          "ppositive": 95.3,
          "npositive": 143,
          "nnegative": 7,
          "cigar": "50=2X48=1I50=",
          "qseq": "ACGT...",
          "sseq": "ACGT...",
          "volume": 0
        }
      ]
    }
  ]
}
```

### SAM/BAM 形式

SAM/BAM 出力には `-mode 3 -stage3_traceback 1` が必要です。`-output_format sam` で SAM、`-output_format bam` で BAM を出力します (BAM は `-o <path>` が必須)。

SAM レコードの構成:
- **QNAME**: qseqid
- **FLAG**: 0 (フォワード) / 16 (リバース)
- **RNAME**: sseqid
- **POS**: アライメント最左の親相対位置 (1-based)。どちらの鎖でもアライメントの subject 座標の**小さい方**であり、逆鎖で `send` と入れ替わる TSV / JSON の `sstart` とは異なります。
- **MAPQ**: 255
- **CIGAR**: 拡張 CIGAR (=/X/I/D 演算子)。常にリファレンスの向きです。アラインメント領域外のクエリ塩基は先頭・末尾のハードクリップ (`H`) として出力されます (SEQ がアラインメント領域しか持たないため)。
- **SEQ**: ギャップなしのアラインメント領域のクエリ配列。リファレンスの向きで出力されるため、FLAG 0x10 が立つ場合は SAM 仕様どおり与えたクエリの逆相補になります。
- **QUAL**: * (利用不可)
- **タグ**: `AS:i` (alnscore), `NM:i` (nnegative), `cs:i` (chainscore), `cv:i` (coverscore)

## デプロイ構成

### 単一マシン

```
ikafssnsearch (スタンドアロン、サーバ不要)
    または
ikafssnserver → ikafssnclient (UNIX ソケット経由)
```

### 複数マシン

```
マシン A (検索サーバ):
  ikafssnserver -ix ./index/mydb -tcp 0.0.0.0:9100

マシン B (HTTP フロントエンド):
  ikafssnhttpd -server_tcp A:9100 -listen :8080
  nginx (TLS 終端、認証、レート制限) → ikafssnhttpd
```

### 複数データベース

1 つの `ikafssnserver` プロセスで複数のデータベースを同時にサーブできます:

```
# 1 プロセスでマルチ DB (推奨)
ikafssnserver -ix ./nt_index -db nt -ix ./rs_index -db refseq_genomic \
    -socket /var/run/ikafssn.sock

ikafssnclient -socket /var/run/ikafssn.sock -ix nt -query query.fasta
ikafssnclient -socket /var/run/ikafssn.sock -ix refseq_genomic -query query.fasta
```

### マルチバックエンド (負荷分散 / マルチサーバ)

1 つの `ikafssnhttpd` で複数の `ikafssnserver` インスタンスに接続できます:

```
# 負荷分散: 同じ DB を 2 台のサーバで提供、1 つの httpd
ikafssnserver -ix ./nt_index -db nt -tcp 0.0.0.0:9100   # サーバ A
ikafssnserver -ix ./nt_index -db nt -tcp 0.0.0.0:9100   # サーバ B
ikafssnhttpd -server_tcp A:9100 -server_tcp B:9100 -listen :8080

# 異なる DB を別サーバで提供、1 つの httpd で統合
ikafssnserver -ix ./nt_index -db nt -socket /var/run/nt.sock
ikafssnserver -ix ./rs_index -db refseq -socket /var/run/rs.sock
ikafssnhttpd -server_socket /var/run/nt.sock -server_socket /var/run/rs.sock -listen :8080
```

同じデータベース名が複数バックエンドに存在する場合、`ikafssnhttpd` は起動時に共通する (db, k) の組み合わせについて合計配列数・合計塩基数の一致を検証します。バックエンド間で k 値セットが異なることは許容されます。統合された能力情報には全バックエンドの k 値グループの和集合が含まれるため、一部のバックエンドのみが持つ k 値へのリクエストは自然にそのバックエンドにルーティングされます。リクエストは有効な空き容量 (スロット空き数と `max_nseq_per_req` の小さい方) がある最高優先度のバックエンドにルーティングされます。容量値 (`max_queue_size`、`queue_depth`、`max_nseq_per_req`) はサーバ単位で共有されるため、同一サーバ上の複数データベース間で共有される点に注意してください。

代替構成: 別プロセスで HTTP パスベースルーティング:

```
ikafssnserver -ix ./nt_index  -socket /var/run/ikafssn_nt.sock
ikafssnserver -ix ./rs_index  -socket /var/run/ikafssn_rs.sock

ikafssnhttpd -server_socket /var/run/ikafssn_nt.sock -listen :8080 -path_prefix /nt
ikafssnhttpd -server_socket /var/run/ikafssn_rs.sock -listen :8081 -path_prefix /rs
# nginx 側で /nt → :8080, /rs → :8081 にルーティング
```

### systemd との統合

サンプルの systemd ユニットファイルが `doc/systemd/` に提供されています。詳細は各ファイルを参照してください。

## インデックスファイル形式

インデックスの**内容**を変えるフラグメント・インデックス用パラメータは、すべてファイル名に符号化されます。これにより、異なるパラメータで構築した複数のインデックスを同一ディレクトリに共存できます。

BLAST DB ボリュームごとに 3 ファイルが生成されます:

```
<vol_basename>.k<k>[.t<t>][.minlen<X>][.minsplit<X>][.ovllen<X>]
              [.maxfreq<X>][.maxexpand<X>][.cod|.opt].(kix|kpx|ksx)
```

DB 単位の共有ファイル（マニフェスト `.kvx`、および `-max_freq_build` 使用時の除外ビットセット `.khx`）は同じ命名規則からボリューム番号を除いた形で生成されます:

```
<db_base>.k<k>[.t<t>][.minlen<X>][.minsplit<X>][.ovllen<X>]
          [.maxfreq<X>][.maxexpand<X>][.cod|.opt].(kvx|khx)
```

サフィックスの省略条件:

| サフィックス | 省略条件 |
|---|---|
| `k<X>`         | 省略しない（常に出力、ゼロ埋めなし） |
| `t<X>`         | `t == 0` |
| `minlen<X>`    | `min_seq_length == 0` |
| `minsplit<X>`  | `min_length_split == 0` |
| `ovllen<X>`    | `overlap_length == 0`（`minsplit` と同期） |
| `maxfreq<X>`   | `max_freq_build == 1`（解決された絶対閾値） |
| `maxexpand<X>` | `max_degen_expand` が `0` または `1` |
| `cod` / `opt`  | `t == 0` |

例:

- デフォルトの `-mode 1`（`-min_seq_length 64 -min_length_split 50000 -overlap_length 500`）:
  - `nt.00.k11.minlen64.minsplit50000.ovllen500.kix`
  - `nt.k11.minlen64.minsplit50000.ovllen500.kvx`
- `-mode 1` + クロスボリュームフィルタ（`-max_freq_build 1000 -max_degen_expand 4`）:
  - `nt.00.k11.minlen64.minsplit50000.ovllen500.maxfreq1000.maxexpand4.kix`
  - `nt.k11.minlen64.minsplit50000.ovllen500.maxfreq1000.maxexpand4.khx`
- `-mode 2` / `-mode 3`（`-min_length_split 0 -overlap_length 0`）:
  - `nt.00.k11.minlen64.kix`
- スペースドシード（`-k 11 -t 16 -template_type coding -mode 1`）:
  - `nt.00.k11.t16.minlen64.minsplit50000.ovllen500.cod.kix`
  - `nt.k11.t16.minlen64.minsplit50000.ovllen500.cod.kvx`

`-max_freq_build` に小数（例: `0.001`）を指定した場合、メタデータ収集パスの終了直後に絶対閾値へ解決され、その**絶対値**がファイル名に出力されます。よって、同じ小数を異なるフラグメント集合に適用すると別のファイル名になります。

`.khx` ファイルは 32 バイトヘッダ（マジック "KMHX"、フォーマットバージョン、k）に続き `ceil(4^k / 8)` バイトのビットセットで構成され、ビット *i* = 1 は k-mer *i* がインデックス構築時に除外されたことを示します。

ID ポスティングと位置ポスティングは別ファイルに格納されるため、Stage 1 フィルタリングが `.kpx` にアクセスすることはなく、ページキャッシュ効率が最大化されます。

**マルチアクセッション defline:** 元の BLAST DB が `makeblastdb -parse_seqids` で構築され、`\x01` (`^A`) 区切りの multi-defline レコード（同一塩基配列を複数アクセッションで登録する NCBI 慣習）を含んでいる場合、`ikafssnindex` は OID ごとに**全てのアクセッション**を保持します。`.ksx` のアクセッション文字列は OID ごとに全アクセッションを `\x01` で連結した形で格納され、検索出力 (`sseqid` カラム / SAM RNAME / FASTA defline / プロトコルの `sseqid` フィールド) も同じ `\x01` 連結形のまま emit されます。受け取り側で `\x01` を区切りとして分割してください。`-seqidlist` フィルタと `ikafssnretrieve` は `\x01` 連結形・個別アクセッションのいずれを指定しても解決できます。

**インデックスフォーマットバージョン:** 現在のインデックスフォーマットは全ファイル (`.kix`、`.kpx`、`.ksx`、`.khx`、`.kvx`) で **v11** です。レイアウト概要:

- **`.kix` / `.kpx` マジック**は `KIX11` / `KPX11`。フラグメント・インデックス用パラメータ（`min_seq_length`、`min_length_split`、`overlap_length`、解決済みの `max_freq_build` / `max_degen_expand`）はヘッダから移動し、ファイル名で表現される設計に変更されました。検索側はボリューム検出時に一度だけ parse して `-min_query_length` / `-max_query_length` のチェックに利用します。
- **`.ksx` 二段レイアウト:** まず各親 OID の `(parent_length, blast_oid, accession)` を親テーブルに記録し、続いてフラグメントテーブルが各内部 SeqId を `(parent_idx, fragment_start, fragment_end)` にマップします。`KsxReader::seq_length` / `accession` などの簡易アクセサは SeqId を受け取り、内部で対応する親に解決します。マジックは `KMSX`。
- **`.kix` ボディ:** Elias-Fano 辞書。各 posting list は `[u32 distinct_count]` に続いて、distinct seq_id のデルタ列を FastPFor `CompositeCodec<SIMDFastPFor<4>, VariableByte>` で符号化したボディが並びます。
- **`.kpx` ボディ:** Elias-Fano 辞書。各 posting list は 2-bit kind map から直接始まり（posting list ごとのヘッダはありません）、partition / short_occ1 / short_occ_ge2 サブバケットが続きます。いずれも frame-of-reference ブロックでビットパックされています。デコードは候補集合ドリブンで、必要な位置を含むブロックだけを展開します。
- **`.kvx`:** マニフェストテキスト形式。`FORMAT_VERSION` 行は 11。

フラグメント分割器は kafssstore の `split_long_sequence_positions` の C++ ポートで、ncbi4na==0xF カットと calcsegment2 公式に従います。`-min_length_split 0` の場合は各親が単一フラグメントとして登録され、ファイル名にも `minsplit` / `ovllen` サフィックスは出力されません。

`format_version` が一致しないインデックスは open 時に以下のようなメッセージで拒否されます:

```
KixReader: index format version mismatch (got 10, expected 11). Please rebuild with the current ikafssnindex.
```

(`.kpx` / `.ksx` / `.khx` / `.kvx` も同様)。現行の `ikafssnindex` で再構築してください:

```
# デフォルト (mode 1 + フラグメント分割有効)
ikafssnindex -db <BLAST_DB_prefix> -k <k> -o <out_dir>

# 分割を無効化 (1 親あたり 1 フラグメント)
ikafssnindex -db <BLAST_DB_prefix> -k <k> -o <out_dir> \
    -min_length_split 0 -overlap_length 0
```

NCBI nt 規模（~700 ボリューム、k=12 t=21）では 32 コアホストで数十時間程度の所要時間を想定してください。`-min_length_split` を非 0 にするとオーバーラップに比例してディスク上のサイズは増えますが、染色体級の親 OID は複数の有限長フラグメントに分割されるため Stage 2 のサブジェクト単位メモリ上限が抑えられます。

## インストール

### macOS 26 Tahoe (Homebrew)

macOS 26 (Tahoe) の Apple Silicon (aarch64) 環境では、Homebrew Tap 経由でインストールできます:

```bash
brew tap astanabe/ikafssn
brew install ikafssn
```

ビルド済み Bottle が利用可能な場合はそれがインストールされ、利用できない場合はソースからビルドされます。

### Conda (チャンネル)

ikafssn は専用の Conda チャンネル `https://conda.ikafssn.org` 経由で配布されます。`linux-64`、`linux-aarch64`、`osx-arm64` 向けのビルド済み `.conda` パッケージが提供されており、Linux 版は glibc 2.34 以上 (RHEL 9 / Ubuntu 22.04 / Debian 12 以降) を、`osx-arm64` 版は macOS 11 以上を対象としています。

既存の miniforge / miniconda 環境に `mamba` (または `conda`) でインストールします。実行時の依存パッケージは `conda-forge` から取得されます:

```bash
mamba install -c https://conda.ikafssn.org -c conda-forge ikafssn
```

専用環境を作成する場合:

```bash
mamba create -n ikafssn -c https://conda.ikafssn.org -c conda-forge ikafssn
```

このチャンネル自体は `repodata.json` インデックスのみを配信しています ([astanabe/conda-ikafssn](https://github.com/astanabe/conda-ikafssn) リポジトリを GitHub Pages で公開)。実際の `.conda` バイナリは [CEP-15 `base_url`](https://github.com/conda/ceps/blob/main/cep-0015.md) 機構により本プロジェクトの [GitHub Releases](https://github.com/astanabe/ikafssn/releases) から直接ダウンロードされます。

要件:
- `conda` 23.7 以上 または `micromamba` 2.0 以上 (これより古いクライアントは `base_url` を無視するため、パッケージダウンロード時に HTTP 404 で失敗します)。

### Ubuntu (APT チャンネル)

ikafssn は署名付き APT チャンネル `https://deb.ikafssn.org` 経由で配布されます。チャンネルは常に最新リリースのみを保持し、過去バージョンは保持されません。

| Suite | Ubuntu リリース | アーキテクチャ |
|---|---|---|
| `jammy` | Ubuntu 22.04 LTS | `amd64`, `arm64` |
| `noble` | Ubuntu 24.04 LTS | `amd64`, `arm64` |
| `resolute` | Ubuntu 26.04 LTS | `amd64`, `arm64` |

初回セットアップ (root または `sudo` で実行):

```bash
sudo install -d -m 0755 /etc/apt/keyrings
curl -fsSL https://deb.ikafssn.org/ikafssn-archive-keyring.asc \
  | sudo gpg --dearmor -o /etc/apt/keyrings/ikafssn-archive-keyring.gpg
echo "deb [signed-by=/etc/apt/keyrings/ikafssn-archive-keyring.gpg] https://deb.ikafssn.org/ $(lsb_release -cs) main" \
  | sudo tee /etc/apt/sources.list.d/ikafssn.list
sudo apt update
sudo apt install ikafssn
```

以降のアップグレード:

```bash
sudo apt update && sudo apt upgrade ikafssn
```

チャンネルは `.deb` バイナリ本体と suite 別の `Packages` / `Release` / `InRelease` / `Release.gpg` メタデータを配信します。公開鍵のフィンガープリントは [astanabe/ikafssn](https://github.com/astanabe/ikafssn) リポジトリ root に commit されている `ikafssn-archive-keyring.asc` でも検証できます。

### Enterprise Linux (DNF / YUM チャンネル)

ikafssn は署名付き DNF チャンネル `https://rpm.ikafssn.org` 経由で配布されます。チャンネルは常に最新リリースのみを保持し、過去バージョンは保持されません。

| `releasever` | ディストリビューション | アーキテクチャ |
|---|---|---|
| `9` | AlmaLinux / Rocky Linux / RHEL 9 | `x86_64`, `aarch64` |
| `10` | AlmaLinux / Rocky Linux / RHEL 10 | `x86_64`, `aarch64` |

初回セットアップ (root または `sudo` で実行):

```bash
sudo tee /etc/yum.repos.d/ikafssn.repo > /dev/null <<'EOF'
[ikafssn]
name=ikafssn
baseurl=https://rpm.ikafssn.org/el$releasever/$basearch/
enabled=1
gpgcheck=1
repo_gpgcheck=1
gpgkey=https://rpm.ikafssn.org/ikafssn-archive-keyring.asc
EOF
sudo dnf install ikafssn
```

以降のアップグレード:

```bash
sudo dnf upgrade ikafssn
```

パッケージ署名 (`gpgcheck=1`) とリポジトリメタデータ署名 (`repo_gpgcheck=1`) の両方がチャンネル公開鍵で検証されます。

### 直接ダウンロード (フォールバック)

Homebrew / Conda / APT / DNF が使用できない環境向けに、同じ `.deb`、`.rpm`、`.conda`、Bottle 成果物を [GitHub Releases](https://github.com/astanabe/ikafssn/releases) ページから直接ダウンロードできます。

パッケージ命名規則:

```
ikafssn_<version>_ubuntu-<ubuntu_ver>_<arch>.deb         # Ubuntu 22.04 / 24.04 / 26.04, amd64 / arm64
ikafssn-<version>.el<el_ver>.<arch>.rpm                   # EL 9 / 10, x86_64 / aarch64
ikafssn_<version>_<conda_subdir>.conda                    # linux-64 / linux-aarch64 / osx-arm64
ikafssn-<version>.arm64_tahoe.bottle.tar.gz               # macOS 26 (Tahoe) arm64 用 Homebrew Bottle
```

ダウンロードした `.deb` のインストール:

```bash
sudo apt install ./ikafssn_<version>_ubuntu-<ubuntu_ver>_<arch>.deb
```

ダウンロードした `.rpm` のインストール:

```bash
sudo dnf install ./ikafssn-<version>.el<el_ver>.<arch>.rpm
```

### インストールの確認

```bash
ikafssnindex -version
```

## ソースからのビルド

### 依存関係

- C++20 コンパイラ (GCC >= 11, Clang >= 13)
- CMake >= 3.16
- NCBI C++ Toolkit (BLAST DB アクセス用)
- Intel TBB (並列化用)
- Parasail >= 2.6 (Stage 3 ペアワイズアライメント用)
- htslib >= 1.17 (SAM/BAM 出力用)
- zlib, libbz2, liblzma, libzstd (透過的な圧縮 I/O 用); libzstd の検出に pkg-config が必要
- Drogon (ikafssnhttpd 用、オプション)
- libcurl (HTTP クライアントモードおよびリモート取得用、オプション)

FastPFor (`.kix` / `.kpx` posting list のコーデック) と、x86_64 では x86-simd-sort は
CMake が configure 時にダウンロードします。CMake の FetchContent キャッシュが未作成の場合は
ビルドホストにネットワーク接続が必要です。これらを事前にインストールする必要はありません。

### CPU 要件

- **x86_64**: SSE4.2 必須 (Intel Nehalem 2008 年以降、AMD Bulldozer 2011 年以降)。実行時 SIMD ディスパッチャーは利用可能な場合に AVX2、AVX-512 BW、AVX-512 VBMI2 を選択します。SSE4.2 非対応の CPU では起動時に `exit(2)` で拒否されます。
- **aarch64**: NEON (ASIMD) 必須 (Armv8.0+)。NEON 非対応の aarch64 CPU は起動時に `exit(2)` で拒否されます。SVE / SVE2 対応 CPU でも FastPFor PForDelta codec は単一の NEON OBJECT library (FastPFor のネイティブ NEON ヘッダを使用) に降格して使用されます (ディスパッチャは tier ≥ NEON を一括して NEON tier object に流す)。各カーネル SIMD ファイル (toupper / ncbi2na unpack / k-mer revcomp / degenerate scan / spaced-seed) は SVE / SVE2 個別のディスパッチを保持します。

### 依存パッケージのインストール

NCBI C++ Toolkit 以外の依存パッケージを以下のコマンドでインストールできます。

**Ubuntu Server 22.04 / 24.04 / 26.04:**

```bash
sudo apt install build-essential cmake pkg-config libtbb-dev liblmdb-dev libsqlite3-dev \
    libcurl4-openssl-dev libjsoncpp-dev
sudo apt install zlib1g-dev libbz2-dev liblzma-dev libzstd-dev libdeflate-dev autoconf \
    libssl-dev uuid-dev
```

`pkg-config` は CMake の `pkg_check_modules` が libzstd を検出するのに必要です。2 行目は Parasail および htslib のソースビルドに必要な依存パッケージに加え、NCBI C++ Toolkit と Drogon のソースビルドに必要な `libssl-dev` と `uuid-dev` です。ikafssnhttpd が不要な場合は `-DBUILD_HTTPD=OFF` を指定してビルドし、`uuid-dev` は省略可能です (`libssl-dev` は NCBI C++ Toolkit が必要とします)。

**AlmaLinux 9 / Rocky Linux 9:**

```bash
sudo dnf config-manager --set-enabled crb
sudo dnf install -y epel-release
sudo dnf group install -y "Development Tools"
sudo dnf install -y cmake gcc-c++ pkgconfig tbb-devel lmdb-devel sqlite-devel \
    libcurl-devel jsoncpp-devel
sudo dnf install -y zlib-devel bzip2-devel xz-devel libzstd-devel libdeflate-devel autoconf
sudo dnf install -y libuuid-devel openssl-devel
```

**Oracle Linux 9:**

```bash
sudo dnf config-manager --set-enabled crb
sudo dnf install -y oracle-epel-release-el9
sudo dnf group install -y "Development Tools"
sudo dnf install -y cmake gcc-c++ pkgconfig tbb-devel lmdb-devel sqlite-devel \
    libcurl-devel jsoncpp-devel
sudo dnf install -y zlib-devel bzip2-devel xz-devel libzstd-devel libdeflate-devel autoconf
sudo dnf install -y libuuid-devel openssl-devel
```

**AlmaLinux 10 / Rocky Linux 10:**

```bash
sudo dnf config-manager --set-enabled crb
sudo dnf group install -y "Development Tools"
sudo dnf install -y cmake gcc-c++ pkgconfig tbb-devel lmdb-devel sqlite-devel \
    libcurl-devel jsoncpp-devel
sudo dnf install -y zlib-devel bzip2-devel xz-devel libzstd-devel libdeflate-devel autoconf
sudo dnf install -y libuuid-devel openssl-devel
```

**Oracle Linux 10:**

```bash
sudo dnf config-manager --set-enabled crb
sudo dnf group install -y "Development Tools"
sudo dnf install -y cmake gcc-c++ pkgconfig tbb-devel lmdb-devel sqlite-devel \
    libcurl-devel jsoncpp-devel
sudo dnf install -y zlib-devel bzip2-devel xz-devel libzstd-devel libdeflate-devel autoconf
sudo dnf install -y libuuid-devel openssl-devel
```

EL9 では `jsoncpp-devel` に EPEL、`lmdb-devel` に CRB リポジトリが必要です。EL10 ではいずれも CRB に収録されているため EPEL は不要です。`pkgconfig` は CMake の `pkg_check_modules` が libzstd を検出するのに必要です。各ブロックの最後から 2 行目は Parasail および htslib のソースビルドに必要な依存パッケージです。最終行は NCBI C++ Toolkit と Drogon のソースビルドに必要な依存パッケージです。ikafssnhttpd が不要な場合は `-DBUILD_HTTPD=OFF` を指定してビルドし、`libuuid-devel` は省略可能です (`openssl-devel` は NCBI C++ Toolkit が必要とします)。

**macOS (Homebrew):**

```bash
brew install cmake pkg-config tbb lmdb sqlite3 curl jsoncpp \
    xz zstd libdeflate autoconf automake libtool openssl@3
```

`openssl@3` は NCBI C++ Toolkit と Drogon のソースビルドに必要です (Drogon は別途ソースからビルドします、後述)。macOS では以下のビルド手順中の `make -j$(nproc)` を `make -j$(sysctl -n hw.ncpu)` に読み替えてください。

### Parasail

ikafssn は Stage 3 ペアワイズアライメントに Parasail ライブラリを使用します。デフォルトではソースルート直下の `./parasail` を参照します。別の場所にインストール済みの場合は `-DPARASAIL_DIR` で指定してください。

Parasail のダウンロード・ビルド・インストールは、ikafssn ソースルートで以下を実行します:

```bash
curl -L -o parasail-2.6.2.tar.gz \
    https://github.com/jeffdaily/parasail/archive/refs/tags/v2.6.2.tar.gz
tar xf parasail-2.6.2.tar.gz
cd parasail-2.6.2
patch -p1 < ../patches/parasail-degmatch-cigar-score.patch
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$(realpath ../..)/parasail" \
    -DBUILD_SHARED_LIBS=OFF
make -j$(nproc)
make install
cd ../..
```

macOS で CMake >= 4.0 を使用する場合、上記の `cmake` コマンドに `-DCMAKE_POLICY_VERSION_MINIMUM=3.5` を追加してください。

### htslib

ikafssn は SAM/BAM 出力に htslib を使用します。デフォルトではソースルート直下の `./htslib` を参照します。別の場所にインストール済みの場合は `-DHTSLIB_DIR` で指定してください。

htslib のダウンロード・ビルド・インストールは、ikafssn ソースルートで以下を実行します:

```bash
curl -L -o htslib-1.23.1.tar.bz2 \
    https://github.com/samtools/htslib/releases/download/1.23.1/htslib-1.23.1.tar.bz2
tar xf htslib-1.23.1.tar.bz2
cd htslib-1.23.1
autoreconf -i
./configure --prefix="$(realpath ..)/htslib" --disable-libcurl --disable-gcs --disable-s3
make -j$(nproc)
make install
cd ..
```

macOS では、configure が xz (lzma) や libdeflate を認識できるよう Homebrew のパスを指定してください:

```bash
./configure --prefix="$(realpath ..)/htslib" --disable-libcurl --disable-gcs --disable-s3 \
    CPPFLAGS="-I$(brew --prefix)/include" LDFLAGS="-L$(brew --prefix)/lib"
```

### NCBI C++ Toolkit

ikafssn のビルドには NCBI C++ Toolkit が必要です。デフォルトではソースルート直下の `./ncbi-cxx-toolkit` を参照します。別の場所にインストール済みの場合は `-DNCBI_TOOLKIT_DIR` で指定してください。

Toolkit 内のビルドサブディレクトリ名 (例: `CMake-GCC1330-Release`) はデフォルトで自動認識されますが、必要に応じて `-DNCBI_TOOLKIT_BUILD_TAG` で変更できます。

Toolkit のダウンロード・ビルド・インストールは、ikafssn ソースルートで以下を実行します:

```bash
curl -L -o ncbi-cxx-toolkit-public-release-30.2.0.tar.gz \
    https://github.com/ncbi/ncbi-cxx-toolkit-public/archive/refs/tags/release/30.2.0.tar.gz
tar xf ncbi-cxx-toolkit-public-release-30.2.0.tar.gz
cd ncbi-cxx-toolkit-public-release-30.2.0
patch -p1 < ../patches/ncbi-cxx-toolkit-seqdb-madvise-random.patch  # BLAST DB mmap のページキャッシュ汚染を防止
export CXXFLAGS="-std=c++20"
./cmake-configure \
    --without-debug \
    --with-projects="objtools/blast/seqdb_reader;objtools/blast/blastdb_format" \
    --with-install="$(realpath ..)/ncbi-cxx-toolkit"
cd CMake-GCC*/build
make -j$(nproc)
make install
cd ../../..
```

ikafssn が必要とするライブラリ (`seqdb`、`blastdb_format` およびその依存) のみをビルドします。Toolkit 全体のビルドは不要です。

macOS では、Homebrew の include パスをコンパイラに認識させる必要があります (`lmdb.h` 用)。また、ビルドディレクトリの glob パターンが異なります:

```bash
export CFLAGS="-I$(brew --prefix)/include"
export CXXFLAGS="-I$(brew --prefix)/include -std=c++20"
patch -p1 < ../patches/ncbi-cxx-toolkit-seqdb-madvise-random.patch  # BLAST DB mmap のページキャッシュ汚染を防止
./cmake-configure \
    --without-debug \
    --with-projects="objtools/blast/seqdb_reader;objtools/blast/blastdb_format" \
    --with-install="$(realpath ..)/ncbi-cxx-toolkit"
cd CMake-Clang*/build
make -j$(sysctl -n hw.ncpu)
make install
cd ../../..
```

### Drogon

ikafssnhttpd は Drogon HTTP フレームワークを使用します。デフォルトではソースルート直下の `./drogon` を参照します。別の場所にインストール済みの場合は `-DDROGON_DIR` で指定してください。`-DBUILD_HTTPD=OFF` を指定する場合、Drogon のビルドは不要です。

Drogon のリリース tarball には trantor サブモジュールの中身が含まれないため、対応する trantor のリリースを別途取得する必要があります (Drogon 1.9.12 が参照するコミットは trantor v1.5.26 です)。Drogon と trantor のダウンロード・ビルド・インストールは、ikafssn ソースルートで以下を実行します:

```bash
curl -L -o drogon-1.9.12.tar.gz \
    https://github.com/drogonframework/drogon/archive/refs/tags/v1.9.12.tar.gz
tar xf drogon-1.9.12.tar.gz
curl -L -o trantor-1.5.26.tar.gz \
    https://github.com/an-tao/trantor/archive/refs/tags/v1.5.26.tar.gz
tar xf trantor-1.5.26.tar.gz
rmdir drogon-1.9.12/trantor
mv trantor-1.5.26 drogon-1.9.12/trantor
cd drogon-1.9.12
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$(realpath ../..)/drogon" \
    -DBUILD_SHARED_LIBS=OFF \
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
    -DBUILD_CTL=OFF \
    -DBUILD_EXAMPLES=OFF \
    -DBUILD_ORM=OFF \
    -DBUILD_BROTLI=OFF \
    -DBUILD_YAML_CONFIG=OFF \
    -DBUILD_DOC=OFF
make -j$(nproc)
make install
cd ../..
```

ORM (PostgreSQL/MySQL/SQLite/Redis)、brotli、yaml-cpp、`drogon_ctl`、examples、ドキュメント生成を無効化することで、リンク時の依存ライブラリを jsoncpp、libuuid、zlib、OpenSSL の 4 つに最小化しています。

macOS では Homebrew の `openssl@3` が keg-only のため、`OPENSSL_ROOT_DIR` を明示する必要があります:

```bash
cmake .. -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$(realpath ../..)/drogon" \
    -DOPENSSL_ROOT_DIR="$(brew --prefix)/opt/openssl@3" \
    -DBUILD_SHARED_LIBS=OFF \
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
    -DBUILD_CTL=OFF \
    -DBUILD_EXAMPLES=OFF \
    -DBUILD_ORM=OFF \
    -DBUILD_BROTLI=OFF \
    -DBUILD_YAML_CONFIG=OFF \
    -DBUILD_DOC=OFF
make -j$(sysctl -n hw.ncpu)
make install
```

### ビルド

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
make test
```

NCBI C++ Toolkit がデフォルト以外の場所にインストールされている場合:

```bash
cmake .. -DCMAKE_BUILD_TYPE=Release \
    -DNCBI_TOOLKIT_DIR=/path/to/ncbi-cxx-toolkit
```

### インストール

```bash
sudo make install
```

デフォルトでは実行ファイルは `/usr/local/bin` に、systemd ユニットファイルは `/usr/local/share/ikafssn/systemd` にインストールされます。インストール先を変更する場合は `CMAKE_INSTALL_PREFIX` を指定します:

```bash
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/opt/ikafssn
make -j$(nproc)
sudo make install
```

この例では実行ファイルが `/opt/ikafssn/bin` にインストールされます。

### CMake オプション

| オプション | デフォルト | 説明 |
|---|---|---|
| `NCBI_TOOLKIT_DIR` | `${CMAKE_SOURCE_DIR}/ncbi-cxx-toolkit` | NCBI C++ Toolkit のインストールルートパス |
| `NCBI_TOOLKIT_BUILD_TAG` | `CMake-GCC1330-Release` | Toolkit ビルドサブディレクトリ名 |
| `PARASAIL_DIR` | `${CMAKE_SOURCE_DIR}/parasail` | Parasail のインストールルートパス |
| `HTSLIB_DIR` | `${CMAKE_SOURCE_DIR}/htslib` | htslib のインストールルートパス |
| `DROGON_DIR` | `${CMAKE_SOURCE_DIR}/drogon` | Drogon のインストールルートパス (`BUILD_HTTPD=ON` のとき使用) |
| `BUILD_HTTPD` | ON | ikafssnhttpd をビルド (Drogon が必要) |
| `BUILD_CLIENT` | ON | ikafssnclient をビルド (HTTP モードで libcurl が必要) |
| `ENABLE_REMOTE_RETRIEVE` | ON | ikafssnretrieve で NCBI efetch を有効化 |
