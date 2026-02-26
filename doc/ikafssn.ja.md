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

### ikafssnindex

BLAST DB から k-mer 転置インデックスを構築します。各ボリュームに対して `.kix` (ID ポスティング)、`.kpx` (位置ポスティング、`-mode 1` の場合は省略)、`.ksx` (配列メタデータ) のファイルを生成します。`-max_freq_build` 使用時は共有 `.khx` (構築時除外ビットセット) も生成されます。`.khx` ファイルは全ボリューム共通 (k 値ごとに 1 つ) です。

```
ikafssnindex [options]

必須:
  -db <path>              BLAST DB プレフィックス
  -k <int>                k-mer 長 (5〜13)
  -o <dir>                出力ディレクトリ

オプション:
  -mode <1|2|3>           インデックスがサポートする検索モード (デフォルト: 2)
                          1 = Stage 1 のみ (.kpx 生成スキップ、ディスク・時間節約)
                          2 = Stage 1+2 (デフォルト)
                          3 = Stage 1+2+3 (インデックス構築では 2 と同一動作)
  -memory_limit <size>    メモリ上限 (デフォルト: 物理メモリの半分)
                          接尾辞 K, M, G を認識
                          パーティション数はこの上限内に収まるよう自動決定
  -max_freq_build <num>   構築時高頻度 k-mer 除外閾値 (ボリューム横断集計)
                          1 以上: 絶対カウント閾値
                          0〜1 未満: 全ボリューム合計 NSEQ に対する割合
                          カウントは全ボリュームで合算後にフィルタリング
                          (デフォルト: 0 = 除外なし)
  -highfreq_filter_threads <int>
                          ボリューム横断高頻度フィルタリングのスレッド数
                          (デフォルト: min(8, threads))
  -openvol <int>          ボリューム同時処理数の上限 (デフォルト: 1)
                          マルチボリューム DB のピークメモリ使用量を制御
  -threads <int>          使用スレッド数 (デフォルト: 利用可能な全コア)
                          計数・パーティションスキャン・ソート・
                          ボリューム処理を並列化
  -v, --verbose           詳細ログ出力
```

**使用例:**

```bash
# 小規模 DB、メモリ豊富
ikafssnindex -db mydb -k 11 -o ./index -memory_limit 128G

# 大規模 DB、メモリ制限、マルチスレッド
ikafssnindex -db nt -k 11 -o ./nt_index -memory_limit 32G -threads 32

# 大規模 DB、ボリューム同時処理数を 2 に制限
ikafssnindex -db nt -k 11 -o ./nt_index -openvol 2

# 高頻度 k-mer を除外して構築 (絶対値指定)
ikafssnindex -db nt -k 11 -o ./nt_index -max_freq_build 50000

# 全ボリューム合計配列数の 1% を超える k-mer を除外
ikafssnindex -db nt -k 11 -o ./nt_index -max_freq_build 0.01

# mode 1 インデックスを構築 (Stage 1 のみ、.kpx なし)
ikafssnindex -db mydb -k 11 -o ./index -mode 1
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
  -threads <int>          並列検索スレッド数 (デフォルト: 利用可能な全コア)
  -mode <1|2|3>           検索モード (デフォルト: 2)
                          1=Stage 1 のみ、2=Stage 1+2、3=Stage 1+2+3
  -db <path>              モード 3 用 BLAST DB パス (デフォルト: -ix と同じ)
  -stage1_score <1|2>     Stage 1 スコア種別 (デフォルト: 1)
                          1=coverscore、2=matchscore
  -stage1_max_freq <num>  高頻度 k-mer スキップ閾値 (デフォルト: 0.5)
                          0〜1 未満: 全ボリューム合計 NSEQ に対する割合
                          1 以上: 絶対カウント閾値; 0 = 自動計算
  -stage1_topn <int>      Stage 1 候補数上限、0=無制限 (デフォルト: 0)
  -stage1_min_score <num> Stage 1 最小スコア閾値 (デフォルト: 0.5)
                          整数 (>= 1): 絶対閾値
                          小数 (0 < P < 1): クエリ k-mer に対する割合
                            クエリごとに ceil(Nqkmer * P) - Nhighfreq に解決
  -stage2_min_score <int> 最小チェインスコア (デフォルト: 0 = 適応的)
                          0 = Stage 1 の解決済み閾値を最小値として使用
                          1 以上: 絶対的な最小チェインスコア
  -stage2_max_gap <int>   チェイニング対角線ずれ許容幅 (デフォルト: 100)
  -stage2_max_lookback <int>  チェイニング DP 探索窓サイズ (デフォルト: 64、0=無制限)
  -stage2_min_diag_hits <int>  対角線フィルタ最小ヒット数 (デフォルト: 1)
  -context <value>        モード 3 のコンテクスト拡張 (デフォルト: 0)
                          整数: 拡張する塩基数; 小数: クエリ長に対する倍率
  -stage3_traceback <0|1> モード 3 でトレースバックを有効化 (デフォルト: 0)
  -stage3_gapopen <int>   モード 3 のギャップオープンペナルティ (デフォルト: 10)
  -stage3_gapext <int>    モード 3 のギャップ伸長ペナルティ (デフォルト: 1)
  -stage3_min_pident <num>  モード 3 の最小配列一致率フィルタ (デフォルト: 0)
  -stage3_min_nident <int>  モード 3 の最小一致塩基数フィルタ (デフォルト: 0)
  -stage3_fetch_threads <int>  モード 3 の BLAST DB 取得スレッド数 (デフォルト: min(8, threads); -threads を超えるとエラー)
  -num_results <int>      最終出力件数、0=無制限 (デフォルト: 0)
  -seqidlist <path>       検索対象を指定アクセッションに限定
  -negative_seqidlist <path>  指定アクセッションを検索対象から除外
  -strand <-1|1|2>       検索する鎖 (デフォルト: 2)
                          1=プラス鎖のみ、-1=マイナス鎖のみ、2=両鎖
  -accept_qdegen <0|1>    縮重塩基を含むクエリを許可 (デフォルト: 1)
  -outfmt <tab|json|sam|bam>  出力形式 (デフォルト: tab)
  -v, --verbose           詳細ログ出力
```

`-ix` オプションには `blastn -db` と同様にインデックスのプレフィックスパス (拡張子なし) を指定します。例えば、マルチボリューム BLAST DB (`nt`、ボリューム `nt.00`、`nt.01`) に対して `ikafssnindex -db nt -k 11 -o /data/index` が以下のファイルを生成した場合:

```
/data/index/nt.00.11mer.kix
/data/index/nt.00.11mer.kpx
/data/index/nt.00.11mer.ksx
/data/index/nt.01.11mer.kix
/data/index/nt.01.11mer.kpx
/data/index/nt.01.11mer.ksx
/data/index/nt.11mer.kvx
```

`-ix /data/index/nt` と指定します。プレフィックス `/data/index/nt` はディレクトリ `/data/index/` とベース名 `nt` に分割され、`.kvx` マニフェストファイル (`nt.11mer.kvx`) からボリュームベースネームの一覧を読み取ってボリュームを検出します。集約データベース (例: `combined` が `foo` と `bar` を束ねたもの) の場合、インデックスファイルは `foo.11mer.kix`、`bar.11mer.kix`、マニフェストは `combined.11mer.kvx` となります。

インデックスディレクトリに複数の k-mer サイズのインデックスが含まれる場合 (例: `nt.09mer.kvx` と `nt.11mer.kvx` の両方が存在する場合)、`-k` で使用するサイズを指定する必要があります。k-mer サイズが 1 種類のみの場合は `-k` を省略できます。

`-accept_qdegen` が 0 の場合、IUPAC 縮重塩基 (R, Y, S, W, K, M, B, D, H, V, N) を含むクエリは警告付きでスキップされ、終了コードは 2 になります。`-accept_qdegen 1` を指定すると縮重塩基を含むクエリも受け付けます。1 文字の縮重塩基を含む k-mer は全可能バリアントに展開して検索に使用されます (例: R→A,G で 2 k-mer、N→A,C,G,T で 4 k-mer)。2 文字以上の縮重塩基を含む k-mer はスキップされます。この場合、クエリごとに 1 回、当該クエリ名とスキップされた旨の警告が標準エラーに出力されます。サーバモード (`ikafssnserver`) ではこの警告がプロトコル経由で `ikafssnclient` に伝播され、クライアント側でも同じメッセージが表示されます。この処理はインデックス構築時のサブジェクト配列に対する処理と同等です。

`-seqidlist` と `-negative_seqidlist` は排他的 (同時指定不可) です。ファイル形式はテキスト (1 行 1 アクセッション) と `blastdb_aliastool -seqid_file_in` で生成されるバイナリ形式の両方を受け付け、先頭のマジックバイトで自動判別します。

**使用例:**

```bash
# 基本的な検索
ikafssnsearch -ix ./index/mydb -query query.fasta -threads 8

# k-mer サイズを指定 (インデックスに複数の k 値が含まれる場合は必須)
ikafssnsearch -ix ./index/mydb -k 11 -query query.fasta

# 感度を上げた検索
ikafssnsearch -ix ./index/mydb -query query.fasta \
    -stage2_min_score 2 -stage1_topn 2000 -stage1_max_freq 50000

# seqidlist で検索対象を限定
ikafssnsearch -ix ./index/mydb -query query.fasta -seqidlist targets.txt

# negative_seqidlist で特定配列を除外
ikafssnsearch -ix ./index/mydb -query query.fasta -negative_seqidlist exclude.txt

# 割合指定の Stage 1 閾値 (クエリ k-mer の 50%)
ikafssnsearch -ix ./index/mydb -query query.fasta -stage1_min_score 0.5

# モード 3: トレースバック付きペアワイズアライメント
ikafssnsearch -ix ./index/mydb -query query.fasta -mode 3 -stage3_traceback 1 -num_results 5

# モード 3: SAM 出力
ikafssnsearch -ix ./index/mydb -query query.fasta -mode 3 -stage3_traceback 1 -outfmt sam -o result.sam

# モード 3: BAM 出力
ikafssnsearch -ix ./index/mydb -query query.fasta -mode 3 -stage3_traceback 1 -outfmt bam -o result.bam

# モード 3: 配列一致率でフィルタ
ikafssnsearch -ix ./index/mydb -query query.fasta -mode 3 -stage3_traceback 1 -stage3_min_pident 90

# モード 3: コンテクスト拡張 (前後各50塩基)
ikafssnsearch -ix ./index/mydb -query query.fasta -mode 3 -context 50 -num_results 5

# モード 3: コンテクスト拡張 (前後各クエリ長の0.1倍)
ikafssnsearch -ix ./index/mydb -query query.fasta -mode 3 -context 0.1 -num_results 5

# パイプラインで ikafssnretrieve に接続
ikafssnsearch -ix ./index/mydb -query query.fasta | ikafssnretrieve -db nt > matches.fasta
```

### ikafssnretrieve

検索結果に基づきマッチした部分配列を抽出します。配列ソースとしてローカル BLAST DB または NCBI E-utilities (efetch) を選択できます。

```
ikafssnretrieve [options]

配列ソース (いずれか必須):
  -db <path>              ローカル BLAST DB プレフィックス
  -remote                 NCBI efetch からリモート取得

入力:
  -results <path>         検索結果ファイル (tab 形式)
  (指定なし)              標準入力から読み込み

共通オプション:
  -o <path>               出力 FASTA ファイル (デフォルト: 標準出力)
  -context <value>        コンテクスト拡張 (デフォルト: 0)
                          整数: マッチ領域の前後に付加する塩基数
                          小数: クエリ長 (q_len) に対する倍率
  -v, --verbose           詳細ログ出力

リモート取得オプション (-remote 時):
  -api_key <key>          NCBI API key (環境変数 NCBI_API_KEY でも可)
  -batch_size <int>       バッチあたりのアクセッション数 (デフォルト: 100)
  -retries <int>          リトライ回数 (デフォルト: 3)
  -timeout <int>          リクエストタイムアウト秒数 (デフォルト: 30)
  -range_threshold <int>  部分取得に切り替える配列長閾値 (デフォルト: 100000)
```

**使用例:**

```bash
# ローカル BLAST DB から抽出 (ファイル入力)
ikafssnsearch -ix ./index/mydb -query query.fasta -o results.tsv
ikafssnretrieve -db nt -results results.tsv -o matches.fasta

# ローカル BLAST DB から抽出 (パイプ)
ikafssnsearch -ix ./index/mydb -query query.fasta | ikafssnretrieve -db nt > matches.fasta

# サーバ経由の検索結果からも同様に利用可能
ikafssnclient -socket /var/run/ikafssn.sock -db nt -query query.fasta | ikafssnretrieve -db nt > matches.fasta

# NCBI efetch からリモート取得
ikafssnsearch -ix ./index/mydb -query query.fasta | ikafssnretrieve -remote > matches.fasta

# NCBI efetch + API key (高スループット)
ikafssnclient -http http://search.example.com:8080 -db nt -query query.fasta \
    | ikafssnretrieve -remote -api_key XXXXXXXX > matches.fasta

# マッチ領域の前後 50bp を含めて抽出
ikafssnretrieve -db nt -results results.tsv -context 50

# コンテクストをクエリ長の倍率で指定 (前後各 0.1 倍)
ikafssnretrieve -db nt -results results.tsv -context 0.1
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

オプション:
  -threads <int>          ワーカースレッド数 (デフォルト: 利用可能な全コア)
  -max_query <int>        同時処理クエリ配列数のグローバル上限 (デフォルト: 1024)
  -max_seqs_per_req <int> 1 リクエストあたりの受理配列数上限 (デフォルト: スレッド数)
  -pid <path>             PID ファイルパス
  -db <path>              モード 3 用 BLAST DB パス (繰り返し指定、-ix と対応;
                          デフォルト: 対応する -ix プレフィックスと同じ)
  -stage1_max_freq <num>  デフォルト高頻度 k-mer スキップ閾値 (デフォルト: 0.5)
                          0〜1 未満: 全ボリューム合計 NSEQ に対する割合
                          1 以上: 絶対カウント閾値; 0 = 自動計算
  -stage1_topn <int>      デフォルト Stage 1 候補数上限 (デフォルト: 0)
  -stage1_min_score <num> デフォルト Stage 1 最小スコア閾値 (デフォルト: 0.5)
                          整数 (>= 1) または小数 (0 < P < 1)
  -stage2_min_score <int> デフォルト最小チェインスコア (デフォルト: 0 = 適応的)
  -stage2_max_gap <int>   デフォルトチェイニング対角線ずれ許容幅 (デフォルト: 100)
  -stage2_max_lookback <int>  デフォルトチェイニング DP 探索窓サイズ (デフォルト: 64、0=無制限)
  -stage2_min_diag_hits <int> デフォルト対角線フィルタ最小ヒット数 (デフォルト: 1)
  -context <value>        デフォルトコンテクスト拡張 (デフォルト: 0)
                          整数: 拡張する塩基数; 小数: クエリ長に対する倍率
  -stage3_traceback <0|1> デフォルトトレースバックモード (デフォルト: 0)
  -stage3_gapopen <int>   デフォルトギャップオープンペナルティ (デフォルト: 10)
  -stage3_gapext <int>    デフォルトギャップ伸長ペナルティ (デフォルト: 1)
  -stage3_min_pident <num>  デフォルト最小配列一致率 (デフォルト: 0)
  -stage3_min_nident <int>  デフォルト最小一致塩基数 (デフォルト: 0)
  -stage3_fetch_threads <int>  BLAST DB 取得スレッド数 (デフォルト: min(8, threads))
  -num_results <int>      デフォルト最終出力件数 (デフォルト: 0)
  -accept_qdegen <0|1>    デフォルト縮重塩基クエリ許可 (デフォルト: 1)
  -shutdown_timeout <int> グレースフルシャットダウンのタイムアウト秒数 (デフォルト: 30)
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

# モード 3 対応: BLAST DB と Stage 3 デフォルトを指定
ikafssnserver -ix ./nt_index -db nt -socket /var/run/ikafssn.sock \
    -stage3_traceback 1 -context 50

# マルチ DB: 1 プロセスで 2 つのデータベースを同時サーブ
ikafssnserver -ix ./nt_index -db nt -ix ./rs_index -db refseq_genomic \
    -socket /var/run/ikafssn.sock -threads 32
```

**運用上の特性:**

- 1 プロセスで複数の BLAST DB インデックスを同時にサーブできます。`-ix` (および必要に応じて `-db`) を複数回指定して複数データベースをロードします。各データベースは `-ix` プレフィックスのベースネーム (パスの最終コンポーネント) で識別され、サーバが複数 DB をホストする場合、クライアントは `-db <name>` でターゲット DB を指定する必要があります。
- `-db` を指定する場合、その数は `-ix` の数と一致する必要があります (順番に対応)。`-db` を省略した DB は `-ix` プレフィックスを BLAST DB パスとして使用します。`-db` パスが未指定の DB はモード 1-2 のみ対応 (max_mode=2)、`-db` を指定するとモード 3 も利用可能 (max_mode=3) になります。
- `-ix` プレフィックスに対応する異なる k-mer サイズのインデックスが存在する場合、全て読み込み、クライアントのリクエストで k を指定できます。
- SIGTERM/SIGINT 受信時はグレースフルシャットダウンを行います。新規接続の受付を停止し、実行中のリクエストの完了を最大 `-shutdown_timeout` 秒待ちます。
- **配列単位の同時実行制御:** サーバは接続単位ではなく、配列単位で同時実行数を制御します。リクエストが到着すると、有効なクエリ配列ごとにパーミットの取得を試みます。グローバル上限 (`-max_query`) に達した場合、超過分の配列はリトライ用に「拒否」としてクライアントに返されます。`-max_seqs_per_req` は 1 リクエストが取得できるパーミット数の上限を設定し、大量配列を含む単一リクエストによるスロットの独占を防ぎます。

### ikafssnhttpd

HTTP REST API デーモンです。1 つ以上の `ikafssnserver` インスタンスに接続し、HTTP REST API として検索サービスを提供します。Drogon フレームワークを使用します。複数のバックエンドを指定することで、マルチデータベース対応や同一 DB のレプリカ間の負荷分散が可能です。

起動時にすべてのバックエンドに接続してサーバ能力情報をキャッシュします (最大 30 秒間、指数バックオフでリトライ)。同じデータベース名が複数のバックエンドに存在する場合、k 値セット・合計配列数・合計塩基数の一致をクロスサーバ検証し、不一致があれば起動エラーとなります。検索リクエストはまず統合された能力情報に対して同期的にチェック (バックエンドへの通信なし) し、明らかに無効なリクエストを即座に拒否します。その後、優先度と空きスロット状況に基づいて最適なバックエンドにルーティングします。

**ルーティングとヘルスチェック:**

- **優先度**: バックエンドは CLI 引数の順序で優先度が決まります (先 = 最高優先度)。
- **選択**: 検索リクエストごとに、最高優先度かつ有効な空き容量があるバックエンドを選択します。有効な空き容量は、スロット空き数 (`max_active_sequences - active_sequences`) とリクエストあたり上限 (`max_seqs_per_req`) の小さい方で決まります。全バックエンドが満杯の場合は、最高優先度のバックエンドが使用されます。
- **事前チェック**: 検索前に選択されたバックエンドへ最新の info リクエストを送信し、接続を確認します。
- **除外**: バックエンドが応答しない場合 (info/search の接続エラー)、`-exclusion_time` 秒間除外されます。除外されたバックエンドはハートビート時に自動的に再チェックされ、到達可能になれば再有効化されます。
- **ハートビート**: バックグラウンドスレッドが `-heartbeat_interval` 秒ごとにすべてのバックエンドの info を更新します。
- **httpd はリトライしない**: 検索リクエストがバックエンド選択後に失敗した場合、エラーはクライアントに返されます。拒否されたクエリのリトライは `ikafssnclient` が行います。

```
ikafssnhttpd [options]

バックエンド接続 (1 つ以上必須、指定順 = 優先度):
  -server_socket <path>      ikafssnserver の UNIX ソケットパス
  -server_tcp <host>:<port>  ikafssnserver の TCP アドレス

オプション:
  -listen <host>:<port>       HTTP リスニングアドレス (デフォルト: 0.0.0.0:8080)
  -path_prefix <prefix>       API パスプレフィックス (例: /nt)
  -threads <int>              Drogon I/O スレッド数 (デフォルト: 利用可能な全コア)
  -heartbeat_interval <int>   ハートビート間隔 (秒、デフォルト: 3600)
  -exclusion_time <int>       バックエンド除外時間 (秒、デフォルト: 3600)
  -pid <path>                 PID ファイルパス
  -v, --verbose               詳細ログ出力
```

**REST API エンドポイント:**

| メソッド | パス | 説明 |
|---|---|---|
| POST | `/api/v1/search` | 検索リクエスト (クエリ配列を JSON ボディで送信) |
| GET | `/api/v1/info` | 全バックエンドの統合インデックス情報 |
| GET | `/api/v1/health` | ヘルスチェック (いずれかのバックエンドが到達可能なら OK) |

`/api/v1/info` レスポンスは全 healthy バックエンドのデータベースを統合して返します。複数バックエンドで提供されるデータベースの場合、容量情報は kmer_group 内の `modes` 配列にモードごとに集約され、全提供バックエンドの `max_active_sequences`、`active_sequences`、および `max_seqs_per_req` (各バックエンドの `min(空きスロット, per_req)` の合計として算出) が表示されます。トップレベルにも全モードの最小値として `max_seqs_per_req` フィールドが出力されます。

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

クライアントコマンドです。`ikafssnserver` に直接ソケット接続するか、`ikafssnhttpd` に HTTP 接続して検索結果を取得します。出力形式は `ikafssnsearch` と同一です。クエリ送信前にサーバの能力情報を取得し、指定されたデータベース名・k-mer サイズ・モードの妥当性を事前検証 (プリフライトチェック) します。無効なパラメータが指定された場合は、クエリデータの送信前に利用可能なデータベース一覧を含むエラーメッセージが表示されます。クライアントはサーバの `max_seqs_per_req` と空きスロット数に基づいてクエリを適切なサイズのバッチに自動分割し、部分的に拒否されるような過大なリクエストを回避します。各バッチ内でサーバが同時実行制限によりクエリ配列を拒否した場合、拒否された配列を指数バックオフ (30 秒、60 秒、120 秒、120 秒、…) で自動リトライし、全配列の処理が完了するまで繰り返します。

```
ikafssnclient [options]

接続先 (いずれか):
  -socket <path>           ikafssnserver の UNIX ソケットパス
  -tcp <host>:<port>       ikafssnserver の TCP アドレス
  -http <url>              ikafssnhttpd の URL (例: http://example.com:8080)

必須:
  -query <path>            クエリ FASTA ファイル (- で標準入力)

オプション:
  -o <path>                出力ファイル (デフォルト: 標準出力)
  -k <int>                 使用する k-mer サイズ (デフォルト: サーバ側デフォルト)
  -mode <1|2|3>            検索モード (デフォルト: サーバ側デフォルト)
  -stage1_score <1|2>      Stage 1 スコア種別 (デフォルト: サーバ側デフォルト)
  -stage1_max_freq <num>   高頻度 k-mer スキップ閾値 (デフォルト: サーバ側デフォルト)
                           0〜1 未満: 全ボリューム合計 NSEQ に対する割合
                           1 以上: 絶対カウント閾値
  -stage1_topn <int>       Stage 1 候補数上限 (デフォルト: サーバ側デフォルト)
  -stage1_min_score <num>  Stage 1 最小スコア閾値 (デフォルト: サーバ側デフォルト)
                           整数 (>= 1) または小数 (0 < P < 1)
  -stage2_min_score <int>  最小チェインスコア (デフォルト: サーバ側デフォルト)
                           0 = 適応的モードを明示的にリクエスト
  -stage2_max_gap <int>    チェイニング対角線ずれ許容幅 (デフォルト: サーバ側デフォルト)
  -stage2_max_lookback <int>  チェイニング DP 探索窓サイズ (デフォルト: サーバ側デフォルト)
  -stage2_min_diag_hits <int> 対角線フィルタ最小ヒット数 (デフォルト: サーバ側デフォルト)
  -context <value>         コンテクスト拡張 (デフォルト: サーバ側デフォルト)
                           整数: 拡張する塩基数; 小数: クエリ長に対する倍率
  -stage3_traceback <0|1>  トレースバック有効化 (デフォルト: サーバ側デフォルト)
  -stage3_gapopen <int>    ギャップオープンペナルティ (デフォルト: サーバ側デフォルト)
  -stage3_gapext <int>     ギャップ伸長ペナルティ (デフォルト: サーバ側デフォルト)
  -stage3_min_pident <num> 最小配列一致率フィルタ (デフォルト: サーバ側デフォルト)
  -stage3_min_nident <int> 最小一致塩基数フィルタ (デフォルト: サーバ側デフォルト)
  -num_results <int>       最終出力件数 (デフォルト: サーバ側デフォルト)
  -db <name>               サーバ上のターゲットデータベース名 (マルチ DB サーバでは必須)
  -seqidlist <path>        検索対象を指定アクセッションに限定
  -negative_seqidlist <path>  指定アクセッションを検索対象から除外
  -strand <-1|1|2>         検索する鎖: 1=プラス、-1=マイナス、2=両鎖 (デフォルト: サーバ側デフォルト)
  -accept_qdegen <0|1>     縮重塩基を含むクエリを許可 (デフォルト: 1)
  -outfmt <tab|json|sam|bam>  出力形式 (デフォルト: tab)
  -v, --verbose            詳細ログ出力

HTTP 認証 (HTTP モード専用):
  --user <user:password>   認証情報 (curl 形式)
  --http-user <USER>       ユーザー名 (wget 形式)
  --http-password <PASS>   パスワード (--http-user と併用)
  --netrc-file <path>      .netrc ファイルのパス
```

**使用例:**

```bash
# UNIX ソケット経由 (ローカル)
ikafssnclient -socket /var/run/ikafssn.sock -db nt -query query.fasta

# TCP 直接接続
ikafssnclient -tcp 10.0.1.5:9100 -db nt -query query.fasta

# HTTP 経由
ikafssnclient -http http://search.example.com:8080 -db nt -query query.fasta

# seqidlist で検索対象を限定
ikafssnclient -socket /var/run/ikafssn.sock -db nt -query query.fasta -seqidlist targets.txt

# パイプラインで ikafssnretrieve に接続
ikafssnclient -socket /var/run/ikafssn.sock -db nt -query query.fasta | ikafssnretrieve -db nt > matches.fasta

# 特定の k-mer サイズを指定
ikafssnclient -socket /var/run/ikafssn.sock -db nt -query query.fasta -k 9

# HTTP Basic 認証 (curl 形式)
ikafssnclient -http http://search.example.com:8080 -db nt -query query.fasta --user admin:secret

# HTTP Basic 認証 (wget 形式)
ikafssnclient -http http://search.example.com:8080 -db nt -query query.fasta --http-user=admin --http-password=secret

# .netrc ファイルによる認証
ikafssnclient -http http://search.example.com:8080 -db nt -query query.fasta --netrc-file ~/.netrc

# モード 3: トレースバック付きペアワイズアライメント
ikafssnclient -socket /var/run/ikafssn.sock -db nt -query query.fasta -mode 3 -stage3_traceback 1

# モード 3: SAM 出力
ikafssnclient -socket /var/run/ikafssn.sock -db nt -query query.fasta -mode 3 -stage3_traceback 1 -outfmt sam -o result.sam
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

ローカルモードオプション:
  -db <path>               BLAST DB プレフィックス (デフォルト: -ix から自動検出)

リモート HTTP 認証:
  --user <user:password>   認証情報 (curl 形式)
  --http-user <USER>       ユーザー名 (wget 形式)
  --http-password <PASS>   パスワード (--http-user と併用)
  --netrc-file <path>      .netrc ファイルのパス

オプション:
  -v, --verbose            詳細出力
```

`-ix` とリモートオプション (`-socket`、`-tcp`、`-http`) は排他的です。リモートオプションは同時に 1 つのみ指定可能です。

**ローカルモード:** インデックスファイルを直接読み取り、詳細な統計情報を表示します。`-db` 未指定の場合、インデックスプレフィックスパスが有効な BLAST DB に対応するかを確認し、自動検出を試みます。

ローカルモードの出力情報:

- k-mer 長 (k) および k-mer 整数型 (uint16/uint32)
- ボリューム数
- 各ボリュームの配列数、総ポスティング数、ファイルサイズ、除外 k-mer 数 (`.khx` 存在時)
- 全体統計: 総配列数、総ポスティング数、総インデックスサイズ、圧縮率
- `-v` 指定時: k-mer 出現頻度分布 (min, max, mean, パーセンタイル)
- `-db` 指定時 (または自動検出時): BLAST DB のタイトル、配列数、総塩基数、ボリューム構成

**リモートモード:** 稼働中のサーバに問い合わせ、サーバの能力情報を表示します。

リモートモードの出力情報:

- アクティブ / 最大配列スロット数
- データベースごとの情報: 名前、デフォルト k、最大モード、k-mer グループ (ボリューム数・配列数・総塩基数・ポスティング数の統計)
- `-v` 指定時: 各 k-mer グループ内のボリュームごとの詳細 (配列数、総塩基数、ポスティング数)

**使用例:**

```bash
# ローカル: 基本的なインデックス情報
ikafssninfo -ix ./index/mydb

# ローカル: BLAST DB 情報も表示
ikafssninfo -ix ./index/mydb -db mydb

# ローカル: 詳細な頻度分布を表示
ikafssninfo -ix ./index/mydb -v

# リモート: UNIX ソケット経由でサーバに問い合わせ
ikafssninfo -socket /var/run/ikafssn.sock

# リモート: TCP 経由でサーバに問い合わせ
ikafssninfo -tcp 10.0.1.5:9100

# リモート: HTTP 経由でサーバに問い合わせ
ikafssninfo -http http://search.example.com:8080

# リモート: 詳細出力 (ボリュームごとの情報を表示)
ikafssninfo -socket /var/run/ikafssn.sock -v

# リモート: HTTP 認証付き
ikafssninfo -http http://search.example.com:8080 --user admin:secret
```

## 検索パイプライン

ikafssn は 3 段階の検索パイプラインを使用します。

デフォルトパラメータはスループットを優先しています。`stage1_topn=0` と `num_results=0` によりソートを省略し、`stage1_min_score=0.5` (割合指定) でクエリ k-mer の 50% 以上のマッチを要求してフィルタリングします。ランク付けされた出力が必要な場合は `-stage1_topn` や `-num_results` に正の値を設定してください。ソートが有効になりますが、結果件数が多い場合は速度が低下する可能性があります。

1. **Stage 1 (候補選択):** クエリの各 k-mer に対して ID ポスティングをスキャンし、配列ごとにスコアを集計します。スコア種別は 2 種類あります: **coverscore** (配列にマッチしたクエリ k-mer の種類数) と **matchscore** (クエリ k-mer と参照配列位置の総マッチ数)。`stage1_min_score` 以上のスコアを持つ配列を候補として選出します。`stage1_topn > 0` の場合はスコア順にソートして切り詰めます。`stage1_topn = 0` (デフォルト) の場合は全候補をソートせずに返します。

2. **Stage 2 (コリニアチェイニング):** 各候補に対して `.kpx` から位置レベルのヒットを収集し、対角線フィルタを適用した後、チェイニング DP により最良のコリニアチェインを求めます。チェインの長さが **chainscore** として報告されます。`chainscore >= stage2_min_score` のチェインが結果に含まれます。DP の内側ループは `-stage2_max_lookback` (デフォルト: 64) で制限され、各ヒットは直前の B 個のヒットのみを前駆候補として参照します。これにより、単一クエリ×サブジェクト間のヒット数が非常に多い場合の最悪計算量を O(n²) から O(n×B) に削減します。0 を指定すると無制限 (従来の O(n²) 動作) になります。

3. **Stage 3 (ペアワイズアライメント):** Stage 2 の各ヒットに対して、BLAST DB からサブジェクト部分配列を取得し (`-context` による拡張オプション付き)、Parasail ライブラリ (nuc44 スコアリングマトリクス) を使って半大域ペアワイズアライメントを実行します。全ヒットに対してアライメントスコア (**alnscore**) が計算されます。`-stage3_traceback 1` を指定すると、CIGAR 文字列、配列一致率、一致塩基数、不一致数、ギャップ付きアライメント配列も計算されます。`-stage3_min_pident` と `-stage3_min_nident` によるフィルタリングが可能です (トレースバックモードのみ)。サブジェクト配列は `-stage3_fetch_threads` で制御されるボリューム並列プリフェッチで取得されます。

**適応的 `-stage2_min_score` (デフォルト):** `-stage2_min_score 0` (デフォルト) の場合、最小チェインスコアはクエリごとに適応的に設定され、解決済みの Stage 1 閾値が使用されます。割合指定の `-stage1_min_score` (例: `0.5`) との組み合わせでは、各クエリの k-mer 構成に基づくクエリごとの適応的閾値が設定されます。絶対値指定の `-stage1_min_score` の場合は、その設定値がそのまま使用されます。固定閾値を使用する場合は `-stage2_min_score` に正の整数を指定してください。

**Mode 1 (Stage 1 のみ):** `-mode 1` を指定すると Stage 2, 3 が省略されます。`.kpx` ファイルへのアクセスが不要となり、I/O とメモリを節約できます。結果には Stage 1 スコアのみが含まれ、位置フィールド (q_start, q_end, s_start, s_end) と chainscore は省略されます。ソート基準は Stage 1 スコアに強制されます。

**Mode 3 (全パイプライン):** `-mode 3` を指定すると全 3 段階が実行されます。BLAST DB が必要です (`-db` で指定、デフォルトはインデックスプレフィックスと同じ)。ソート基準は alnscore に自動設定されます。SAM/BAM 出力には `-mode 3` と `-stage3_traceback 1` の両方が必要です。

デフォルトではクエリのフォワード鎖とリバースコンプリメント鎖の両方を検索します。`-strand 1` でプラス (フォワード) 鎖のみ、`-strand -1` でマイナス (リバースコンプリメント) 鎖のみの検索に制限できます。

### 高頻度 k-mer フィルタリング

高頻度 k-mer フィルタリングはボリューム単位のループに入る前に全ボリューム横断でグローバルに行われます。全ボリュームの k-mer カウントを合算し、`stage1_max_freq` を超える k-mer はクエリから一度だけ除去されます。これにより、データのボリューム分割方法に関わらず一貫したフィルタリングが保証されます。ビルド時除外情報 (`.khx`) もグローバルにチェックされます。

`-stage1_max_freq` のデフォルト値は `0.5` で、全ボリューム合計配列数の 50% を超えて出現する k-mer がスキップされます。より一般に、小数値 (0 < x < 1) を指定すると、閾値は `ceil(x * total_NSEQ)` に解決されます (total_NSEQ は全ボリュームの配列数の合計)。整数値 (1 以上) はそのまま絶対カウント閾値として使用されます。

`-stage1_max_freq 0` を明示的に指定すると、ボリュームごとに以下の式で自動計算されます:

```
max_freq = mean_count * 10    ([1000, 100000] に制限)
ここで mean_count = total_postings / 4^k
```

この自動計算モードではボリュームごとに `.kix` ヘッダから値が算出されます。

**構築時除外** (`-max_freq_build`): `-max_freq_build` を指定してインデックスを構築すると、高頻度 k-mer がインデックスから完全に除外されます。k-mer カウントは全ボリュームで合算された後に閾値と比較されるため、各ボリュームでは閾値未満だが合計では閾値を超える k-mer も正しく除外されます。除外された k-mer は共有 `.khx` ファイル (k 値ごとに 1 つ、ボリュームごとではない) に記録されます。小数値 (0 < x < 1) を指定した場合、閾値は全ボリューム合計 NSEQ に基づいて解決されます (`-stage1_max_freq` と同じ方式)。検索時に割合指定の `-stage1_min_score` を使用する場合、構築時に除外された k-mer が `.khx` ファイルから認識され、閾値計算から差し引かれます。

### 割合指定の Stage 1 閾値

`-stage1_min_score` を小数 (0 < P < 1) で指定すると、閾値はクエリごとに以下の式で解決されます:

```
threshold = ceil(Nqkmer * P) - Nhighfreq
```

各変数の意味:
- **Nqkmer**: クエリ k-mer の数 (coverscore では種類数、matchscore では総位置数)
- **Nhighfreq**: 除外されるクエリ k-mer の数。以下を合算:
  - 検索時除外: カウント > `max_freq` の k-mer
  - 構築時除外: `.khx` に記録された k-mer (存在する場合)

解決された閾値が 0 以下の場合、その鎖は警告付きでスキップされます。

### スコア種別

ikafssn は 3 種類のスコアを計算します。

| スコア | 説明 | 計算ステージ |
|---|---|---|
| **coverscore** | 参照配列にマッチしたクエリ k-mer の種類数。各クエリ k-mer は参照配列あたり最大 1 回カウントされます (複数位置にマッチしても重複計上されません)。 | Stage 1 |
| **matchscore** | (クエリ k-mer, 参照配列位置) の総マッチ数。1 つのクエリ k-mer が参照配列の複数位置にマッチした場合、その分だけ加算されます。 | Stage 1 |
| **chainscore** | チェイニング DP が求めた最良コリニアチェインの長さ (k-mer ヒット数)。`.kpx` の位置データを使用します。 | Stage 2 |
| **alnscore** | Parasail による半大域ペアワイズアライメントスコア (nuc44 マトリクス)。BLAST DB からのサブジェクト配列取得が必要です。 | Stage 3 |

- `-stage1_score` で Stage 1 が使用するスコア種別を選択します (1=coverscore, 2=matchscore)。候補のランキングと出力される Stage 1 スコアに影響します。
- ソート基準はモードにより自動決定: mode 1 は Stage 1 スコア、mode 2 は chainscore、mode 3 は alnscore。
- `-mode 1` では Stage 1 スコアのみが利用可能で、chainscore と alnscore は計算されません。

## 出力形式

### Tab 形式 (デフォルト)

**Mode 2** (デフォルト):

タブ区切りのカラムです。`-stage1_score 2` のとき `coverscore` が `matchscore` に置き換わります。

```
# query_id  accession  strand  q_start  q_end  q_len  s_start  s_end  s_len  coverscore  chainscore  volume
```

**Mode 1** (`-mode 1`):

```
# query_id  accession  strand  q_len  s_len  coverscore  volume
```

**Mode 3, traceback=0** (`-mode 3`):

```
# query_id  accession  strand  q_end  q_len  s_end  s_len  coverscore  chainscore  alnscore  volume
```

注: トレースバックなしでは正確なアライメント開始位置が得られないため、`q_start` と `s_start` は省略されます。

**Mode 3, traceback=1** (`-mode 3 -stage3_traceback 1`):

```
# query_id  accession  strand  q_start  q_end  q_len  s_start  s_end  s_len  coverscore  chainscore  alnscore  pident  nident  nmismatch  cigar  q_seq  s_seq  volume
```

### JSON 形式

**Mode 2** (デフォルト):

```json
{
  "results": [
    {
      "query_id": "query1",
      "hits": [
        {
          "accession": "NC_001234.5",
          "strand": "+",
          "q_start": 0,
          "q_end": 150,
          "q_len": 200,
          "s_start": 1000,
          "s_end": 1150,
          "s_len": 5000,
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
      "query_id": "query1",
      "hits": [
        {
          "accession": "NC_001234.5",
          "strand": "+",
          "q_len": 200,
          "s_len": 5000,
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
      "query_id": "query1",
      "hits": [
        {
          "accession": "NC_001234.5",
          "strand": "+",
          "q_end": 150,
          "q_len": 200,
          "s_end": 1150,
          "s_len": 5000,
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
      "query_id": "query1",
      "hits": [
        {
          "accession": "NC_001234.5",
          "strand": "+",
          "q_start": 0,
          "q_end": 150,
          "q_len": 200,
          "s_start": 1000,
          "s_end": 1150,
          "s_len": 5000,
          "coverscore": 8,
          "chainscore": 12,
          "alnscore": 240,
          "pident": 95.3,
          "nident": 143,
          "nmismatch": 7,
          "cigar": "50=2X48=1I50=",
          "q_seq": "ACGT...",
          "s_seq": "ACGT...",
          "volume": 0
        }
      ]
    }
  ]
}
```

### SAM/BAM 形式

SAM/BAM 出力には `-mode 3 -stage3_traceback 1` が必要です。`-outfmt sam` で SAM、`-outfmt bam` で BAM を出力します (BAM は `-o <path>` が必須)。

SAM レコードの構成:
- **QNAME**: query_id
- **FLAG**: 0 (フォワード) / 16 (リバース)
- **RNAME**: accession
- **POS**: s_start + 1 (1-based)
- **MAPQ**: 255
- **CIGAR**: 拡張 CIGAR (=/X/I/D 演算子)
- **SEQ**: ギャップなしクエリ配列
- **QUAL**: * (利用不可)
- **タグ**: `AS:i` (alnscore), `NM:i` (nmismatch), `cs:i` (chainscore), `s1:i` (stage1_score)

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

ikafssnclient -socket /var/run/ikafssn.sock -db nt -query query.fasta
ikafssnclient -socket /var/run/ikafssn.sock -db refseq_genomic -query query.fasta
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

同じデータベース名が複数バックエンドに存在する場合、`ikafssnhttpd` は起動時に k 値セット・合計配列数・合計塩基数の一致を検証します。リクエストは有効な空き容量 (スロット空き数と `max_seqs_per_req` の小さい方) がある最高優先度のバックエンドにルーティングされます。容量値 (`max_active_sequences`、`active_sequences`、`max_seqs_per_req`) はサーバ単位で共有されるため、同一サーバ上の複数データベース間で共有される点に注意してください。

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

BLAST DB ボリュームごとに、ボリューム自身のベースネームを使って 3 つのファイルが生成されます:

```
<vol_basename>.<kk>mer.kix   — ID ポスティング (直接アドレステーブル + デルタ圧縮)
<vol_basename>.<kk>mer.kpx   — 位置ポスティング (デルタ圧縮)
<vol_basename>.<kk>mer.ksx   — 配列メタデータ (配列長 + アクセッション)
```

ボリューム検出用の `.kvx` マニフェストファイルが常に生成されます:

```
<db_base>.<kk>mer.kvx        — ボリュームマニフェスト (テキスト形式、ボリュームベースネーム一覧)
```

`-max_freq_build` 使用時は全ボリューム共有の除外ビットセットファイルも生成されます (k 値ごとに 1 つ):

```
<db_base>.<kk>mer.khx        — 構築時除外ビットセット (全ボリューム共有)
```

例:
- 標準マルチボリューム (`nt`、ボリューム `nt.00`、`nt.01`): `nt.00.11mer.kix`、`nt.01.11mer.kpx`、`nt.11mer.kvx`、`nt.11mer.khx`
- 集約 DB (`combined`、ボリューム `foo`、`bar`): `foo.11mer.kix`、`bar.11mer.kix`、`combined.11mer.kvx`

`.khx` ファイルは 32 バイトヘッダ (マジック "KMHX"、フォーマットバージョン、k) に続き、`ceil(4^k / 8)` バイトのビットセットで構成されます。ビット *i* = 1 は k-mer *i* がボリューム横断の合算カウントに基づきインデックス構築時に除外されたことを示します。

ID ポスティングと位置ポスティングは別ファイルに格納されるため、Stage 1 フィルタリングが `.kpx` にアクセスすることはなく、ページキャッシュ効率が最大化されます。

## ソースからのビルド

### 依存関係

- C++17 コンパイラ (GCC >= 10, Clang >= 12)
- CMake >= 3.16
- NCBI C++ Toolkit (BLAST DB アクセス用)
- Intel TBB (並列化用)
- Parasail >= 2.6 (Stage 3 ペアワイズアライメント用)
- htslib >= 1.17 (SAM/BAM 出力用)
- Drogon (ikafssnhttpd 用、オプション)
- libcurl (HTTP クライアントモードおよびリモート取得用、オプション)

### 依存パッケージのインストール

NCBI C++ Toolkit 以外の依存パッケージを以下のコマンドでインストールできます。

**Ubuntu Server 24.04:**

```bash
sudo apt install build-essential cmake libtbb-dev liblmdb-dev libsqlite3-dev \
    libcurl4-openssl-dev libjsoncpp-dev
sudo apt install zlib1g-dev libbz2-dev liblzma-dev libdeflate-dev autoconf
sudo apt install libdrogon-dev uuid-dev libmariadb-dev libyaml-cpp-dev libbrotli-dev libhiredis-dev
```

2 行目は Parasail および htslib のソースビルドに必要な依存パッケージです。3 行目は Drogon および Ubuntu で `libdrogon-dev` が自動的に導入しない追加依存パッケージです。ikafssnhttpd が不要な場合は 3 行目を省略し、ビルド時に `-DBUILD_HTTPD=OFF` を指定してください。

**AlmaLinux 9 / Rocky Linux 9:**

```bash
sudo dnf config-manager --set-enabled crb
sudo dnf install -y epel-release
sudo dnf group install -y "Development Tools"
sudo dnf install -y cmake gcc-c++ tbb-devel lmdb-devel sqlite-devel \
    libcurl-devel jsoncpp-devel
sudo dnf install -y zlib-devel bzip2-devel xz-devel libdeflate-devel autoconf
sudo dnf install -y libuuid-devel openssl-devel
```

**Oracle Linux 9:**

```bash
sudo dnf config-manager --set-enabled crb
sudo dnf install -y oracle-epel-release-el9
sudo dnf group install -y "Development Tools"
sudo dnf install -y cmake gcc-c++ tbb-devel lmdb-devel sqlite-devel \
    libcurl-devel jsoncpp-devel
sudo dnf install -y zlib-devel bzip2-devel xz-devel libdeflate-devel autoconf
sudo dnf install -y libuuid-devel openssl-devel
```

**AlmaLinux 10 / Rocky Linux 10:**

```bash
sudo dnf config-manager --set-enabled crb
sudo dnf group install -y "Development Tools"
sudo dnf install -y cmake gcc-c++ tbb-devel lmdb-devel sqlite-devel \
    libcurl-devel jsoncpp-devel
sudo dnf install -y zlib-devel bzip2-devel xz-devel libdeflate-devel autoconf
sudo dnf install -y libuuid-devel openssl-devel
```

**Oracle Linux 10:**

```bash
sudo dnf config-manager --set-enabled crb
sudo dnf group install -y "Development Tools"
sudo dnf install -y cmake gcc-c++ tbb-devel lmdb-devel sqlite-devel \
    libcurl-devel jsoncpp-devel
sudo dnf install -y zlib-devel bzip2-devel xz-devel libdeflate-devel autoconf
sudo dnf install -y libuuid-devel openssl-devel
```

EL9 では `jsoncpp-devel` に EPEL、`lmdb-devel` に CRB リポジトリが必要です。EL10 ではいずれも CRB に収録されているため EPEL は不要です。各ブロックの最後から 2 行目は Parasail および htslib のソースビルドに必要な依存パッケージです。最終行は Drogon のソースビルドに必要な依存パッケージです。ikafssnhttpd が不要な場合は最終行を省略し、`-DBUILD_HTTPD=OFF` でビルドしてください。

### Parasail

ikafssn は Stage 3 ペアワイズアライメントに Parasail ライブラリを使用します。デフォルトではソースルート直下の `./parasail` を参照します。別の場所にインストール済みの場合は `-DPARASAIL_DIR` で指定してください。

Parasail のダウンロード・ビルド・インストールは、ikafssn ソースルートで以下を実行します:

```bash
curl -L -o parasail-2.6.2.tar.gz \
    https://github.com/jeffdaily/parasail/archive/refs/tags/v2.6.2.tar.gz
tar xf parasail-2.6.2.tar.gz
cd parasail-2.6.2
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$(realpath ../..)/parasail" \
    -DBUILD_SHARED_LIBS=OFF
make -j$(nproc)
make install
cd ../..
```

### htslib

ikafssn は SAM/BAM 出力に htslib を使用します。デフォルトではソースルート直下の `./htslib` を参照します。別の場所にインストール済みの場合は `-DHTSLIB_DIR` で指定してください。

htslib のダウンロード・ビルド・インストールは、ikafssn ソースルートで以下を実行します:

```bash
curl -L -o htslib-1.23.tar.bz2 \
    https://github.com/samtools/htslib/releases/download/1.23/htslib-1.23.tar.bz2
tar xf htslib-1.23.tar.bz2
cd htslib-1.23
autoreconf -i
./configure --prefix="$(realpath ..)/htslib" --disable-libcurl --disable-gcs --disable-s3
make -j$(nproc)
make install
cd ..
```

### NCBI C++ Toolkit

ikafssn のビルドには NCBI C++ Toolkit が必要です。デフォルトではソースルート直下の `./ncbi-cxx-toolkit` を参照します。別の場所にインストール済みの場合は `-DNCBI_TOOLKIT_DIR` で指定してください。

Toolkit 内のビルドサブディレクトリ名 (例: `CMake-GCC1330-Release`) はデフォルトで自動認識されますが、必要に応じて `-DNCBI_TOOLKIT_BUILD_TAG` で変更できます。

Toolkit のダウンロード・ビルド・インストールは、ikafssn ソースルートで以下を実行します:

```bash
curl -L -o ncbi-cxx-toolkit-public-release-30.0.0.tar.gz \
    https://github.com/ncbi/ncbi-cxx-toolkit-public/archive/refs/tags/release/30.0.0.tar.gz
tar xf ncbi-cxx-toolkit-public-release-30.0.0.tar.gz
cd ncbi-cxx-toolkit-public-release-30.0.0
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
| `BUILD_HTTPD` | ON | ikafssnhttpd をビルド (Drogon が必要) |
| `BUILD_CLIENT` | ON | ikafssnclient をビルド (HTTP モードで libcurl が必要) |
| `ENABLE_REMOTE_RETRIEVE` | ON | ikafssnretrieve で NCBI efetch を有効化 |
