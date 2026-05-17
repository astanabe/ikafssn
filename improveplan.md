# ikafssnsearch mode 1 性能改善プラン

本ドキュメントは `ikafssnsearch -mode 1` を NCBI nt 規模 (約 6.4 TiB の `.kix`) で実行した際の性能調査結果と、それに基づく改善案をまとめたものである。別のマシン上で本プランに沿って実装を進められるよう、現状コードへの参照と変更方針を具体的に記述する。

---

## 1. 対象実行と環境

### 1.1 対象コマンド

インデックス構築:
```
ikafssnindex -mode 1 -db /work/shimotsuki/ncbi_nt_20260107/nt \
  -o /work/shimotsuki/ncbi_nt_20260107 \
  -k 11 -t 21 -template_type both \
  -min_seq_length 64 -min_length_split 50000 -overlap_length 500 \
  -nthread 32
```

検索 (調査対象):
```
ikafssnsearch -mode 1 -ix /work/shimotsuki/ncbi_nt_20260107/nt \
  -query /home/shimotsuki/ikafssntest/mifish500.fasta \
  -o /home/shimotsuki/ikafssntest/nt_mifish500_ikafssn_k11t21_mode1.tsv \
  -output_format tsv -k 11 -t 21 -template_type both \
  -min_query_length 64 -stage1_min_score 0.1 \
  -strand 1 -nthread 32 -memory_limit 170G -v
```

### 1.2 計測値 (2026-05-16, PID 124514, 実行開始から約 7:44 時点)

| 項目 | 値 |
|---|---|
| ボリューム数 | 287 |
| 1 ボリュームの `.kix` サイズ | cod ≈ 16 GiB, opt ≈ 16 GiB (合計 ≈ 32 GiB/vol) |
| 全 `.kix` 合計 | 約 6.36 TiB (cod 3.18 + opt 3.18) |
| ステージ1 バッチ計画 (現状) | `posting_budget = 170 GiB / 32 GiB ≈ 5 vol/バッチ × ~57 バッチ` (ボリューム全体の `.kix` サイズを単位にバッチ数を決めている) |
| ステージ1 バッチ計画 (C1 後の見積) | 実需要 ≈ 200 MB/vol × 287 vol ≈ 60 GiB → 170 GiB の予算に **1 バッチで全ボリューム収容可能** (詳細 §3 C1-b) |
| クエリ | 500 配列, 平均 171 bp, 合計 85,745 bp |
| プロセス CPU 使用率 | 2,720 % / 3,200 % (≈ 85 % CPU bound) |
| スレッド状態 | 29 R + 3 `folio_wait_bit_common` |
| RSS / Cached / Swap | 3.7 GiB / 124 GiB / SwapPss 5 GiB |
| 累積 `read_bytes` | 4.68 TB (≒ 10 GB/s 持続), `syscr=144` (mmap 経由) |
| 進捗 | vol 159–165 を処理中 (≈ 58 %) |
| 推定総実行時間 | 約 13–14 分 |

### 1.3 環境メモ

- 物理メモリは 197 GB。`-memory_limit 170G` は申告値が大きすぎ、実効 file cache 上限 (約 125 GB) を超えた `MADV_WILLNEED` 分は次バッチ前にページ回収されている。
- ディスクは ~10 GB/s 持続できる NVMe。CPU 律速。
- 32 スレッド構成、AVX-512BW 以上対応 CPU。

### 1.4 実行経路 (本ベンチで通る code path)

- `-template_type both` 指定のため Stage 1 の active path は `stage1_one_strand_both` (`src/search/parallel_search.cpp:160-243`)。cod / opt の query stream を q_pos 整列でマージし、各 q_pos ブロックで `stage1_filter_accumulate` を 1〜2 回呼ぶ。1 ext_job あたり kix を **2 本** (cod + opt) 読みに行くため、B2 のヒープ確保コストは単一テンプレ経路の最大 2 倍。
- `-mode 1` 指定のため `search_orchestrator.cpp:361-377` の fast path に入り、**Stage 2A / Stage 2B はスキップ** される。よって本ベンチで効くのは Stage 1 関連 (B1〜B5, B7, C1〜C4, C6〜C8) と mode 1 fold (B6, C5) のみ。Stage 2 系の改善は別ベンチで評価する。
- ext_jobs 数 = 500 query × 287 vol × 1 strand = **143,500**。`states` ベクタも同サイズ。

---

## 2. ボトルネック分析

### B1. `apply_madvise_full(true)` による posting body 全面 WILLNEED — 読出量 80× 過剰

**該当**:
- `src/search/search_orchestrator.cpp:148-154` (`apply_stage1_madvise`)
- `src/index/kix_reader.cpp:91-105` (`KixReader::apply_madvise_full`) — `willneed=true` 時に `dict_size` 以降を `MADV_WILLNEED`
- `src/index/kix_reader.cpp:107-116` (`KixReader::apply_madvise_posting_random`) — 比較対象、posting body は `MADV_RANDOM`

Stage 1 の Tier4 経路では `MADV_WILLNEED` をボリュームの `.kix` 全体 (cod+opt で 32 GiB/バッチ) に対して発行している。一方 mode 1 でステージ1が実際に触る k-mer 数は

```
500 query × 約 322 k-mer (fwd both-template) ≈ 161,000 probe/volume
```

で、重複除去後でも数万。dictionary (約 32 MB) + (probe 数 × 平均 posting list ≈ 4 KB) ≈ **200 MB/volume** しか実需要がない。

**結果**: 16 GB の posting file に対し 200 MB しか触らないため、**約 80× の read amplification**。実測 4.68 TB はこの過剰先読みでほぼ説明できる (本来必要なのは 287 vol × 200 MB ≈ 60 GB)。

CPU 律速のためレイテンシ上は隠蔽されているが、I/O 帯域が遅い環境や `-memory_limit` が小さい環境ではこの項が支配項に変わる。

**Stage 2A への波及 (mode 2/3 への影響)**: `apply_stage2a_madvise` (`src/search/search_orchestrator.cpp:163-182`) は Tier4 で `kix.apply_madvise_full(true)` と `kpx.apply_madvise_full(true)` を、Tier3 で `kix.apply_madvise_full(true)`、Tier2 で `kpx.apply_madvise_full(true)` を呼ぶ。`KpxReader::apply_madvise_full(true)` も同形の posting body 全面 WILLNEED を発行する。mode 1 では実行されないが、mode 2/3 では同じ過剰先読みが kpx 側でも発生し、しかも posting list サイズが kpx > kix のため影響が大きい。本プランの C1 は実装時に Stage 2A 側にも適用すること。

### B2. `open_stream_kix` がプローブごとにヒープ確保 + 完全マテリアライズ

**該当**: `src/index/pfd_codec_tier.cpp:481-518` (`open_stream_kix`)
**呼出元**: `src/search/seq_id_decoder.hpp:31-37` (`SeqIdDecoder::ensure_decoded`)
**さらに呼出元**: `src/search/stage1_filter.cpp:115` (`SeqIdDecoder decoder(...)` が qi ループ内のスタックローカル) — `stage1_filter_accumulate_impl` の `for (size_t qi = 0; qi < n; qi++)` 内で k-mer 1 個につき 1 個の `SeqIdDecoder` が `StreamCtx` ごと作られて捨てられる
**both-mode 倍率**: `stage1_one_strand_both` (`src/search/parallel_search.cpp:160-243`) は cod / opt 双方で `stage1_filter_accumulate` を呼ぶため、1 ext_job あたりの確保回数は単一テンプレ経路の **2 倍**

現状コードは k-mer プローブ 1 回ごとに:
```cpp
std::vector<std::uint32_t> codec_in(body_words);  // ★ ヒープ確保 1
std::memcpy(codec_in.data(), posting_list + 4, body_bytes_avail);
ctx.decoded.resize(count);                         // ★ ヒープ確保 2 (再利用時は capacity 内で済む)
kix_codec().decodeArray(codec_in.data(), body_words,
                        ctx.decoded.data(), nvalue);
for (std::uint32_t i = 1; i < count; i++) {
    ctx.decoded[i] += ctx.decoded[i - 1];          // スカラ累積和
}
```
を実行する。

加えて `SeqIdDecoder` は qi ループ内のスタックローカルなので、`StreamCtx` (内部に `std::vector<uint32_t> decoded`) が **k-mer プローブごとに新規生成/破棄** される。生成のたびに `decoded` の capacity が 0 から始まるため、`ctx.decoded.resize(count)` は (再利用 thread_local 化されていれば不要な) 確実なヒープ確保になる。

なお `kix_codec()` 自体は `thread_local` (`src/index/pfd_codec_tier.cpp:63-66`) なので codec 内部 scratch は再利用されている。残るヒープ確保ホットスポットは `codec_in` (line 500) と `ctx.decoded.resize` (line 503) の 2 箇所のみ。

**規模**: 約 800,000 プローブ/バッチ × 57 バッチ × 2 (cod+opt) = 約 9,000 万回の動的確保。mifish のような conserved primer は **長い posting list (数十万〜数百万エントリ)** を引くため、1 回の `decoded.resize` で 1–10 MB のヒープ確保 + cumulative sum が走る。アロケータ pressure と L2/L3 ミスの主因。

### B3. AVX-512 scatter が num_seqs 規模のバッファを叩く

**該当**: `src/search/stage1_filter_simd.cpp:43-125` (`flush_avx512_t32`)

`_mm512_mask_i32gather_epi32` (lines 63-64, 70-71) / `_mm512_mask_i32scatter_epi32` (lines 82-83) が `scores[sid]` / `last_pos[sid]` (各 thread あたり `sizeof(ScoreT/PosT) × num_seqs` バイト = T32 では `4 × num_seqs`) にランダム書込みする。nt の 1 volume は `.ksx` が 100–200 MB あり fragment 数は数千万、Stage1Buffer は **SoA レイアウト** (score と last_pos が別の `AlignedBuffer`、`src/search/stage1_filter.hpp:29-47`) なので 1 thread あたり 2 本の数十 MB バッファ。

**結果**: L3 (一般的に 32–96 MB) に乗り切らないため、posting list が長い k-mer ほど scatter は全レーン L3 ミスとなり、`vpscatterdd` (microcoded) の実効スループットが大きく落ちる。CPU 時間を食う主要因の一つ。`tier == T32` 経路でのみ AVX-512 が選ばれ、T8/T16 は `flush_scalar_impl` にフォールバックする (`src/search/stage1_filter_simd.cpp:140-156`)。本ベンチでは num_seqs (fragment 数) が桁違いに大きいため T32 が選ばれる (`src/search/stage1_filter.hpp:67-72`)。

### B4. 同一 k-mer を query 数だけ重複デコード

**該当**:
- `src/search/parallel_search.cpp:434-463` (`run_stage1_jobs`) — `tbb::parallel_for(blocked_range(0, batch_indices.size()))` で `(query × volume × strand)` ext_job 軸を回す
- `src/search/parallel_search.cpp:297-360` (`dispatch_stage1`) — single / both をディスパッチ
- `src/search/stage1_filter.cpp:108-136` — qi ループ本体 (single 経路の `stage1_filter_accumulate_impl`)
- `src/search/query_preprocessor.hpp:15-29` — `QueryKmerData` は **query-major** な SoA (`fwd_positions[]`, `fwd_kmer_values[]`)、kmer-major へのグルーピングは持っていない

同一バッチ内の 500 query は同じ k-mer を多数共有する (例: mifish 12S プライマーは conserved region)。それでも `dispatch_stage1` は query (= ext_job 1 個) ごとに `posting_list_range → SeqIdDecoder → flush_batch_simd` を回すため、**1 つの common k-mer の posting list を 100 倍以上重複デコード**している。

`collect_position_hits` (`src/search/parallel_search.cpp:57-96`) も同じ構造で、kpx の `for_each_candidate` を query ごとに走らせる (mode 1 では非経路だが mode 2/3 で同じ問題)。なお `collect_position_hits` には既に `tls_pos_scratch()` (`src/search/parallel_search.cpp:29-32`) というスレッドローカル scratch があるが、kix の `SeqIdDecoder` 側に同等の仕組みは無い。

### B5. `dirty` リストの 2 度走査 + ランダム書き戻し

**該当**: `src/search/stage1_filter.cpp:58-70` (`Stage1Buffer::clear_dirty_typed`), `src/search/stage1_filter.cpp:141-167` (`stage1_filter_finish_impl`)

`stage1_filter_finish_impl` は `dirty` をなめて閾値以上の sid を candidate に積み (line 146-150)、その直後 `clear_dirty_typed<Tier>()` (line 152) が **再度** `dirty` をなめて `scores[sid]=0`, `last_pos[sid]=sentinel` を書き戻す (line 65-68)。

なお `dirty` に push されるのは `flush_batch_simd` で `scores[sid] == 0` だった瞬間 (= sid の初回観測) のみ (`src/search/stage1_filter_simd.cpp:72-93, 104, 119`)、その直後に必ず `scores[sid]++` が起きるため、`stage1_filter_finish_impl` に到達した時点で **`scores[sid] >= 1` がほぼ保証** される (sentinel が q_pos に一致する病的ケース以外)。したがって閾値以上の判定は `dirty` を走査する以外に手段が無い。

低 threshold (`-stage1_min_score 0.1` 等) では query 1 個あたりの `dirty.size()` が数百万に到達する場合があり、**ランダムインデックスへの書き戻し**が ext_job (= バッチ内 2,500 ジョブ程度) ごとに発生する。`dirty.size() > capacity/8` 程度になれば `memset(scores, 0, ...)` + `last_pos` sentinel の sequential reset の方が安い。

### B6. mode 1 結果集約がシリアル + バッチ末で解放されない

**該当**:
- `src/search/search_orchestrator.cpp:361-377` (mode 1 fast path)
- `src/search/search_orchestrator.cpp:480-486` (Stage 2 バッチ末解放) — mode 1 fast path は `return results;` で line 376 で抜けるため、**この解放ループは mode 1 では一切実行されない**
- `src/search/parallel_search.cpp:127-131` (`stage1_one_strand_single`), `:228-232` (`stage1_one_strand_both`) — mode 1 では `state.mode1_only = true; state.mode1_results = stage1_only_results(...)` が Stage 1 完了時にセットされる

```cpp
if (in.config.mode == 1) {
    std::vector<OrchestratorHit> results;
    for (size_t ji = 0; ji < ext_jobs.size(); ++ji) {
        ...
        for (auto& cr : st.mode1_results) {
            OrchestratorHit oh;
            oh.query_idx = ej.qi;
            oh.volume_idx = ej.vi;
            oh.volume_index = volume_indices[ej.vi];
            oh.cr = std::move(cr);
            results.push_back(std::move(oh));   // シングルスレッド fold
        }
        st.mode1_results.clear();   // ← size=0 にするだけ、capacity は残る
    }
    return results;
}
```

すべてのバッチ終了後、143,500 ext_job (= 500 query × 287 vol × 1 strand) 分の `mode1_results` を 1 スレッドで連結している。低 threshold では数 GB の `OrchestratorHit` を vector 拡張しながら直列に積むことになる。

さらに `JobState::mode1_results` は**Stage 1 バッチ末に解放されない**。`std::vector::clear()` は capacity を保持するため、Stage 1 が終わって fast path に入る時点で全 ext_job の `mode1_results` が累積保持されており、peak 時の anonymous heap が数 GB に達する。観測された 5 GB の SwapPss と整合する。

### B7. dict 領域の `MADV_HUGEPAGE` が file-backed mmap では効かない可能性

**該当**: `src/index/kix_reader.cpp:96-98` (`apply_madvise_full`), `:111-113` (`apply_madvise_posting_random`)

Linux の Transparent Huge Page は **file-backed mmap 領域** に対しては `MADV_HUGEPAGE` だけでは適用されないことが多い (CONFIG/khugepaged のグローバル設定依存)。32 MB/kix の dictionary 領域は実態として 4 KB ページのまま walk されており、`access_pair` 1 回ごとに TLB ミス可能性がある。実害は B2/B3 と比べれば小さいが、副次的に効く。

---

## 3. 改善案

### 3.0 前提となる設計決定 (確定済み)

詳細計画立案前に確定した方針。改善案 C1–C8 はこの前提に従って実装する。

- **C4 スコープ**: per-thread `Stage1Buffer` を維持したまま、外ループを unique kmer に転置 (内ループで `(qi, q_pos)` を回す)。per-query Stage1Buffer 案 (500 query × num_seqs × 4 B で数百 GB) は採用しない。per-kmer bitset 等のさらに深い再設計も本タスクのスコープ外。
- **C1 `needed_bytes[vi]` 計算タイミング**: バッチ planning 時に全 287 ボリュームを順に open して `posting_list_range` を全 unique k-mer に対して呼び、事前計算してからバッチを組む。スタートアップに数十秒〜数分の追加コストが乗るが、batch 計画を確定できるトレードオフを取る。
- **C1-d (KPX 側)**: C1 本体と同時実装。コード対称性を保ち、mode 2/3 ベンチの実施時にすぐ効果が出るようにする。
- **完了判定**: 目標値は事前に設定しない。各改善 (C1〜C8) を順次入れて §7 の計測値を残し、効果を実測で評価する。ある時点で「これ以上は効果が小さい」と判断したらそこで止める。

---

### C1. WILLNEED 戦略の見直し — 実需要バイト単位の prefetch + バッチ計画

**対象**:
- `src/search/search_orchestrator.cpp:59-97` (`plan_stage1_batches` — バッチ計画ロジック)
- `src/search/search_orchestrator.cpp:99-146` (`plan_stage2_batches` — Stage 2 用バッチ計画)
- `src/search/search_orchestrator.cpp:148-154` (`apply_stage1_madvise`)
- `src/search/search_orchestrator.cpp:163-199` (`apply_stage2a_madvise` / `release_stage2a_madvise`) — mode 2/3 への波及対応
- `src/index/kix_reader.cpp:91-116` (`apply_madvise_full` / `apply_madvise_posting_random`)
- `src/index/kix_reader.hpp:65-72` (`posting_list_range`) — レンジ取得 API
- `src/index/kpx_reader.{hpp,cpp}` (`KpxReader::apply_madvise_full`, posting range API) — Stage 2A Tier4/Tier2 で同じ問題

**既存パターンの再利用**: `plan_stage1_batches` は既に「1 ボリュームのコストが `posting_budget` を超えたら単独バッチで `Stage1Tier4::kTier1` (= `apply_madvise_posting_random` 経路) に落とす」フォールバックを持っている (line 87-93)。`plan_stage2_batches` も同様 (line 127-142, Tier2/Tier3/Tier1 を kix/kpx の収まり方で選ぶ)。C1 はこのコスト計算を「ボリューム全体サイズ」から「実需要バイト」に置き換える改造である。フォールバック分岐の骨格は流用できる。

**前提 API の確認**: `src/io/mmap_file.hpp:29` に `bool advise(size_t offset, size_t length, int advice);` が既にあるため、レンジ単位の `madvise` 発行は追加実装なしで可能。

**基本方針**: ボリューム全体の `.kix`/`.kpx` を `MADV_WILLNEED` で前ロードするのを完全にやめ、各バッチが**実際に触る k-mer の posting list バイト合計**を予算単位として、ボリュームを 1 個ずつ詰め込んでバッチを構成する。1 ボリューム単独でも予算を超える場合は、そのボリュームだけ `MADV_RANDOM` 経路にフォールバックして処理する。

#### C1-a: バッチ planner を「実需要バイト」ベースに作り直す

現状の `plan_stage1_batches` / `plan_stage2_batches` は「ボリュームの `.kix`/`.kpx` 合計サイズ」を予算単位にしている。これを「そのボリュームでバッチ全 query が触る posting list バイトの合計 + dict サイズ」に置き換える。

**事前計算 (3.0 の決定により確定)**: 全 query から `kmers_unique` (ユニーク k-mer 集合) を 1 度作り (ボリューム非依存)、planning フェーズで全 287 ボリュームを順に open → `posting_list_range` を全 unique k-mer に対して呼んで `needed_bytes[vi]` を確定する。dict pages は probe で fault されるが、その後 batch 実行時に再利用されるため無駄にはならない。スタートアップに追加コストが乗るがバッチ計画を確定できる。

**ボリューム単位の needed_bytes**:
```
open kix[vi]                              // dict pages はこの後 fault される
needed_bytes[vi]
   = dict_size(kix[vi])
   + sum( kix[vi].posting_list_range(k).len for k in kmers_unique )
   + (both-mode なら同様に kix_opt[vi] 側も加算)
   + (mode 2/3 なら kpx[vi] / kpx_opt[vi] 側も加算 — Stage 2 batch 用)
```

`posting_list_range` は dict の `access_pair` を 1 度叩くだけなので、k-mer 集合サイズ (本ベンチで ~数万) × 287 vol の dict probe は計算コスト的に問題ない。

**バッチ greedy 詰め込み**:
```
current_batch = []
current_bytes = 0
for vi in 0..num_volumes:
    if needed_bytes[vi] > posting_budget:
        // 単独でも溢れる: C1-c のフォールバック経路へ
        if current_batch: flush_batch(current_batch, tier=RangePrefetch)
        emit_batch([vi], tier=PostingRandom)
        current_batch = []; current_bytes = 0
    elif current_bytes + needed_bytes[vi] > posting_budget:
        flush_batch(current_batch, tier=RangePrefetch)
        current_batch = [vi]; current_bytes = needed_bytes[vi]
    else:
        current_batch.append(vi)
        current_bytes += needed_bytes[vi]
if current_batch: flush_batch(current_batch, tier=RangePrefetch)
```

本ベンチでは `needed_bytes[vi] ≈ 200 MB` (両テンプレ合計) なので、170 GiB の予算に **287 vol が 1 バッチで収まる** 見込み。バッチ間 open/close オーバーヘッドも 1 回に縮む。

**注意**: `posting_budget` は申告値であって OS の page cache 上限ではない。物理 RAM (197 GB) と実効 file cache 上限 (≈125 GB) を考慮して `-memory_limit` を実態より小さめに設定する運用は引き続き必要。

#### C1-b: バッチ実行時の範囲特化 prefetch

C1-a で組まれたバッチを実行する直前に、そのバッチが触る posting list のレンジだけを `MADV_WILLNEED` で先読みする。`apply_stage1_madvise` / `apply_stage2a_madvise` は `apply_madvise_full(true)` を呼ぶのをやめ、以下に置き換える:

```cpp
// pseudo (run_stage1_jobs の前段に挿入、または apply_stage1_madvise を作り直す)
// dict 部分: 全 query が dict probe するので一括 WILLNEED
kix[vi].mmap_.advise(0, dict_size, MADV_WILLNEED);
// posting body: バッチが触るレンジだけ
std::vector<std::pair<uint64_t,uint64_t>> ranges;
for (kmer in kmers_unique) {
    auto [off, len] = kix[vi].posting_list_range(kmer);
    if (len > 0) ranges.push_back({off, len});
}
sort(ranges);
coalesce_adjacent(ranges);  // 隣接 / 重複レンジを統合
for (auto [off, len] : ranges) kix[vi].mmap_.advise(off, len, MADV_WILLNEED);
// 残りの posting body は MADV_RANDOM (readahead 抑止)
```

C1-a の `needed_bytes` 計算で既に `posting_list_range` を呼んでいるので、計算結果を planner と共有してホットパスで二度叩かないようにする (planner → execution に `unique_ranges_per_volume` を渡す)。

#### C1-c: 単一ボリューム予算超過時のフォールバック (PostingRandom tier)

`needed_bytes[vi] > posting_budget` のボリュームは、そのボリューム単独でバッチを組み、`apply_madvise_posting_random()` で処理する (現状の Tier1 経路相当)。posting body へは `MADV_WILLNEED` を出さず `MADV_RANDOM` のみ、dict は WILLNEED + HUGEPAGE。

このフォールバックは「ベンチ環境でクエリが極端に多く、1 ボリュームの実需要が予算を超える」例外ケースのため。本ベンチではまず発生しない。

#### C1-d: KPX 側にも同戦略を適用 (mode 2/3 向け、C1 本体と同時実装)

`apply_stage2a_madvise` の Tier4/Tier2 で `kpx.apply_madvise_full(true)` が呼ばれているため、同じ posting body 全面 WILLNEED が kpx 側でも発生する。kpx の posting list サイズは kix より大きいため、本来の影響は kix より大きい。

Stage 2A バッチでは、Stage 1 で決まった `candidate_sids` と query k-mer 集合から触る kpx posting list の範囲が決まる。C1-a/C1-b と同形のレンジ収集 + 範囲特化 madvise を `KpxReader` 側にも入れる。mode 1 では実行されないため本ベンチでの効果は無いが、コード対称性と mode 2/3 ベンチでの効果のため、3.0 の決定により **C1-a/b/c と同時に実装する**。

**C1 全体の期待効果**: 実 I/O 量は 4.68 TB → 数十 GB オーダに縮む (理論下限は 287 vol × 200 MB ≈ 60 GB)。バッチ数が 57 → 1 程度に減ることで、`folio_wait_bit_common` で 3 スレッドが待たされる現象も消える。CPU 律速の現状では elapsed の短縮幅は限定的だが、低速ディスク環境や `-memory_limit` が物理 RAM より大きい環境では決定的に効く。

### C2. `SeqIdDecoder` / `StreamCtx` をスレッドローカル再利用に + memcpy 排除

**対象**:
- `src/index/pfd_codec_tier.cpp:481-518` (`open_stream_kix`)
- `src/search/seq_id_decoder.hpp:21-102` (`SeqIdDecoder` クラス)
- `src/search/stage1_filter.cpp:115` (qi ループ内 `SeqIdDecoder decoder(...)`)
- `src/search/parallel_search.cpp:76-77` (`collect_position_hits` 内の同様パターン — mode 2/3 影響)

#### C2-a: `StreamCtx` を thread_local 化
`tls_pos_scratch()` (`src/search/parallel_search.cpp:29-32`) と同じ思想で `pfd::StreamCtx` (内部の `decoded` バッファ) を thread_local 化する。または `SeqIdDecoder` を qi ループの外側で生成し、各イテレーションで `reset(data, end)` する API を追加する。

具体案:
- `SeqIdDecoder::reset(const uint8_t* data, const uint8_t* end)` を追加 (`decoded_` フラグだけ false に戻す)
- ループ外で `static thread_local SeqIdDecoder dec;` を持ち、各 qi で `dec.reset(...)` を呼ぶ
- `StreamCtx::decoded` の capacity は ensure_capacity 的に成長維持
- 既に `kix_codec()` (`src/index/pfd_codec_tier.cpp:63-66`) は thread_local なので、本変更で「open_stream_kix の毎呼出 2 ヒープ確保」が完全に消える

#### C2-b: `codec_in` の memcpy を排除
EF 辞書の構造上、posting list は word-aligned で配置可能なはず。`posting_list + 4` が 4 byte aligned かを実行時にアサート / 静的保証した上で、`(const uint32_t*)(posting_list + 4)` を `decodeArray` に直接渡す。アラインメントが保証できないケースだけ scratch buffer (thread_local `std::vector<uint32_t>`) にコピーするフォールバックを残す。

注意: `FastPForLib::SIMDFastPFor<4>::decodeArray` は入力ポインタを `__m128i` でロードするためアラインメント要件は実装依存。実装時に FastPFor 側の `loadu` 使用を確認すること。

**期待効果**: ホットパスのヒープ確保が毎呼出 2 回 → 0 回。CPU 数十 % 削減見込み。both-mode では cod / opt 双方で発火するため、効果は単一テンプレ経路の 2 倍。

### C3. PFD 結果をブロック単位で stream decode

**対象**: `src/index/pfd_codec_tier.cpp:481-518` + 上位の `SeqIdDecoder` API + `src/index/pfd_codec.hpp:108-123` (`StreamCtx` / `open_stream_kix` 宣言)

mode 1 のステージ1 は「各 sid に対し `last_pos != q_pos` なら score++」を行うだけで、posting list 全体をメモリ上に持つ必要がない。

方針:
- FastPFor の `SIMDFastPFor<4>` は 128 要素単位 (`src/index/pfd_codec.hpp:68` の `kPfdBlockSize = 128`、`src/index/pfd_codec_tier.cpp:56` の `kBlockSize = 128`)。posting list を 128 要素ブロックずつ decode する API (`decode_block_kix(...)` 等) を `pfd_codec_tier.cpp` に追加
- `SeqIdDecoder` を「現在ブロックの 128 要素を持つ + 必要に応じて次ブロックを decode」のストリーミング型に作り替え
- 累積和 (`decoded[i] += decoded[i-1]`) を SIMD prefix-sum kernel 化 (AVX2/AVX-512 で 4–8× 高速化可能)
- 注意: 現在の `open_stream_kix` は「[abs_first, d1, d2, ...]」を 1 回でデコードしてから累積和を取る。ブロック化する場合、cross-block の prev 値伝搬を `SeqIdDecoder` 側で保持して、次ブロック先頭の累積和起点に使う必要がある

**期待効果**: working set が 128 × 4 B = 512 B となり L1 で完結。長い posting list でアロケーション量も削減。

### C4. k-mer 単位への転置 (query-major → kmer-major)

**対象**:
- `src/search/parallel_search.cpp:434-463` (`run_stage1_jobs`) のジョブ単位を再設計
- `src/search/parallel_search.cpp:160-243` (`stage1_one_strand_both`) — both-mode の q_pos マージ構造も影響
- `src/search/query_preprocessor.hpp:15-29` (`QueryKmerData`) — kmer→(qi, q_pos)\* の inverted view を追加するか、preprocess_query で同時生成するか要設計判断

**確定方針 (3.0 による)**: per-thread `Stage1Buffer` を維持したまま、k-mer 単位の外ループで posting list を 1 度だけ decode し、内ループで `(qi, q_pos)` ペアを回す hybrid 構成を採用する。per-query Stage1Buffer (500 query × num_seqs × 4 B = 数百 GB) は採用しない。

バッチに入る全 query から `(qi, q_pos, kmer)` を集めて kmer でグループ化:

```
バッチ前処理:
  triples = [(qi, q_pos, kmer) for qi in queries for (q_pos, kmer) in query]
  group_by_kmer(triples)  // ハッシュ or 基数ソート

ステージ1 ループ (per-thread Stage1Buffer + kmer-major 外ループ):
  parallel_for each kmer_idx in unique_kmers:
    open_stream_kix(kmer_idx) once → SeqIdDecoder (thread_local)
    while decoder.has_more():
      sid_batch = decoder.next_batch(...)
      for each (qi, q_pos) referencing kmer_idx:
        flush_batch_simd<Tier>(sid_batch, count, q_pos,
                               scores[qi], last_pos[qi], dirty[qi])
```

検討事項:
- **スコアバッファの軸**: per-thread Stage1Buffer を維持するが、論理的には (qi × sid) 二軸が必要になる。3.0 の確定方針では Stage1Buffer 1 個を qi で切り替え、kmer ループの外側 (またはバッチごと) で qi を回す擬似 hybrid にする。具体的な切り替え方は plan モードで詰める。
- **both-mode との統合**: `stage1_one_strand_both` の cod / opt 統合は kmer-major 化と相性が良い (q_pos マージループ自体が不要になる)。ただし `stage1_filter_accumulate` の per-position dedup (`last_pos[sid] != q_pos`) が cod / opt をまたいでも保たれるよう、kmer iteration 順を q_pos 昇順に保つ必要がある。
- mifish のような common k-mer 集中型クエリで効果大。

**期待効果**: 同一 k-mer の重複 decode を排除 (mifish のような conserved primer で 10–50×)。

### C5. mode 1 結果集約の並列化 + バッチ末解放

**対象**: `src/search/search_orchestrator.cpp:361-377` (mode 1 fast path), `:480-486` (Stage 2 バッチ末解放 — mode 1 では未実行)

#### C5-a: Stage 1 バッチ末 flush + 解放
Stage 1 の各バッチ末 (`src/search/search_orchestrator.cpp:339-348` のクローズ処理直後、または `run_stage1_jobs` 戻り直後) で、そのバッチの `batch_indices` に含まれる ext_job について:
1. `mode1_results` をスレッド分割可能な per-batch 結果領域へ flush
2. `std::vector<ChainResult>().swap(st.mode1_results)` で capacity ごと解放 (`.clear()` だけでは不足)

これで peak 占有が `posting_budget` 相当に収まり、SwapPss が消える。

#### C5-b: 結果集約の並列化
3 通りの選択肢を本プランで提示する (詳細は §4 参照)。実装としては **C5-b-iii の 2-pass prefix-sum 方式** を推奨する。C5-a を先に入れて per-batch fold にしてしまえば、各バッチの fold が小さくなるため C5-b の並列化要件も緩む。

### C6. `clear_dirty_typed` の閾値ベース fallback

**対象**: `src/search/stage1_filter.cpp:58-70` (`clear_dirty_typed`), `:25-33` (`reset_all_typed`), `:141-167` (`stage1_filter_finish_impl`)

```cpp
template <Stage1Tier Tier>
void Stage1Buffer::clear_dirty_typed() {
    if (dirty.size() * 8 > capacity) {
        // dirty 比率が高い → sequential reset の方が安い
        reset_all_typed<Tier>(*this);
    } else {
        // 現状の random write 経路
        ...
    }
    dirty.clear();
}
```

加えて `stage1_filter_finish_impl` の dirty 走査と reset を 1 パスに統合し、`scores[sid] >= threshold` を判定したその場で `scores[sid]=0`, `last_pos[sid]=sentinel` を書き戻す形にする (cache を 2 度叩かない)。`stage1_topn > 0` の場合は nth_element / sort を別途行うので、1 パス統合はあくまで `min_stage1_score` フィルタリングと reset の融合に限る。

T8/T16 では `scores` と `last_pos` が同じ cache line に収まるよう SoA → AoS 変更も検討余地あり (ただし SIMD scatter のパターンが変わるため要ベンチ)。本ベンチ (num_seqs 数千万) は T32 が選ばれるため、本案は T8/T16 ベンチ環境でのみ効く可能性がある。

### C7. dict 領域の huge page 戦略を明示化

**対象**: `src/index/kix_reader.cpp:91-116`

`MADV_HUGEPAGE` は file-backed mmap では効きにくいため、以下のいずれかに変える:

- 起動時に dict 領域を `MAP_POPULATE` で一括 faultin させ、`MADV_RANDOM` で readahead を抑えるだけにする。`MADV_HUGEPAGE` の行は削除して挙動の予測性を上げる。
- もし dictionary を anonymous mmap 上にコピーしてから使う構成が許容できるなら、そちらに transparent huge page が乗る (メモリは 32 MB × num_vol ≈ 9 GB 増)。

予測性重視で前者を推奨。

### C8. OidFilter 無し時の hot loop 特殊化

**対象**: `src/search/oid_filter.hpp:33-37` (`pass`), `src/search/stage1_filter.cpp:80-137` (qi ループ内 `filter.pass(sid)` line 122)

`stage1_filter_accumulate_impl` を `filter_mode` でテンプレート特殊化し、`OidFilterMode::kNone` 時は `filter.pass(sid)` 呼出を完全に消す。ステージ1 hot loop の per-sid 分岐が 1 つ減る。`OidFilter::mode()` (`src/search/oid_filter.hpp:39`) で実行時に判定できる。

実装メモ:
```cpp
template <typename KmerInt, Stage1Tier Tier, bool HasFilter>
static void stage1_filter_accumulate_impl(...);
```
として、呼出側で `filter.mode() == OidFilterMode::kNone` で分岐ディスパッチする。本ベンチでは `-seqidlist` 未指定のため `kNone` 経路に入る。

---

## 4. 結果集約と TBB tree reduction

### 4.1 `tbb::parallel_reduce` は確かに tree reduction

`tbb::parallel_reduce` は range を分割して葉で部分結果を作り、ペアでマージしながら上方向にたどる典型的な tree reduction を実装している。ご指摘の通り、本ケース (`std::vector<OrchestratorHit>` を集める) に適用可能。

ただし vector 連結は 1 回の join が `O(|a|+|b|)` のコピー/move を伴うため、tree 全体では **`O(N log P)`** の総作業量になる。

### 4.2 候補となる 3 つの実装パターン

#### (i) `tbb::combinable` + `combine_each` (最小変更)
既に `run_stage2b_jobs` で使われている方式 (`src/search/parallel_search.cpp:538-565`)。スレッドローカルに溜め、最後にメインスレッドで連結。
- 実装が最も簡単
- combine_each 自体はシリアル走査だが、各スレッドのローカル vector は独立に積まれるので並列化の効果はある
- 最終 fold は依然 `O(N)` シリアル

#### (ii) `tbb::parallel_reduce` (tree reduction)
range を `[0, ext_jobs.size())` で取り、各 leaf が `std::vector<OrchestratorHit>` を返し、join で連結。
- 真の tree reduction だが連結コストが `O(N log P)`
- vector の move 連結を join 内で実装する必要がある

```cpp
std::vector<OrchestratorHit> results = tbb::parallel_reduce(
    tbb::blocked_range<size_t>(0, ext_jobs.size()),
    std::vector<OrchestratorHit>{},
    [&](const tbb::blocked_range<size_t>& r, std::vector<OrchestratorHit> acc) {
        for (size_t ji = r.begin(); ji != r.end(); ++ji) {
            for (auto& cr : states[ji].mode1_results) {
                OrchestratorHit oh{ ej.qi, ej.vi, vol_idx[ej.vi], std::move(cr) };
                acc.push_back(std::move(oh));
            }
            std::vector<ChainResult>().swap(states[ji].mode1_results);
        }
        return acc;
    },
    [](std::vector<OrchestratorHit> a, std::vector<OrchestratorHit> b) {
        if (a.empty()) return b;
        a.insert(a.end(), std::make_move_iterator(b.begin()),
                          std::make_move_iterator(b.end()));
        return a;
    });
```

#### (iii) 2-pass + `tbb::parallel_scan` (推奨)
最速。総作業量 `O(N)`、並列度 `P`。

```
pass 1 (sizes):   sizes[ji] = states[ji].mode1_results.size();
pass 2 (scan):    offsets[ji] = prefix_sum(sizes)  // tbb::parallel_scan
pass 3 (write):   results.resize(total);
                  parallel_for(ji) {
                      for (size_t k = 0; k < sizes[ji]; ++k)
                          results[offsets[ji] + k] = ...;  // ロックなし
                  }
```

- 1 度しかメモリ確保しない (`results.resize(total)`)
- 書込は同期不要 (オフセット既知)
- vector 連結のオーバーヘッドなし

**本プランでは (iii) を推奨**。実装難度は (i) > (iii) > (ii) の順で (iii) が中間だが、性能差は (iii) > (ii) > (i) になる見込み。最初に (i) で動作確認 → (iii) で性能化する 2 段構えでも良い。

### 4.3 適用箇所
- `src/search/search_orchestrator.cpp:361-377` の mode 1 fast path
- `src/search/parallel_search.cpp:538-565` (`run_stage2b_jobs` 内の `tls_results.combine_each` line 561-565) — Stage 2/3 を使う場合の同様の改善対象。ただし stage2b はすでに per-batch fold (`out_results` を引数で受け取り `insert` する) になっており、mode 1 fast path より状況はマシ

---

## 5. 優先度と期待効果

| 優先度 | 案 | 期待効果 | 実装難度 |
|---|---|---|---|
| 高 | **C2** thread_local StreamCtx + memcpy 排除 | プローブごとの確保 0 化、CPU 数十 % 減 (both-mode で 2 倍効く) | 中 |
| 高 | **C1** バッチ planner を実需要バイト単位に + 範囲特化 prefetch | I/O 4.68 TB → 数十 GB (理論下限 ~60 GB)、バッチ数 57 → 1、`folio_wait_bit_common` 待ち解消、低速ディスクで決定的 | 中〜高 |
| 高 | **C4** k-mer 単位への転置 | conserved query で decode 10–50× 減 | 高 (再設計) |
| 中 | **C5** mode 1 fold 並列化 + Stage 1 バッチ末解放 (§4 の (iii) 推奨) | swap 流入解消、テール時間短縮 | 中 |
| 中 | **C3** ブロック単位 stream decode | L1 局在化、cumsum SIMD 化 | 中 |
| 中 | **C6** dirty reset の閾値切替 + 走査 1 パス化 | 低 threshold 時の reset O(n) → O(1) | 低 |
| 低 | **C7** dict huge page 戦略の予測性向上 | 副次的、TLB ミス減 | 低 |
| 低 | **C8** OidFilter 無し時の hot loop 特殊化 | 分岐 1 個削減 | 低 |

CPU 律速の現状を最も大きく改善するのは **C2 → C4 → C5** の順 (見込み)。I/O が支配項になる別環境では **C1** が最大効果。本ベンチは `-template_type both` のため C2 の効果は実質 2 倍と見込まれる (cod / opt 双方の qi ループで毎回 `SeqIdDecoder` を構築している)。

ただし表内の「期待効果」はあくまで分析からの推定であり、数値目標ではない (3.0 の決定参照)。各改善後に §7 の計測を行って実値を残す運用。

---

## 6. 実装順序の推奨

1. **C6** (最小変更で安全に効果確認)
2. **C8** (テンプレ特殊化のみ、リグレッション可能性低)
3. **C2** (`SeqIdDecoder` API 拡張 + thread_local 化, alignment 確認 → memcpy 排除) — 本ベンチで最も大きく効くと予想
4. **C5** (まず C5-a で Stage 1 バッチ末解放 → C5-b-iii で 2-pass 並列化)
5. **C1** (バッチ planner を実需要バイト単位に作り直し + 範囲特化 prefetch を一括実装)
   - C1-a (planner 改造) と C1-b (範囲特化 madvise) は密結合 (planner で収集した `unique_ranges_per_volume` を execution に渡す) のため同時実装
   - C1-c (単一ボリューム fallback) は既存 Tier1 分岐をそのまま流用、必要に応じて閾値計算だけ差し替え
   - C1-d (KPX 側) はコード対称性のため同時に投入
6. **C3** (FastPFor の block 単位 decode API 追加, prefix-sum SIMD 化)
7. **C4** (k-mer 単位への転置, ステージ1 再設計) — 最も影響範囲が広いので最後

各ステップ後に下記 §7 の検証手順で効果を測る。C2 と C5 は独立に効くため、どちらを先に入れても他方の効果測定を歪めない。C1 は I/O 量変化が主指標なので CPU 律速環境では elapsed が大きく動かない可能性がある — その場合は `/proc/$PID/io` の `read_bytes` 差分とバッチ数 (`logger->info("Stage 1 batch plan: ...")` `src/search/search_orchestrator.cpp:306-308`) の変化で判断する。

---

## 7. 検証手順

**前提 (3.0 の決定により)**: 数値的な完了目標は事前に設定しない。各改善 (C1〜C8) を順次入れて以下の計測を行い、効果の有無と大きさを実測で記録する。「これ以上は効果が小さい」と判断したらそこで止めるという運用にする。

### 7.1 マイクロベンチ
- `bench/` 配下に Stage 1 hot loop 単体ベンチを追加 (まだ無ければ)。1 ボリューム固定で代表的 posting list 長 (短/中/長) のサンプルを用意し、`open_stream_kix` 単体および `flush_batch_simd` 単体の per-call サイクル数を測る。
- `perf stat -e cycles,instructions,L1-dcache-load-misses,LLC-load-misses,dTLB-load-misses` で命令あたりサイクル / L3 ミス / TLB ミスを変化させて比較。

### 7.2 全体ベンチ (本書の対象実行で)

すべて「改善前との差分」を記録する目的の計測。数値目標は持たない。

- 実行前後で `/usr/bin/time -v` の Elapsed, Maximum resident set size, Major/Minor page faults を記録・比較。
- `/proc/$PID/io` の `read_bytes` を周期サンプル (例: `while sleep 5; do cat /proc/$PID/io; done`) し、I/O 量の絶対値と単位時間あたり量を記録。C1 の主指標。
- `vmstat 1` で `si/so` (swap I/O) と `bi/bo` (block I/O) を観察。C5-a の主指標は SwapPss の減少。
- `cat /proc/$PID/status | grep -E 'VmRSS|VmPeak'` で peak RSS を確認。
- `-v` ログの `Stage 1 batch plan: N batch(es) over M volume(s) (tier4=X, tier1=Y)` 行 (`src/search/search_orchestrator.cpp:306-308`) を見て、C1 後のバッチ数の変化と `tier1` (= 単一ボリューム fallback、C1-c の発火回数) を記録。
- スレッド状態を `for t in /proc/$PID/task/*; do cat $t/wchan; echo; done | sort | uniq -c` で集計し、`folio_wait_bit_common` で待っているスレッド数の変化を観察 (C1 の効果指標)。

### 7.3 結果同値性
- 改善前後で出力 TSV が **完全一致** することを `diff` で確認。順序が変わる場合は `sort -k1,1 -k2,2 -k3,3 | sha256sum` で同値性を担保。
- `test/` 配下の既存テストすべてが pass すること。特に:
  - `test/ssu_test_fixture.hpp` を使う mode 1 / mode 2 / mode 3 テスト
  - both-mode (`-template_type both`) を含むテスト — `stage1_one_strand_both` の挙動が C2/C4 で変わるため必須
  - `-seqidlist` を使うテスト — C8 の `OidFilter::pass` 経路を網羅
- C4 (kmer-major 化) は per-position dedup の不変条件 (`last_pos[sid] != q_pos` で +1) を壊しやすいので、both-mode テストの数値同値性を最重点で確認する。

### 7.4 ログメトリクス追加 (任意)
- 既存 `logger->info` (`src/search/search_orchestrator.cpp:306-308, :354-355`) に加え、バッチごとの: ユニーク probed k-mer 数 / 実 madvise バイト数 / decode 回数 / dirty.size() ピーク を `-v` 時に出力させると、各改善の効果が定量で見える。

---

## 8. やってはいけないこと / 制約

- **インデックスフォーマット (`.kix`/`.kpx`/`.ksx`/`.khx`) の変更は本プランの範囲外**。すべてリーダ側 (decoder / orchestrator) の改善で完結する。フォーマットを変更する必要が生じた場合は `format_version` bump を伴うため、別タスクとして切る。
- **canonical vocabulary を守る** (`CLAUDE.md` の Terminology 節)。新規追加コード / コメント / ログでも `dictionary`, `posting list`, `posting file`, `posting list header/body`, `parent`, `fragment`, `parent-relative coordinates` を使う。`offset table`, `posting data region`, `payload`, `document`, `OID-level sequence`, `subsequence/chunk/window/tile`, `OID-relative coordinates` などの禁止語を導入しない。既存箇所にこれらが残っていたら改善作業の一環で書き換える。
- **`IKAFSSN_VERSION` の更新を忘れない**: `src/`/`include/`/`test/` を編集してコミットする場合は `CMakeLists.txt` の `IKAFSSN_VERSION` をコミット日 (`0.1.YYYY.MM.DD`) に揃える (`CLAUDE.md` の Version Management 節)。
- **"Phase X" 表記をコードやドキュメントに残さない** (`CLAUDE.md` の Code Comments 節)。コミットメッセージや PR 本文では使って構わない。
- **コメントは現在の実装を説明するもののみ**。"以前は X だった / v10 → v11" のような変更履歴ナレーションをコードに残さない。
- **mode 2 / mode 3 経路を壊さない**: 本ドキュメントの主対象は mode 1 だが、`SeqIdDecoder` / `Stage1Buffer` / WILLNEED 計画 / 結果集約は mode 2/3 でも共有される。改善のたびに既存のステージ2A/2B/3 テストを必ず通す。
- **NUMA / アーキ依存をハードコードしない**: AVX-512 を前提とせず、`util/simd_dispatch.hpp` 経由で tier dispatch する既存方針を保つ。

---

## 9. 参考: 関連ファイル一覧

ステージ1 / mode 1 関連:
- `src/ikafssnsearch/main.cpp` — エントリポイント、preprocess、orchestrator 呼出、TSV 出力
- `src/search/search_orchestrator.{hpp,cpp}` — ステージ1/2 のバッチ計画とドライバ。`apply_stage1_madvise`/`apply_stage2a_madvise` (line 148-199)、mode 1 fast path (line 361-377)、Stage 2 バッチ末解放 (line 480-486)
- `src/search/parallel_search.{hpp,cpp}` — ステージ1/2A/2B のディスパッチと TBB 並列化。`tls_pos_scratch` (line 29-32)、`collect_position_hits` (line 57-96)、`stage1_one_strand_single` (line 98-138)、`stage1_one_strand_both` (line 160-243)、`dispatch_stage1` (line 297-360)、`run_stage1_jobs` (line 434-463)、`run_stage2b_jobs` (line 503-566)
- `src/search/stage1_filter.{hpp,cpp}` — `Stage1Buffer` (SoA レイアウト、tier 別 ScoreT/PosT)、`stage1_filter_accumulate_impl` (qi ループ line 80-137)、`stage1_filter_finish_impl` (line 141-167)、`clear_dirty_typed` (line 58-70)
- `src/search/stage1_filter_simd.{hpp,cpp}` — AVX-512 scatter kernel `flush_avx512_t32` (line 43-125)、ディスパッチ `flush_batch_simd` (line 130-156)
- `src/search/seq_id_decoder.hpp` — `.kix` posting list の sid イテレータ (line 21-102)
- `src/search/posting_decoder.hpp` — `.kpx` candidate-driven decoder (mode 2/3)
- `src/search/oid_filter.{hpp,cpp}` — accession 包含/除外フィルタ。`pass()` は line 33-37
- `src/search/query_preprocessor.{hpp,cpp}` — k-mer 抽出、閾値解決、`QueryKmerData` SoA (line 15-29)
- `src/search/result_writer_mode1.cpp` — mode 1 TSV/JSON 並列ライタ

インデックス I/O:
- `src/index/kix_reader.{hpp,cpp}` — `.kix` mmap, EF dictionary 経由の `posting_list_range` (hpp line 65-72)。`apply_madvise_full` (cpp line 91-105)、`apply_madvise_posting_random` (cpp line 107-116)
- `src/index/kpx_reader.{hpp,cpp}` — `.kpx` mmap。`apply_madvise_full` は同形パターン
- `src/index/ksx_reader.{hpp,cpp}` — fragment/parent メタ
- `src/index/pfd_codec.{hpp,cpp}` — PFD codec dispatcher。`StreamCtx` (line 112-116)、`open_stream_kix` 宣言 (line 123)、`kPfdBlockSize = 128` (line 68)、`PosDecodeScratch` (line 131-138)
- `src/index/pfd_codec_tier.cpp` — tier 別 FastPFor wrapper。`kix_codec()` thread_local (line 63-66)、`open_stream_kix` 実装 (line 481-518)、`kBlockSize = 128` (line 56)

下位:
- `src/io/mmap_file.{hpp,cpp}` — `madvise` ラッパ。`advise(offset, length, advice)` (hpp line 29)、`advise_dict_posting` (hpp line 35)
- `src/index/ef_codec.{cpp,hpp}` — Elias-Fano dictionary

ベンチ:
- `bench/` 配下 (Stage 1 単体ベンチを必要に応じて追加)

---

## 進捗状況

### Phase 1: C6 — dirty 走査の閾値切替 + finish 1 パス融合
- 日付: 2026-05-17
- 主な変更:
  - `src/search/stage1_filter.cpp` の `Stage1Buffer::clear_dirty_typed<Tier>` に `dirty.size() * 8 > capacity` の閾値分岐を追加し、dense ケースで `reset_all_typed<Tier>` に切り替え
  - `stage1_filter_finish_impl<Tier>` の `stage1_topn == 0` 経路で dirty 走査と reset を 1 パスに融合 (高 dirty 比率では bulk reset、低 dirty 比率では per-index に統合)
  - `bench/bench_stage1.cpp` に `BM_Stage1ClearDirtyT32` を追加 (num_seqs × dirty 比率 1/64–1/2 の格子)
  - `test/test_stage1.cpp` に `test_clear_dirty_bulk_reset` (bulk reset 経路) と `test_stage1_finish_no_side_effects` (融合 finish の冪等性) を追加
  - `CMakeLists.txt` の `IKAFSSN_VERSION` を `0.1.2026.05.17` へ bump
- テスト: `ctest --test-dir build` 全 158 件 pass
- 観察: SIMD 全変種 (`test_stage1_force_*`) も含めて回帰なし。マイクロベンチ結果と本番ベンチは別マシン実施待ち
- 残課題 / 次フェーズへの引き継ぎ: なし。Phase 2 (C8) へ進む

### Phase 2: C8 — OidFilter HasFilter NTTP specialization
- 日付: 2026-05-17
- 主な変更:
  - `src/search/stage1_filter.cpp` の `stage1_filter_accumulate_impl` に `bool HasFilter` NTTP を追加し、`HasFilter == false` では `filter.pass(sid)` 分岐を `if constexpr` で完全に除去
  - 新ヘルパ `dispatch_accumulate_by_filter<KmerInt, Tier>` で `filter.mode() == OidFilterMode::kNone` を見て 2 つの特殊化に振り分け
  - `stage1_filter_impl` / `stage1_filter_accumulate` は両方この新ディスパッチを経由
  - 公開シンボル (`stage1_filter<KmerInt>` / `stage1_filter_accumulate<KmerInt>`) のシグネチャ・external instantiation は不変
  - `test/test_stage1.cpp` に `test_stage1_filter_none_vs_pass_all` を追加 (kNone 経路と全 OID 通過 kInclude 経路の数値同値性)
- テスト: `ctest --test-dir build` 全 158 件 pass。`test_oid_filter` の include/exclude 経路も pass
- 観察: `filter.pass()` 呼出が hot loop から消えるため per-sid 分岐が 1 個減る。本番ベンチでの elapsed 変化は別マシン測定待ち
- 残課題 / 次フェーズへの引き継ぎ: なし。Phase 3 (C2) へ進む

### Phase 3: C2 — SeqIdDecoder thread_local + open_stream_kix heap alloc 排除
- 日付: 2026-05-17
- 主な変更:
  - `src/search/seq_id_decoder.hpp` に `reset(const uint8_t*, const uint8_t*)` API を追加 (StreamCtx::decoded の capacity を保持したまま状態リセット)。未使用の `explicit SeqIdDecoder(const uint8_t*)` ctor を削除
  - `src/index/pfd_codec_tier.cpp::open_stream_kix` の `std::vector<uint32_t> codec_in(body_words)` を排除。`posting_list + 4` が u32-aligned (常に成立) なら FastPFor へ直接 pass、不揃いなら `thread_local std::vector<uint32_t> codec_in_scratch` に fall back (どちらも per-call heap alloc なし)
  - `src/search/stage1_filter.cpp::stage1_filter_accumulate_impl` の qi ループ内で `SeqIdDecoder decoder(...)` を構築していたのを、関数先頭の `thread_local SeqIdDecoder decoder` を `decoder.reset(...)` で再利用する形に hoist
  - `src/search/parallel_search.cpp::collect_position_hits` (mode 2/3 共有経路) でも同形に hoist
  - `test/test_stage1.cpp` に `test_seq_id_decoder_reset_across_runs` を追加 (3 種類の query を 1 スレッド上で連続呼出 → 独立呼出結果と同値であることを確認)
- テスト: `ctest --test-dir build` 全 158 件 pass。全 SIMD 強制バリアントも pass
- 観察: hot loop での per-probe heap alloc 2 箇所 (`std::vector<uint32_t> codec_in` と `StreamCtx::decoded` の再 grow) が両方とも消える。both-mode では cod / opt 双方で発火していたためコスト削減幅は単一テンプレ経路の 2 倍見込み。本番ベンチでの elapsed / IPC / L3 ミス変化は別マシン測定待ち
- 残課題 / 次フェーズへの引き継ぎ: なし。Phase 4a (C5-a) へ進む

### Phase 4a: C5-a — Stage 1 バッチ末 mode1_results flush + capacity 解放
- 日付: 2026-05-17
- 主な変更:
  - `src/search/search_orchestrator.cpp`: Stage 1 バッチループの先頭で関数ローカルな `std::vector<OrchestratorHit> mode1_results` を導入。`run_stage1_jobs` の戻り直後に、`in.config.mode == 1` 限定で `idxs` をなめて `states[ji].mode1_results` を `OrchestratorHit` に変換しながら `mode1_results` に move-append し、直後 `std::vector<ChainResult>().swap(st.mode1_results)` で capacity ごと解放
  - 旧 mode 1 fast path (Stage 2 直前で全 ext_jobs を 1 スレッドで fold していたループ) はバッチ末で既に drain 済みのため、`return mode1_results;` のみに簡略化
  - `test/test_multivolume.cpp` に `test_mode1_batched_equals_single` を追加。`run_search<uint16_t>` を直接呼び、posting_budget を大 (両ボリュームが 1 バッチに収まる) / 1 (各ボリュームが単独 Tier1 バッチに落ちる) で 2 回実行し、結果のキー集合 (qi, vi, volume_index, seq_id, stage1_score) が完全一致することを確認
- テスト: `ctest --test-dir build` 全 158 件 pass。`test_multivolume` の新ケースも pass
- 観察: バッチ末で per-ji の `mode1_results` capacity が即解放されるため、Stage 1 全体での anonymous heap が `posting_budget` 相当に bound される (本番ベンチで観察された SwapPss 5 GiB は本フェーズで消える見込み — 別マシン測定待ち)
- 残課題 / 次フェーズへの引き継ぎ: 現状の per-batch fold は単一スレッド逐次。バッチ数縮小 (Phase 5 後) に伴い 1 バッチ内 candidate 数が増えるため、Phase 4b (parallel_scan) は Phase 5 後に効果を再計測する価値あり

### Phase 4b: C5-b-iii — mode1 fold を tbb::parallel_scan で並列化
- 日付: 2026-05-17
- 主な変更:
  - `src/search/search_orchestrator.cpp`: Phase 4a で導入した per-batch fold を 2 パス並列化:
    1. `tbb::parallel_scan` で `offsets[i] = prefix_sum(states[idxs[i]].mode1_results.size())` を計算 (`offsets[n_batch]` がバッチ合計)
    2. `mode1_results.resize(base + batch_total)` を 1 回行い、`tbb::parallel_for` で各 `i` が `mode1_results[base + offsets[i] + k]` にロック不要で書き込み
    3. 各 `ji` のキャパシティ解放 (`std::vector<ChainResult>().swap(...)`) も同じ parallel_for で実施
    - 全てのバッチサイズで対応 (空バッチや batch_total==0 でも安全に解放)
  - `tbb::task_arena arena` (line 286 で既に構築) 内で `arena.execute([&] { ... })` 経由で発火し、ネストした並列化からの worker fairness を保証
  - `test/test_multivolume.cpp::test_mode1_batched_equals_single` を拡張: 16 query × 2 vol を 4 thread で 2 回連続実行し、ソート済みキー集合が完全一致することを確認 (parallel_scan / parallel_for の決定論性検証)
- テスト: `ctest --test-dir build` 全 158 件 pass
- 観察: per-batch の単一スレッド fold (O(batch_candidates)) が並列に置き換わる。バッチ数が少ない (Phase 5 で 57 → 1 想定) ケースほど fold サイズが大きくなるため、本フェーズの効果は Phase 5 後により顕著になる見込み
- 残課題 / 次フェーズへの引き継ぎ: なし。Phase 4.5 (both-mode integration テスト追加) へ進む

### Phase 4.5: both-mode integration テストの追加
- 日付: 2026-05-17
- 主な変更:
  - `test/test_search_integration.cpp`: ヘルパ `build_both_template_indexes(k, t, ...)` を追加 (cod / opt インデックスペアを `g_test_dir` に build。再ビルド回避のため `.kix` 存在チェックで memoize)
  - 新ケース `test_search_mode1_both_template()` を追加 — mode 1 + both-template の e2e。`search_volume_both` を `config.mode = 1` で実行し、FJ876973.1 が検出され、ChainResult の chain フィールド (`q_start/q_end/s_start/s_end/chainscore`) がすべて 0、`stage1_score` が `[1, Nqkmer]` に収まることを検証
  - 新ケース `test_search_stage1_both_template()` を追加 — `stage1_one_strand_both<uint32_t>` を直接呼出 (fresh `JobState` + `Stage1Buffer`)。FJ876973.1 が candidates に含まれ、cod / opt 双方ヒットしても per-(sid, q_pos) dedup により score が `Nqkmer` を超えないことを検証
  - `setup_ssu_testdata.sh` への追加は不要 (テスト内で cod / opt をビルドする既存パターンを再利用)
- テスト: `ctest --test-dir build` 全 158 件 pass。新規 2 ケースも pass
- 観察: 既存 `test_search_both_template` は mode 2 経路のみカバーしていたため、mode 1 + both と stage1_one_strand_both の直接呼出を今回追加。Phase 5 (C1) 以降の改造で both-mode 経路の同値性を回帰検出できる
- 残課題 / 次フェーズへの引き継ぎ: なし。Phase 5 (C1) へ進む

