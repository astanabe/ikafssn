# ikafssn SIMD ベクトル化計画 概要

本書は、ikafssn (C++20、Linux x86_64 / aarch64) に SIMD ベクトル化を導入するための **概要計画書** である。次セッションではコンテキストをクリアした状態でこのファイルを起点に詳細計画を立て、実装を進める。

---

## 0. 進捗状況 (live)

| Phase | 状態 | ブランチ / コミット | 詳細計画ファイル |
|---|---|---|---|
| **Phase 0** — 基盤整備 | **完了 (2026-04-27)** | `vectorize` 515a8b3 | `/home/shimotsuki/.claude/plans/vectorize-plan-md-ikafssn-kind-key.md` |
| **Phase 1 §0 訂正** (Google Benchmark 取り込み) | **完了 (2026-04-27)** | `vectorize` 165a19a | `/home/shimotsuki/.claude/plans/vectorize-plan-md-phase-0-home-shimotsu-imperative-alpaca.md §0` |
| **Phase 1a** — H-4 (FASTA toupper) | **完了 (2026-04-27)** | `vectorize` 56b8afb | 同上 §1 |
| **Phase 1b** — H-1 (ncbi2na unpack) | **完了 (2026-04-27)** | `vectorize` 909fbff | 同上 §2 |
| **Phase 2** — H-3 (Stage 1 SoA) + M-4 (parallel_sort) | 未着手 | — | `/home/shimotsuki/.claude/plans/vectorize-phase-2-stage1-soa-parallel-sort.md` |
| Phase 3 — H-2 (varint バッチ) + M-1 (k-mer revcomp) | 未着手 | — | 未作成 |
| Phase 4 — M-2 (degenerate scan) + M-3 (spaced seed 並列) | 未着手 | — | 未作成 |

### 次セッションでの作業起点

1. **本ファイル (`vectorize_plan.md`) を読む** — 計画全体像
2. **`vectorize` ブランチに切り替える** — `git checkout vectorize` (Phase 0/1 commit を含む。Phase 1b commit は `909fbff`)
3. **次に取り組む Phase の詳細計画** を本書 §4 のロードマップから抜き出して、新しいプランファイル (`.md`) を作成しユーザに提示する
4. **過去の Phase 詳細計画を参照する場合**:
   - Phase 0: `/home/shimotsuki/.claude/plans/vectorize-plan-md-ikafssn-kind-key.md`
   - Phase 1: `/home/shimotsuki/.claude/plans/vectorize-plan-md-phase-0-home-shimotsu-imperative-alpaca.md`
   - Phase 2: `/home/shimotsuki/.claude/plans/vectorize-phase-2-stage1-soa-parallel-sort.md`
   - 以前 Phase で確定した API・テストマクロ (`SimdCap` enum、`IKAFSSN_FORCE_SIMD`、`add_ikafssn_simd_test`、`check_required_tier_or_skip()`、Google Benchmark fixture `apply_force_tier`/`label_tier`) は次 Phase でも流用する
5. **計画立案は Plan モードでユーザの明示的指示後に行う** (CLAUDE.md の Plan Mode Rules 準拠)

### Phase 0 で確定した API (Phase 1 以降が依存)

`src/util/simd_dispatch.hpp` で定義済み。詳細仕様は Phase 0 詳細計画 §0.1-0.2 を参照。

- `SimdCap` enum (アーキ別値域分離: x86=10..50, arm=100..140)
- `init_simd_dispatch(Logger*)` — `std::call_once` で起動時 1 回呼ぶ。7 main で実装済み
- `current_simd_cap()` / `auto_detected_simd_cap()` / `current_simd_flags()` (slow_bmi2 / has_avx512f_only)
- `force_simd_cap(SimdCap)` — downgrade のみ、超過要求は Scalar に降格
- `parse_force_simd_env(const char*)` — `IKAFSSN_FORCE_SIMD` env パーサ (大小・区切り無視)
- `check_required_tier_or_skip()` — `IKAFSSN_TEST_REQUIRE_TIER` 非対応で `exit(77)`
- `reset_simd_dispatch_for_testing()` — header 側のみ `IKAFSSN_TESTING` ガード、cpp は常時コンパイル

### Phase 2 以降で考慮すべき積み残し決定事項

以下は Phase 0/1 で据え置いた項目、または Phase 1 で確定した方針。Phase 2 以降の詳細計画で再判断する。

- SIMD 実装の配置 (Phase 1 で確定: 新規ファイル分離 — `src/io/text_simd.{hpp,cpp}`、`src/core/ncbi2na_unpack.{hpp,cpp}`。Phase 2 以降も同方針を継続)
- Google Benchmark 導入 — **完了** (Phase 1 §0 訂正パッチ `165a19a` で `IKAFSSN_BUILD_BENCH` option + `FetchContent v1.9.4` を取り込み済み)
- `.kix`/`.kpx` v4 (Stream VByte) — Phase 3a 効果測定後に判断
- AVX-512 CI 対応 (self-hosted / Intel SDE / 開発者ローカル) — Phase 2 完了時に再評価

### 検証環境メモ (Phase 0 時点)

- 開発機: Ryzen 9 9950X (Zen5, family 0x1A) — `auto_detected_simd_cap = AVX512VBMI2`、`slow_bmi2 = false`
- CTest: 48/48 pass、x86 6 tier 実行 + aarch64 3 tier skip
- aarch64 実機検証は Phase 1 完了後に Graviton / Apple Silicon で別途実施

---

## 1. 目的とスコープ

- **目的**: ikafssn の主要ホットパスをアーキテクチャ別 SIMD で加速し、index 構築時間と検索レイテンシを短縮する
- **対象アーキ**: x86_64 (SSE4.2 〜 AVX-512 VBMI2)、aarch64 (NEON 〜 SME2)
- **対象外**:
  - **AMX / SME / SME2** — 整数タイル積に特化したアクセラレータ。ikafssn の主要ホットパス (bit 操作・prefix sum・gather) と噛み合わない。Stage 3 を 8-bit 量子化 batch alignment に置換する場合にのみ将来検討
  - **FMA** — 浮動小数点積和。ikafssn は整数・ビット演算中心のため不要 (Parasail は元々 SIMD 化済み・本計画の対象外)
  - **AVX512F のみ (BW/VBMI 非対応) CPU** — 該当は Knights Landing / Knights Mill のみで EOL。byte 粒度処理が中心の ikafssn では AVX2 経路で十分
  - **AVX512VNNI** (`vpdpbusd` / `vpdpwssd`) — 量子化ニューラルネット推論専用 (int8/int16 の dense 整数積和)。ikafssn の処理構造 (bit 操作・byte シャッフル・gather/scatter 中心) と完全に噛み合わない。スコアリングも単純加算で積和構造を持たない
  - **AVX512VP2INTERSECT** (`vp2intersectd` / `vp2intersectq`) — 集合交差専用設計だが、唯一の出荷シリコン Tiger Lake で**マイクロコード実装 (~40 cycle)** 、`vpbroadcastd` + `vpcmpeqd` の手動エミュレーション (10〜20 cycle) より遅い。Intel は Sapphire Rapids 以降で削除しており将来性なし。集合交差が必要な箇所 (将来の OID filter SIMD 化等) は **broadcast + compare + `vpcompressd` でエミュレート** する

---

## 2. ディスパッチ機構 — `pg_kmersearch` 方式の採用

`../pg_kmersearch` で採用されている **起動時 CPU 検出 + ラッパー関数による if-else 分岐ディスパッチ** 方式をそのまま ikafssn に適用する。

### 2.1 採用理由

- C++ でも `__attribute__((target("avx2")))` 等が g++/clang++ で完全サポートされ、C 言語実装と同じ書き方が成立する
- 単一バイナリで複数 ISA レベルに対応でき、配布・パッケージング (conda/homebrew) が単純化される
- ifunc / `target_clones` (FMV) より単純で、デバッグ容易・インライン展開を阻害しない
- 関数ポインタ経由ではなく **ループ外でディスパッチ → 内部は SIMD ベタ書き** とすることで per-iteration の分岐コストをゼロにできる

### 2.2 検出と初期化

- `src/util/simd_dispatch.{hpp,cpp}` を新設
- `src/util/common_init.hpp` の起動シーケンス (各バイナリの `main()` 冒頭で呼ばれる) から `init_simd_dispatch()` を呼ぶ
- `std::once_flag` で初期化の冪等性を保証 (テストから複数回呼ばれても安全)
- 検出方式:
  - **x86_64**: `__builtin_cpu_supports("avx2")` 等。古いコンパイラ向け補助に CPUID 直接読み (leaf 7, ECX bit 1 = VBMI、bit 6 = VBMI2)
  - **aarch64**: `getauxval(AT_HWCAP)` / `AT_HWCAP2` (Linux 限定なのでシグナル試行は不要)
- 検出結果は `extern SimdCap g_simd_cap;` に格納し、初期化後は **読み取り専用**

### 2.3 SimdCap 列挙型 (アーキ別に値域分離)

```cpp
enum class SimdCap : int {
    Scalar      = 0,

    // ── x86_64 系 (10〜99) ──
    SSE42       = 10,  // SSE4.2 + POPCNT             (Nehalem, 2008-)
    AVX2        = 20,  // AVX2 + FMA + BMI1 + BMI2    (Haswell, 2013-)
    AVX512BW    = 30,  // AVX-512 F+BW+DQ+CD+VL       (Skylake-X, 2017-)
    AVX512VBMI  = 40,  // + VBMI + IFMA               (Ice Lake client, 2019-)
    AVX512VBMI2 = 50,  // + VBMI2+BITALG+VPOPCNTDQ+GFNI (Ice Lake server, 2021-)

    // ── aarch64 系 (100〜199) ──
    NEON        = 100, // Armv8.0 ベースライン
    SVE         = 110, // Armv8.2-A + SVE             (A64FX, Neoverse V1/N2)
    SVE2        = 120, // Armv9-A                      (Neoverse V2, Cortex-X2+)
    SME         = 130, // Scalable Matrix Extension   (Apple M4, 一部 Armv9.2) — 当面不使用
    SME2        = 140, //                                                          — 当面不使用
};
```

- x64 binary では `Scalar` または `10〜50` のみ、ARM64 binary では `Scalar` または `100〜140` のみが入る
- ifdef ガード内でのみ `>=` 比較するため、アーキ跨ぎの誤比較は構造的に発生しない
- **AVX2 tier には BMI1/BMI2 を同梱** (CPU 世代が一致するため。ただし BMI2 は AMD Zen1/2 の `pdep`/`pext` が遅いためベンダ判定で scalar fallback)
- **AVX512F のみ (BW なし) は SSE42 の次に上げず、AVX2 と同等扱い** (KNL/KNM 専用、市場ほぼ消滅、byte 粒度処理で利点なし)

### 2.4 BMI2 ブラックリスト (重要)

- AMD Zen1 / Zen2 の `pdep` / `pext` はマイクロコード実装で約 100 サイクル
- これらの世代で BMI2 経路を選ぶと scalar より遅くなる
- `init_simd_dispatch()` 内で CPUID で vendor=AMD かつ family ≦ Zen2 を判定し、AVX2 tier で `slow_pext=true` フラグを立て、BMI2 を必須とする実装は scalar 経路を選択する
- Zen3 以降 (`family >= 0x19`) と Intel 全般は問題なし

### 2.5 ラッパー関数のテンプレート

```cpp
__attribute__((target("avx512vbmi2,avx512vbmi,avx512bw,avx512f,avx2,bmi2")))
void unpack_ncbi2na_avx512vbmi2(const uint8_t* in, uint8_t* out, size_t n);

__attribute__((target("avx512vbmi,avx512bw,avx512f,avx2,bmi2")))
void unpack_ncbi2na_avx512vbmi (const uint8_t* in, uint8_t* out, size_t n);

__attribute__((target("avx512bw,avx512f,avx2,bmi2")))
void unpack_ncbi2na_avx512bw   (const uint8_t* in, uint8_t* out, size_t n);

__attribute__((target("avx2,bmi2")))
void unpack_ncbi2na_avx2_bmi2  (const uint8_t* in, uint8_t* out, size_t n);

void unpack_ncbi2na_scalar     (const uint8_t* in, uint8_t* out, size_t n);

void unpack_ncbi2na(const uint8_t* in, uint8_t* out, size_t n) {
#ifdef __x86_64__
    if (g_simd_cap >= SimdCap::AVX512VBMI2 && n >= 512) { unpack_ncbi2na_avx512vbmi2(in,out,n); return; }
    if (g_simd_cap >= SimdCap::AVX512VBMI  && n >= 256) { unpack_ncbi2na_avx512vbmi (in,out,n); return; }
    if (g_simd_cap >= SimdCap::AVX512BW    && n >= 128) { unpack_ncbi2na_avx512bw   (in,out,n); return; }
    if (g_simd_cap >= SimdCap::AVX2        && n >=  64) { unpack_ncbi2na_avx2_bmi2  (in,out,n); return; }
#elif defined(__aarch64__)
    if (g_simd_cap >= SimdCap::SVE2        && n >=  64) { unpack_ncbi2na_sve2(in,out,n); return; }
    if (g_simd_cap >= SimdCap::SVE         && n >=  64) { unpack_ncbi2na_sve (in,out,n); return; }
    if (g_simd_cap >= SimdCap::NEON        && n >=  32) { unpack_ncbi2na_neon(in,out,n); return; }
#endif
    unpack_ncbi2na_scalar(in, out, n);
}
```

- **サイズ閾値併用**: 小さい入力では SIMD setup コストが回収できないため最小サイズを設ける
- **テンプレート関数には `target` 属性を付けない** (NTTP 特殊化と相性が悪いため)。SIMD 経路は非テンプレートのプレーン関数として書き、テンプレート由来の scalar 経路と並置する
- **target 属性は累積指定** (`avx512bw,avx512f,avx2,bmi2`) しておくと安心 (コンパイラが上位機能から下位機能を自動有効化しないケースに対応)
- **ホットループでは関数呼び出し境界を跨がない**: Stage 1 のような数千万回ループ内では、ループ外で関数全体を SIMD 化したものをディスパッチする (per-call ではなく per-task で分岐)

### 2.6 ビルドシステム (CMake)

- 全 SIMD 実装を 1 つの `.cpp` ファイル内に `__attribute__((target))` で個別指定 → CMake 側のファイル単位フラグ付与は不要
- ベース ISA はコンパイラ既定 (x86_64 → SSE2、aarch64 → NEON) のままにし、`-march=native` 等は **使わない** (配布バイナリの ABI 互換性確保)
- runtime 強制用 env (`IKAFSSN_FORCE_SIMD=scalar` 等) を実装し、テストで全経路を検証可能にする

---

## 3. ベクトル化候補一覧 (優先度順)

### 優先度: 高

#### H-1. ncbi2na パック形式の 2-bit 塩基抽出 (sliding scan)

- **場所**: `src/core/packed_kmer_scanner.hpp:16` (`ncbi2na_base_at`)、同 `PackedKmerScanner::scan` / `scan_spaced` (約 80〜400 行)
- **処理**: ncbi2na (MSB-first 4 base/byte) から連続 2-bit 塩基を抽出し、k-mer を sliding window 更新
- **想定 SIMD**:
  - x86_64: AVX-512 `vpermb` (VBMI) で 64 base 一括展開、AVX2 `vpshufb` で 32 base
  - aarch64: SVE2 `tbx` + `lsr`、NEON `tbl` で 16 base 一括
- **想定 speedup**: 4〜16x (lane 数比例)。index 構築全体で 1.5〜3x
- **注意**: 32〜64 base ずつ一時バッファに展開し、k-mer 構築 (`kmer = (kmer<<2 | base) & mask`) は既存ループを流用。ambig 範囲は scalar fallback。tail handling も scalar。

#### H-2. LEB128 (varint) バッチデコード

- **場所**: `src/core/varint.hpp:22` (`varint_decode`)、利用元 = `src/index/kix_reader.cpp` posting decode、`src/search/stage1_filter.cpp` の `SeqIdDecoder`、`src/index/kpx_reader.cpp`
- **処理**: delta 圧縮された posting 列の varint 展開。Stage 1 ホットパス。
- **想定 SIMD**:
  - x86_64: **AVX-512 VBMI2** の `vpcompressb` + `vpmovb2m` で複数 varint 同時展開 (Stream VByte 系手法)。VBMI2 が決定的に効く tier
  - aarch64: SVE2 `bext` + `match`
- **想定 speedup**: 2〜5x。Stage 1 全体で 1.3〜1.8x
- **注意**: 本格対応は `.kix`/`.kpx` フォーマット v4 (Stream VByte) への移行も視野。既存 v3 互換性を保つため reader 内に両経路を持つ。

#### H-3. delta デコード (prefix sum) と Stage 1 集計

- **場所**: `src/search/stage1_filter.cpp` Stage 1 メインループ、`src/index/kix_reader.cpp` posting prefix sum
- **処理**: delta → 絶対値復元 (prefix sum)、`entries[sid].score++` (last_pos 比較条件付き)
- **想定 SIMD**:
  - x86_64: prefix sum = AVX-512 Sklansky pattern。集計 = `vpgatherdd` → `vpcmpeqd` → masked `vpaddd` → `vpscatterdd` (+ 衝突検出 `vpconflictd`)
  - aarch64: SVE2 gather/scatter (`ld1d`/`st1d`) + 述語化マスク
- **想定 speedup**: prefix sum 単独 2〜4x、集計含む 1.3〜1.8x
- **注意**: 現状の `Stage1Entry<Tier>` AoS を **SoA (`score[]`, `last_pos[]`)** に変更する設計変更を伴う。tier 切替 (T8→T16→T32 昇格) ロジックは scalar 維持。同一 sid 衝突は AVX-512 CD で検出。

#### H-4. FASTA / シーケンス入力の toupper・改行スキャン

- **場所**: `src/io/fasta_reader.cpp:16` 付近 (toupper)、同 `:30` 以降 (line scanning)
- **処理**: ASCII テキストの大文字化、ヘッダー `>` 検出、改行除去
- **想定 SIMD**:
  - x86_64: AVX2 `vpcmpgtb` で `'a'..'z'` 範囲判定 → `vpsubb` で 0x20 引き (LUT 不要)。AVX-512 BW で 64 byte/cycle
  - aarch64: NEON `vcgeq_u8` + `vsubq_u8`、SVE2 `match` 命令
- **想定 speedup**: 8〜16x (純粋 byte 処理)
- **注意**: ASCII 限定であることを明示 (ロケール非依存化)。改行検出は `_mm256_movemask_epi8` + `tzcnt`、または既に AVX2 化されている glibc `memchr('\n')` ベースで足りる場合あり。

### 優先度: 中

#### M-1. K-mer reverse complement のバッチ処理

- **場所**: `src/core/kmer_encoding.hpp:49` (`kmer_revcomp`)、利用元 = `src/search/query_preprocessor.cpp`
- **処理**: 2-bit packed k-mer (uint16/32) のビット反転 + 2-bit ペア順序反転 + シフト
- **想定 SIMD**:
  - x86_64: AVX2 `vpshufb` で 16 個同時、AVX-512 VBMI `vpermb` + `vpternlogd`
  - aarch64: NEON `rbit` + `vrev16q_u8` (ほぼ即値)、SVE2 述語化版
- **想定 speedup**: 4〜8x (バッチ化前提)
- **注意**: query は数十〜数千 k-mer 程度。閾値 (例: 16 k-mer 以上) で SIMD 経路に切替。spaced seed モードの文字列 RC (`reverse_complement_string`) も別経路で SIMD 化候補。

#### M-2. degenerate (曖昧) 塩基スキャン

- **場所**: `src/core/kmer_encoding.hpp:103` (`contains_degenerate_base`)
- **処理**: 文字列を走査し `R/Y/S/W/K/M/B/D/H/V/N` を検出
- **想定 SIMD**:
  - x86_64: AVX2 `vpcmpeqb` の OR 集約 + `_mm256_movemask_epi8`、AVX-512 `vpermi2b` で LUT lookup
  - aarch64: **SVE2 `match` 命令** (16-byte set membership 専用設計、本処理に直球)
- **想定 speedup**: 8〜16x

#### M-3. spaced seed の複数マスク並列抽出

- **場所**: `src/core/spaced_seed.hpp:177` (`extract_for_mask`, `extract_kmer_ct`)、`PackedKmerScanner::scan_spaced`
- **処理**: 64-bit accumulator から複数マスクで bit gather
- **想定 SIMD**:
  - x86_64: BMI2 scalar `pext` (Intel/Zen3+ で高速)、AVX-512 VBMI2 `vpcompressb` で SIMD 化
  - aarch64: **SVE2 `bdep` / `bext`** (用途専用命令、Neoverse V2 等で対応)
- **想定 speedup**: template_type=both 時 1.5〜2x、それ以外 1〜1.2x
- **注意**: 現状の compile-time 特殊化 (NTTP) は維持。SIMD 経路は「マスク数が多く効果が見込める CPU」で選択。BMI2 は AMD Zen1/2 ブラックリスト適用。

#### M-4. インデックス構築の `tbb::parallel_sort` SIMD 化

- **場所**: `src/index/index_builder.cpp` の posting buffer ソート
- **処理**: posting buffer を k-mer 値で sort
- **想定 SIMD**: Intel x86-simd-sort (vqsort、AVX-512) を導入。32/64-bit キーで 5〜10x。ARM 版 VQSort もあり。
- **想定 speedup**: ソート単体 3〜10x、index 構築全体 1.1〜1.3x
- **注意**: 外部ライブラリ (BSD-3) を vendoring または submodule 取り込み。`tbb::parallel_sort` のコンパレータ位置に局所置換可能。

### 優先度: 低

| # | 候補 | 場所 | 想定 speedup | 不採用/後回し理由 |
|---|---|---|---|---|
| L-1 | Stage 2 chaining DP | `src/search/stage2_chaining.cpp:50` | 1.2〜1.5x | `dp[j]` 依存チェーンで律速。prefetch + scalar unroll の方が ROI 高 |
| L-2 | diagonal filter ヒストグラム | `src/search/diagonal_filter.cpp:7` | 1.2〜1.5x | hash collision のため SIMD ヒストグラム困難。AVX-512 CD で限定的に可 |
| L-3 | nth_element / 最終ソート | `src/search/stage1_filter.cpp:82` | 1.05〜1.2x | 候補数が小規模 (10〜数千) でオーバーヘッドが見合わない |
| L-4 | OID filter bitset 参照 | `src/search/oid_filter.hpp:33` | — | ランダムアクセス + ヒットあたり 1 回。単独 SIMD 化の意味薄い |
| L-5 | ambig 展開 | `src/core/kmer_encoding.hpp:166` | — | 再帰 branching で SIMD 不適合。固定 max_expansion なら LUT 化が代替案 |
| L-6 | **Stage 2/3 のハミング距離検証フィルタ** (将来候補) | `src/search/stage2_chaining.cpp` / `src/search/stage3_align.cpp` | 1.1〜1.5x (Stage 3 呼び出し削減効果) | 既に得られている candidate (q_pos, s_pos) hit ペアに対し、subject 側の k-mer を取り出して XOR + popcount で塩基不一致数を 1 命令算出。Stage 3 (Parasail) を呼ぶ前のローコスト品質フィルタ、または chaining 時のスコア重み付けに利用。**AVX512VPOPCNTDQ + BITALG** (`vpopcntq` / `vpopcntb`)、aarch64 は NEON `cnt` + `addv`、SVE2 `cnt`。導入する場合は `SimdCap` に `AVX512VBMI2` 上位 tier として `AVX512BITALG` を追加 (Ice Lake server で VBMI2 と同梱だが命名は分離)。**注**: Stage 1 の fuzzy lookup には適用不可 (下記補足参照) |

> **Fuzzy lookup については popcount で解決しない**:
> ハミング距離 D ≤ 1 でも近傍 k-mer 数は k=11 で 34、k=12 で 37、D=2 では数百個に膨れる。direct-address table によるルックアップ回数がそのまま倍増するため、本質的なボトルネックは「距離計算」ではなく「近傍列挙」。popcount は 2 つの既知 k-mer ペアの距離計算には効くが、4^11 = 400 万 k-mer 空間からの近傍探索は加速できない。**ikafssn における fuzzy 機能は既存の spaced seed (discontiguous megablast 風テンプレート) で達成されており**、これは "1" 位置だけ exact 一致を要求するため exact lookup の速度を保つ。L-6 はあくまで「ヒット済みペアの後段検証」用途に限定する。

> **§1 「対象外」で除外した命令の補足**:
> - **AVX512VNNI**: dense int8/int16 積和専用。ikafssn 構造と非整合 → 一切採用しない
> - **AVX512VP2INTERSECT**: 集合交差専用だが Tiger Lake で microcoded (~40 cycle)、Sapphire Rapids 以降削除。集合交差は `vpbroadcastd` + `vpcmpeqd` + `vpcompressd` のエミュレーションで代替する (10〜20 cycle)

---

## 4. 実装ロードマップ

### Phase 0 — 基盤整備 (これがすべての前提) ✅ **完了 (2026-04-27, commit 515a8b3 on `vectorize`)**

詳細計画: `/home/shimotsuki/.claude/plans/vectorize-plan-md-ikafssn-kind-key.md`

1. ✅ `src/util/simd_dispatch.{hpp,cpp}` 新設
   - `SimdCap` enum (アーキ別値域分離: x86=10..50, arm=100..140)
   - `init_simd_dispatch(Logger*)` (CPUID / `__builtin_cpu_supports` / `getauxval`)
   - AMD Zen1/2 BMI2 ブラックリスト (vendor=`AuthenticAMD` + family 0x15/0x16/0x17 判定)
   - `IKAFSSN_FORCE_SIMD` env での経路強制 (downgrade-only、超過要求は Scalar 降格 + WARN)
   - **`extern` グローバルではなくアクセサ関数** (`current_simd_cap()` 等) に変更 (test reset 制御のため)
2. ✅ 7 main から `init_simd_dispatch(&logger)` を `std::call_once` 経由で呼ぶ (`common_init.hpp` への統合は見送り、main 個別構成を尊重)
3. ✅ `add_ikafssn_simd_test` CMake マクロ追加。9 tier (x86 6 + arm 3) variant + `SKIP_RETURN_CODE 77`。`test_simd_dispatch` (新規) と `test_kmer_encoding` (既存) を variant 化
4. ✅ CMake オプション `-DIKAFSSN_ENABLE_SIMD=ON` (default ON) 追加。OFF 時は `[INFO] simd: scalar (build-time disabled)`

### Phase 1 — 低リスク・即効性 (Phase 0 完了後)

5. **H-4** (FASTA toupper / 改行スキャン) — 最小変更で 8〜16x、回帰リスク最小
6. **H-1** (ncbi2na 2-bit unpack) — index 構築の主要ボトルネック解消

### Phase 2 — 設計変更を伴うが効果大

7. **H-3** (Stage 1 prefix sum + SoA 化) — `Stage1Buffer` の SoA 化を伴う
8. **M-4** (parallel_sort の SIMD 化) — 外部ライブラリ導入のみ

### Phase 3 — フォーマット変更含む高リスク

9. **H-2** (varint バッチデコード) — `.kix`/`.kpx` フォーマット v4 (Stream VByte) も視野
10. **M-1** (k-mer revcomp バッチ) — query 前処理の SoA 化

### Phase 4 — 限定効果・特定 CPU 依存

11. **M-2** (degenerate scan、SVE2 `match`)
12. **M-3** (spaced seed 並列マスク、BMI2 + SVE2 `bdep`/`bext`)

---

## 5. 共通の実装方針・注意点

- **アライメント**: SoA 配列は 64 byte (AVX-512 cache line) アライン推奨。`Stage1Buffer` は thread-local なので alloc 戦略の変更のみ
- **テール処理**: ループ末尾は scalar に落とす素朴版で十分 (ベクトル長は十分大きい)
- **little-endian 前提**: CLAUDE.md で明記済み。SIMD コードでも前提として活用
- **コンパイラ前提**: GCC 13.3 (`__builtin_cpu_supports("avx512vbmi")` は GCC 8+、`avx512vbmi2` は GCC 11+ で利用可)。古い場合は CPUID 直接読みで補完
- **計測**: 各 Phase 完了時に `db/SSU_eukaryote_rRNA` を使った既存テスト + マイクロベンチで前後比較を必須化
- **テンプレートと target 属性の分離**: NTTP 特殊化されている `extract_kmer_ct` 等のテンプレートには `target` 属性を付けず、SIMD 経路は非テンプレートの並置関数として実装
- **per-call 分岐の最小化**: ホットループ内で関数呼び出しごとに `g_simd_cap` を見ない設計 (関数全体を SIMD 化してループ外でディスパッチ)
- **ベンダ判定 (BMI2)**: AMD CPUID Vendor=`AuthenticAMD` かつ Family `0x15`/`0x16`/`0x17` (Zen1/Zen2 含む) では BMI2 ベース経路を選ばない
- **AVX512F-only CPU の扱い**: 検出された場合は AVX2 tier に降格 (`SimdCap::AVX2` を割り当て)、AVX512BW 経路には進ませない

---

## 6. アーキテクチャ別 命令の使いどころサマリ

| 用途 | x86_64 best | aarch64 best |
|---|---|---|
| 2-bit unpack (H-1) | AVX-512 VBMI `vpermb` | SVE2 `tbx` / NEON `tbl` |
| varint バッチデコード (H-2) | **AVX-512 VBMI2 `vpcompressb` + `vpmovb2m`** | SVE2 `bext` + `match` |
| prefix sum (H-3) | AVX-512 `vpermd` 経由 Sklansky | SVE2 述語化 reduce |
| gather/scatter 集計 (H-3) | AVX-512 `vpgatherdd` + `vpscatterdd` (+ `vpconflictd`) | SVE2 `ld1d`/`st1d` gather/scatter |
| toupper / scan (H-4, M-2) | AVX-512 BW `vpcmpgtb` + `vpsubb` | **SVE2 `match`**、NEON `vcgeq_u8` |
| bit-permutation (M-3) | BMI2 `pext` / `pdep` (Intel/Zen3+) | **SVE2 `bdep` / `bext`** |
| 比較ソート (M-4) | x86-simd-sort (AVX-512) | ARM 版 VQSort |

---

## 7. 次セッションへの引き継ぎ事項

次セッション開始時、コンテキストクリア後の作業手順:

1. **本ファイル (`vectorize_plan.md`) を最初に読む** — これが計画の起点。§0 進捗状況で次に取り組む Phase を確認
2. **CLAUDE.md を読む** — プロジェクト規約 (日本語応答、IKAFSSN_VERSION 更新ルール、Plan モード規約) を再確認
3. **`vectorize` ブランチに切り替える** — `git checkout vectorize`。Phase 0 の commit `515a8b3` が含まれている
4. **次 Phase の詳細計画を作成する** — 本書 §4 ロードマップの当該 Phase 節に書かれた要点をベースに、実装着手可能な粒度の `.md` プランファイルを新規作成しユーザに提示する
   - 過去 Phase の詳細計画ファイルは `~/.claude/plans/` 配下にある (Phase 0: `/home/shimotsuki/.claude/plans/vectorize-plan-md-ikafssn-kind-key.md`)。Phase 0 で確定した API・テストマクロを参照可能
   - Plan モードへの進入は **ユーザの明示的指示があってから** (CLAUDE.md の Plan Mode Rules 準拠)
5. **実装は Phase 単位で進める** — 各 Phase 完了時:
   - ✅ `IKAFSSN_FORCE_SIMD=scalar` で前 Phase マージ時点と bit-exact 一致を CTest で確認
   - ✅ 実 CPU が対応する全 SIMD variant が scalar と bit-exact 一致 (9950X 機なら x86 6 tier)
   - ✅ ベンチ受け入れ基準を満たす
   - ✅ `IKAFSSN_VERSION` を本日付に更新してコミット (CLAUDE.md Version Management 準拠、ソースコード変更時のみ)
   - ✅ 本ファイル §0 進捗表を更新
6. **Phase 0 で据え置いた検討事項** は §0 の「Phase 1 以降で考慮すべき積み残し決定事項」を参照 (SIMD 実装配置、Google Benchmark 導入、`.kix`/`.kpx` v4 判断、AVX-512 CI 対応)
