#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <dirent.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>

#include "ikafssnclient/checkpoint.hpp"
#include "io/result_writer.hpp"
#include "protocol/messages.hpp"
#include "util/logger.hpp"

using namespace ikafssn;

static std::string test_dir;

static void setup_test_dir() {
    char tmpl[] = "/tmp/ikafssn_ckpt_test_XXXXXX";
    char* dir = ::mkdtemp(tmpl);
    assert(dir != nullptr);
    test_dir = dir;
}

static void remove_recursive(const std::string& path) {
    struct stat st;
    if (::stat(path.c_str(), &st) != 0) return;
    if (S_ISDIR(st.st_mode)) {
        DIR* d = ::opendir(path.c_str());
        if (d) {
            struct dirent* ent;
            while ((ent = ::readdir(d)) != nullptr) {
                std::string name(ent->d_name);
                if (name == "." || name == "..") continue;
                remove_recursive(path + "/" + name);
            }
            ::closedir(d);
        }
        ::rmdir(path.c_str());
    } else {
        ::unlink(path.c_str());
    }
}

static void cleanup_test_dir() {
    if (!test_dir.empty()) {
        remove_recursive(test_dir);
    }
}

// ---- SHA256 tests ----

static void test_sha256_string() {
    std::cout << "test_sha256_string... ";
    // SHA256("hello") = 2cf24dba5fb0a30e26e83b2ac5b9e29e1b161e5c1fa7425e73043362938b9824
    std::string hash = sha256_string("hello");
    assert(hash == "2cf24dba5fb0a30e26e83b2ac5b9e29e1b161e5c1fa7425e73043362938b9824");

    // Empty string
    std::string empty_hash = sha256_string("");
    assert(empty_hash == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");

    std::cout << "OK\n";
}

static void test_sha256_file() {
    std::cout << "test_sha256_file... ";
    std::string path = test_dir + "/sha256_test.txt";
    {
        std::ofstream out(path);
        out << "hello";
    }
    std::string hash = sha256_file(path);
    assert(hash == "2cf24dba5fb0a30e26e83b2ac5b9e29e1b161e5c1fa7425e73043362938b9824");

    // Non-existent file
    std::string empty = sha256_file(test_dir + "/nonexistent");
    assert(empty.empty());

    std::cout << "OK\n";
}

// ---- build_options_text tests ----

static void test_build_options_text_deterministic() {
    std::cout << "test_build_options_text_deterministic... ";

    SearchRequest req;
    req.k = 9;
    req.mode = 2;
    req.db = "testdb";
    req.stage1_score = 1;

    DbStats stats;
    stats.db_name = "testdb";
    stats.total_sequences = 1000;
    stats.total_bases = 500000;

    std::string text1 = build_options_text(req, stats, 9, OutputFormat::kTab,
                                            "abc123", "");
    std::string text2 = build_options_text(req, stats, 9, OutputFormat::kTab,
                                            "abc123", "");
    assert(text1 == text2);
    assert(!text1.empty());

    // Different parameter produces different text
    req.mode = 3;
    std::string text3 = build_options_text(req, stats, 9, OutputFormat::kTab,
                                            "abc123", "");
    assert(text3 != text1);

    std::cout << "OK\n";
}

// ---- LockGuard tests ----

static void test_lock_mechanism() {
    std::cout << "test_lock_mechanism... ";
    std::string lock_dir = test_dir + "/test_lock";

    // Acquire lock
    {
        LockGuard guard(lock_dir);
        assert(guard.locked());

        // Double acquire should fail
        LockGuard guard2(lock_dir);
        assert(!guard2.locked());
    }
    // After scope exit, lock should be released

    // Re-acquire should succeed
    LockGuard guard3(lock_dir);
    assert(guard3.locked());
    guard3.release();

    std::cout << "OK\n";
}

// ---- Checkpoint round-trip test ----

static void test_checkpoint_roundtrip() {
    std::cout << "test_checkpoint_roundtrip... ";

    Logger logger(Logger::kDebug);

    // Create a temporary input file
    std::string input_path = test_dir + "/input.fasta";
    {
        std::ofstream out(input_path);
        out << ">seq1\nACGTACGT\n>seq2\nTTTTAAAA\n>seq3\nGGGGCCCC\n";
    }

    SearchRequest req;
    req.k = 9;
    req.mode = 2;
    req.db = "testdb";
    req.stage1_score = 1;

    DbStats stats;
    stats.db_name = "testdb";
    stats.total_sequences = 100;
    stats.total_bases = 50000;

    std::string options_text = build_options_text(req, stats, 9,
                                                   OutputFormat::kTab, "", "");
    std::string input_sha = sha256_file(input_path);

    // Use test_dir as working directory prefix
    Checkpoint::Config cfg;
    cfg.output_path = test_dir + "/output.txt";
    cfg.input_path = input_path;
    cfg.ix_name = "testdb";
    cfg.resolved_k = 9;
    cfg.outfmt = OutputFormat::kTab;

    Checkpoint ckpt(cfg, logger);

    // Initialize
    assert(!ckpt.exists());
    assert(ckpt.initialize(options_text, input_sha, ""));
    assert(ckpt.exists());

    // Acquire lock
    LockGuard lock;
    assert(ckpt.acquire_lock(lock));

    // Write batch 0
    std::vector<std::string> batch0_seqids = {"seq1", "seq2"};
    assert(ckpt.write_batch_seqids(0, batch0_seqids));

    std::vector<OutputHit> batch0_hits;
    {
        OutputHit h;
        h.qseqid = "seq1";
        h.sseqid = "sub1";
        h.sstrand = '+';
        h.qstart = 0; h.qend = 8; h.qlen = 8;
        h.sstart = 10; h.send = 18; h.slen = 100;
        h.coverscore = 5; h.chainscore = 10;
        h.volume = 0;
        batch0_hits.push_back(h);
    }
    assert(ckpt.write_batch_results(0, batch0_hits, 2, 1, false));

    // Write response meta
    assert(ckpt.write_response_meta(2, 1, false));

    // Release lock and simulate resume
    lock.release();

    // Create new checkpoint instance (simulating restart)
    Checkpoint ckpt2(cfg, logger);
    assert(ckpt2.exists());

    LockGuard lock2;
    assert(ckpt2.acquire_lock(lock2));

    std::unordered_set<std::string> completed;
    int next_batch = 0;
    assert(ckpt2.resume(options_text, input_sha, completed, next_batch));
    assert(completed.size() == 2);
    assert(completed.count("seq1") == 1);
    assert(completed.count("seq2") == 1);
    assert(next_batch == 1);

    // Write batch 1
    std::vector<std::string> batch1_seqids = {"seq3"};
    assert(ckpt2.write_batch_seqids(1, batch1_seqids));
    std::vector<OutputHit> batch1_hits;
    assert(ckpt2.write_batch_results(1, batch1_hits, 2, 1, false));

    // Read response meta
    uint8_t mode, s1score;
    bool traceback;
    assert(ckpt2.read_response_meta(mode, s1score, traceback));
    assert(mode == 2);
    assert(s1score == 1);
    assert(!traceback);

    // Merge results
    assert(ckpt2.merge_results(cfg.output_path, 2, 1, false));

    // Verify output file exists and has content
    {
        std::ifstream in(cfg.output_path);
        assert(in.is_open());
        std::string line;
        assert(std::getline(in, line)); // header
        assert(line[0] == '#');
        assert(std::getline(in, line)); // batch 0 hit
        assert(line.find("seq1") != std::string::npos);
    }

    // Cleanup
    lock2.release();
    ckpt2.cleanup();
    assert(!ckpt2.exists());

    // Remove output file
    ::unlink(cfg.output_path.c_str());

    std::cout << "OK\n";
}

// ---- Tab merge test (header dedup) ----

static void test_tab_merge_header_dedup() {
    std::cout << "test_tab_merge_header_dedup... ";

    Logger logger(Logger::kInfo);

    Checkpoint::Config cfg;
    cfg.output_path = test_dir + "/tab_merge_output.txt";
    cfg.input_path = test_dir + "/dummy.fasta";
    cfg.ix_name = "db";
    cfg.resolved_k = 9;
    cfg.outfmt = OutputFormat::kTab;

    // Create dummy input
    {
        std::ofstream out(cfg.input_path);
        out << ">q1\nACGT\n>q2\nTTTT\n";
    }

    SearchRequest req;
    req.k = 9;
    req.mode = 2;
    req.db = "db";

    DbStats stats;
    stats.db_name = "db";
    stats.total_sequences = 10;
    stats.total_bases = 1000;

    std::string options_text = build_options_text(req, stats, 9,
                                                   OutputFormat::kTab, "", "");
    std::string input_sha = sha256_file(cfg.input_path);

    Checkpoint ckpt(cfg, logger);
    assert(ckpt.initialize(options_text, input_sha, ""));

    // Write two batches with hits
    {
        OutputHit h;
        h.qseqid = "q1"; h.sseqid = "s1"; h.sstrand = '+';
        h.qstart = 0; h.qend = 4; h.qlen = 4;
        h.sstart = 0; h.send = 4; h.slen = 100;
        h.coverscore = 3; h.chainscore = 5; h.volume = 0;
        std::vector<OutputHit> hits = {h};
        ckpt.write_batch_seqids(0, {"q1"});
        ckpt.write_batch_results(0, hits, 2, 1, false);
    }
    {
        OutputHit h;
        h.qseqid = "q2"; h.sseqid = "s2"; h.sstrand = '+';
        h.qstart = 0; h.qend = 4; h.qlen = 4;
        h.sstart = 0; h.send = 4; h.slen = 200;
        h.coverscore = 2; h.chainscore = 4; h.volume = 0;
        std::vector<OutputHit> hits = {h};
        ckpt.write_batch_seqids(1, {"q2"});
        ckpt.write_batch_results(1, hits, 2, 1, false);
    }

    // Merge
    assert(ckpt.merge_results(cfg.output_path, 2, 1, false));

    // Read and verify: should have exactly one header line
    {
        std::ifstream in(cfg.output_path);
        std::string content((std::istreambuf_iterator<char>(in)),
                             std::istreambuf_iterator<char>());
        int header_count = 0;
        int data_count = 0;
        std::istringstream iss(content);
        std::string line;
        while (std::getline(iss, line)) {
            if (line.empty()) continue;
            if (line[0] == '#') header_count++;
            else data_count++;
        }
        assert(header_count == 1);
        assert(data_count == 2);
    }

    ckpt.cleanup();
    ::unlink(cfg.output_path.c_str());
    ::unlink(cfg.input_path.c_str());

    std::cout << "OK\n";
}

// ---- JSON fragment merge test ----

static void test_json_fragment_merge() {
    std::cout << "test_json_fragment_merge... ";

    Logger logger(Logger::kInfo);

    Checkpoint::Config cfg;
    cfg.output_path = test_dir + "/json_merge_output.json";
    cfg.input_path = test_dir + "/dummy2.fasta";
    cfg.ix_name = "db";
    cfg.resolved_k = 9;
    cfg.outfmt = OutputFormat::kJson;

    // Create dummy input
    {
        std::ofstream out(cfg.input_path);
        out << ">q1\nACGT\n>q2\nTTTT\n";
    }

    SearchRequest req;
    req.k = 9;
    req.mode = 2;
    req.db = "db";

    DbStats stats;
    stats.db_name = "db";
    stats.total_sequences = 10;
    stats.total_bases = 1000;

    std::string options_text = build_options_text(req, stats, 9,
                                                   OutputFormat::kJson, "", "");
    std::string input_sha = sha256_file(cfg.input_path);

    Checkpoint ckpt(cfg, logger);
    assert(ckpt.initialize(options_text, input_sha, ""));

    // Write two batches
    {
        OutputHit h;
        h.qseqid = "q1"; h.sseqid = "s1"; h.sstrand = '+';
        h.qstart = 0; h.qend = 4; h.qlen = 4;
        h.sstart = 0; h.send = 4; h.slen = 100;
        h.coverscore = 3; h.chainscore = 5; h.volume = 0;
        std::vector<OutputHit> hits = {h};
        ckpt.write_batch_seqids(0, {"q1"});
        ckpt.write_batch_results(0, hits, 2, 1, false);
    }
    {
        OutputHit h;
        h.qseqid = "q2"; h.sseqid = "s2"; h.sstrand = '+';
        h.qstart = 0; h.qend = 4; h.qlen = 4;
        h.sstart = 0; h.send = 4; h.slen = 200;
        h.coverscore = 2; h.chainscore = 4; h.volume = 0;
        std::vector<OutputHit> hits = {h};
        ckpt.write_batch_seqids(1, {"q2"});
        ckpt.write_batch_results(1, hits, 2, 1, false);
    }

    // Merge
    assert(ckpt.merge_results(cfg.output_path, 2, 1, false));

    // Verify JSON structure
    {
        std::ifstream in(cfg.output_path);
        std::string content((std::istreambuf_iterator<char>(in)),
                             std::istreambuf_iterator<char>());
        // Should start with {"results": [ and end with ]}
        assert(content.find("{") != std::string::npos);
        assert(content.find("\"results\"") != std::string::npos);
        assert(content.find("\"q1\"") != std::string::npos);
        assert(content.find("\"q2\"") != std::string::npos);
        // Should not have trailing comma before ]
        auto close_bracket = content.rfind(']');
        assert(close_bracket != std::string::npos);
        // Find last non-whitespace before ]
        size_t pos = close_bracket - 1;
        while (pos > 0 && (content[pos] == ' ' || content[pos] == '\n'))
            pos--;
        assert(content[pos] != ',');
    }

    ckpt.cleanup();
    ::unlink(cfg.output_path.c_str());
    ::unlink(cfg.input_path.c_str());

    std::cout << "OK\n";
}

// ---- resolve_db_stats test ----

static void test_resolve_db_stats() {
    std::cout << "test_resolve_db_stats... ";

    InfoResponse info;
    DatabaseInfo db;
    db.name = "nt";
    db.default_k = 11;
    KmerGroupInfo grp;
    grp.k = 11;
    grp.kmer_type = 1;
    VolumeInfo v1;
    v1.volume_index = 0;
    v1.num_sequences = 100;
    v1.total_bases = 50000;
    VolumeInfo v2;
    v2.volume_index = 1;
    v2.num_sequences = 200;
    v2.total_bases = 80000;
    grp.volumes = {v1, v2};
    db.groups = {grp};
    info.databases = {db};

    DbStats stats = resolve_db_stats(info, "nt", 11);
    assert(stats.db_name == "nt");
    assert(stats.total_sequences == 300);
    assert(stats.total_bases == 130000);

    // Non-matching k
    DbStats stats2 = resolve_db_stats(info, "nt", 9);
    assert(stats2.total_sequences == 0);

    // Non-matching db
    DbStats stats3 = resolve_db_stats(info, "other", 11);
    assert(stats3.total_sequences == 0);

    std::cout << "OK\n";
}

int main() {
    setup_test_dir();

    test_sha256_string();
    test_sha256_file();
    test_build_options_text_deterministic();
    test_lock_mechanism();
    test_checkpoint_roundtrip();
    test_tab_merge_header_dedup();
    test_json_fragment_merge();
    test_resolve_db_stats();

    cleanup_test_dir();

    std::cout << "\nAll checkpoint tests passed!\n";
    return 0;
}
