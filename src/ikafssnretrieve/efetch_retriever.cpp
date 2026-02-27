#include "ikafssnretrieve/efetch_retriever.hpp"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <sstream>

#ifdef IKAFSSN_ENABLE_REMOTE
#include <chrono>
#include <thread>
#include <unordered_map>
#include <curl/curl.h>
#endif

namespace ikafssn {

static const char* EFETCH_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi";

std::string build_efetch_url_batch(const std::vector<std::string>& accessions,
                                   const std::string& api_key) {
    std::ostringstream url;
    url << EFETCH_BASE << "?db=nuccore&rettype=fasta&retmode=text&id=";
    for (size_t i = 0; i < accessions.size(); i++) {
        if (i > 0) url << ',';
        url << accessions[i];
    }
    if (!api_key.empty()) {
        url << "&api_key=" << api_key;
    }
    return url.str();
}

std::string build_efetch_url_range(const std::string& accession,
                                   uint32_t seq_start, uint32_t seq_stop,
                                   const std::string& api_key) {
    std::ostringstream url;
    url << EFETCH_BASE << "?db=nuccore&rettype=fasta&retmode=text&id=" << accession
        << "&seq_start=" << seq_start
        << "&seq_stop=" << seq_stop;
    if (!api_key.empty()) {
        url << "&api_key=" << api_key;
    }
    return url.str();
}

// Extract accession from an efetch FASTA defline.
// Deflines look like: >AB123456.1 Description text
// or: >gi|12345|gb|AB123456.1| Description
// We extract the versioned accession.
static std::string extract_accession_from_defline(const std::string& defline) {
    if (defline.empty() || defline[0] != '>') return {};

    // Skip '>'
    size_t start = 1;
    while (start < defline.size() && std::isspace(static_cast<unsigned char>(defline[start])))
        start++;

    // Check for gi| format: >gi|...|gb|ACCESSION.V|
    if (defline.substr(start, 3) == "gi|") {
        // Find the accession part after db|
        size_t pipe_count = 0;
        size_t pos = start;
        while (pos < defline.size() && pipe_count < 3) {
            if (defline[pos] == '|') pipe_count++;
            pos++;
        }
        if (pipe_count == 3) {
            size_t acc_end = defline.find('|', pos);
            if (acc_end == std::string::npos)
                acc_end = defline.find(' ', pos);
            if (acc_end == std::string::npos)
                acc_end = defline.size();
            std::string acc = defline.substr(pos, acc_end - pos);
            // Strip version suffix for matching
            auto dot = acc.find('.');
            if (dot != std::string::npos) {
                return acc.substr(0, dot);
            }
            return acc;
        }
    }

    // Simple format: >ACCESSION.V Description
    size_t end = start;
    while (end < defline.size() && !std::isspace(static_cast<unsigned char>(defline[end])))
        end++;
    std::string acc = defline.substr(start, end - start);
    // Strip version suffix
    auto dot = acc.find('.');
    if (dot != std::string::npos) {
        return acc.substr(0, dot);
    }
    return acc;
}

std::vector<EfetchRecord> parse_efetch_response(const std::string& response) {
    std::vector<EfetchRecord> records;
    std::istringstream iss(response);
    std::string line;
    std::string cur_acc;
    std::string cur_seq;

    auto finish_record = [&]() {
        if (!cur_acc.empty() && !cur_seq.empty()) {
            // Convert to uppercase
            for (auto& c : cur_seq)
                c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
            records.push_back({std::move(cur_acc), std::move(cur_seq)});
        }
        cur_acc.clear();
        cur_seq.clear();
    };

    while (std::getline(iss, line)) {
        if (!line.empty() && line.back() == '\r')
            line.pop_back();

        if (line.empty()) continue;

        if (line[0] == '>') {
            finish_record();
            cur_acc = extract_accession_from_defline(line);
        } else {
            cur_seq += line;
        }
    }
    finish_record();

    return records;
}

uint32_t rate_limit_sleep_ms(bool has_api_key) {
    return has_api_key ? 100 : 334;
}

bool is_retryable_http_status(long status_code) {
    return status_code == 429 || status_code == 503;
}

bool is_skip_http_status(long status_code) {
    return status_code == 400 || status_code == 404;
}

#ifdef IKAFSSN_ENABLE_REMOTE

// libcurl write callback
static size_t write_callback(char* ptr, size_t size, size_t nmemb, void* userdata) {
    auto* buf = static_cast<std::string*>(userdata);
    size_t total = size * nmemb;
    buf->append(ptr, total);
    return total;
}

// Perform HTTP GET with retries and rate limiting.
static bool http_get(const std::string& url, std::string& response_body,
                     long& http_status, const EfetchOptions& opts) {
    CURL* curl = curl_easy_init();
    if (!curl) {
        std::fprintf(stderr, "efetch: curl_easy_init failed\n");
        return false;
    }

    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_callback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response_body);
    curl_easy_setopt(curl, CURLOPT_TIMEOUT, static_cast<long>(opts.timeout_sec));
    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(curl, CURLOPT_USERAGENT, "ikafssn/0.1");

    uint32_t backoff_ms = 1000;

    for (uint32_t attempt = 0; attempt <= opts.retries; attempt++) {
        response_body.clear();
        CURLcode res = curl_easy_perform(curl);

        if (res != CURLE_OK) {
            std::fprintf(stderr, "efetch: request failed: %s\n", curl_easy_strerror(res));
            if (attempt < opts.retries) {
                std::fprintf(stderr, "efetch: retrying in %u ms (attempt %u/%u)\n",
                             backoff_ms, attempt + 1, opts.retries);
                std::this_thread::sleep_for(std::chrono::milliseconds(backoff_ms));
                backoff_ms *= 2;
                continue;
            }
            curl_easy_cleanup(curl);
            return false;
        }

        curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &http_status);

        if (http_status == 200) {
            curl_easy_cleanup(curl);
            return true;
        }

        if (is_retryable_http_status(http_status) && attempt < opts.retries) {
            std::fprintf(stderr, "efetch: HTTP %ld, retrying in %u ms (attempt %u/%u)\n",
                         http_status, backoff_ms, attempt + 1, opts.retries);
            std::this_thread::sleep_for(std::chrono::milliseconds(backoff_ms));
            backoff_ms *= 2;
            continue;
        }

        // Non-retryable or out of retries
        curl_easy_cleanup(curl);
        return false;
    }

    curl_easy_cleanup(curl);
    return false;
}

// Reverse complement a DNA string in-place.
static void reverse_complement(std::string& seq) {
    std::reverse(seq.begin(), seq.end());
    for (auto& c : seq) {
        switch (c) {
            case 'A': c = 'T'; break;
            case 'T': c = 'A'; break;
            case 'C': c = 'G'; break;
            case 'G': c = 'C'; break;
        }
    }
}

// Write a single FASTA record with 70-char line wrapping.
static void write_fasta_record(std::ostream& out, const std::string& header,
                               const std::string& sequence) {
    out << '>' << header << '\n';
    for (size_t i = 0; i < sequence.size(); i += 70) {
        size_t len = std::min<size_t>(70, sequence.size() - i);
        out.write(sequence.data() + i, static_cast<std::streamsize>(len));
        out << '\n';
    }
}

uint32_t retrieve_remote(const std::vector<OutputHit>& hits,
                         const EfetchOptions& opts,
                         std::ostream& out) {
    if (hits.empty()) return 0;

    bool has_key = !opts.api_key.empty();
    uint32_t sleep_ms = rate_limit_sleep_ms(has_key);

    if (!has_key) {
        std::fprintf(stderr, "WARNING: No NCBI API key set. "
                     "Rate limited to 3 requests/sec. "
                     "Set -api_key or NCBI_API_KEY for higher throughput.\n");
    }

    // Group hits: separate into batch-eligible (short seqs) and individual (long seqs).
    // For batch: group by accession to avoid duplicate fetches; extract locally from fetched full seq.
    // For individual: fetch only the needed range.

    struct HitRef {
        size_t hit_index;
        uint32_t ext_start;
        uint32_t ext_end;
    };

    // Collect all hits per accession and determine fetch strategy
    struct AccessionInfo {
        std::vector<HitRef> hit_refs;
        uint32_t max_s_end = 0;  // for sequence length estimation
    };
    std::unordered_map<std::string, AccessionInfo> acc_info;
    for (size_t i = 0; i < hits.size(); i++) {
        auto& info = acc_info[hits[i].accession];
        uint32_t ext_start = hits[i].s_start;
        uint32_t ext_end = hits[i].s_end;
        if (opts.context > 0) {
            ext_start = (ext_start >= opts.context) ? ext_start - opts.context : 0;
            ext_end += opts.context;  // may exceed actual seq length; efetch handles this
        }
        info.hit_refs.push_back({i, ext_start, ext_end});
        if (hits[i].s_end > info.max_s_end)
            info.max_s_end = hits[i].s_end;
    }

    // Separate into batch and individual lists
    std::vector<std::string> batch_accessions;
    std::vector<std::string> individual_accessions;

    for (const auto& [acc, info] : acc_info) {
        // Estimate sequence length from max s_end (add margin)
        uint32_t est_len = info.max_s_end + 1;
        if (est_len > opts.range_threshold) {
            individual_accessions.push_back(acc);
        } else {
            batch_accessions.push_back(acc);
        }
    }

    uint32_t retrieved = 0;

    // --- Batch retrieval ---
    for (size_t i = 0; i < batch_accessions.size(); i += opts.batch_size) {
        size_t end = std::min(i + opts.batch_size, batch_accessions.size());
        std::vector<std::string> batch(batch_accessions.begin() + i,
                                       batch_accessions.begin() + end);

        std::string url = build_efetch_url_batch(batch, opts.api_key);
        std::string response;
        long http_status = 0;

        if (!http_get(url, response, http_status, opts)) {
            if (is_skip_http_status(http_status)) {
                std::fprintf(stderr, "efetch: HTTP %ld for batch, skipping %zu accessions\n",
                             http_status, batch.size());
            } else {
                std::fprintf(stderr, "efetch: batch request failed (HTTP %ld)\n", http_status);
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(sleep_ms));
            continue;
        }

        // Parse response and match to hits
        auto records = parse_efetch_response(response);
        std::unordered_map<std::string, std::string> fetched_seqs;
        for (auto& rec : records) {
            fetched_seqs[rec.accession] = std::move(rec.sequence);
        }

        for (const auto& acc : batch) {
            auto seq_it = fetched_seqs.find(acc);
            if (seq_it == fetched_seqs.end()) {
                std::fprintf(stderr, "efetch: accession '%s' not in response\n", acc.c_str());
                continue;
            }
            const auto& full_seq = seq_it->second;
            uint32_t actual_len = static_cast<uint32_t>(full_seq.size());

            for (const auto& hr : acc_info[acc].hit_refs) {
                const auto& hit = hits[hr.hit_index];
                uint32_t ext_start = hr.ext_start;
                uint32_t ext_end = std::min(hr.ext_end, actual_len > 0 ? actual_len - 1 : 0);
                if (ext_start > ext_end) continue;

                std::string subseq = full_seq.substr(ext_start, ext_end - ext_start + 1);
                if (hit.strand == '-') {
                    reverse_complement(subseq);
                }

                std::ostringstream header;
                header << hit.accession
                       << " query=" << hit.query_id
                       << " strand=" << hit.strand
                       << " range=" << ext_start << '-' << ext_end
                       << " score=" << hit.chainscore;
                write_fasta_record(out, header.str(), subseq);
                retrieved++;
            }
        }

        // Rate limit sleep
        std::this_thread::sleep_for(std::chrono::milliseconds(sleep_ms));
    }

    // --- Individual retrieval (long sequences) ---
    for (const auto& acc : individual_accessions) {
        const auto& info = acc_info[acc];

        for (const auto& hr : info.hit_refs) {
            const auto& hit = hits[hr.hit_index];
            // efetch uses 1-based inclusive coordinates
            uint32_t seq_start = hr.ext_start + 1;
            uint32_t seq_stop = hr.ext_end + 1;

            std::string url = build_efetch_url_range(acc, seq_start, seq_stop, opts.api_key);
            std::string response;
            long http_status = 0;

            if (!http_get(url, response, http_status, opts)) {
                if (is_skip_http_status(http_status)) {
                    std::fprintf(stderr, "efetch: HTTP %ld for '%s', skipping\n",
                                 http_status, acc.c_str());
                } else {
                    std::fprintf(stderr, "efetch: individual request failed for '%s' (HTTP %ld)\n",
                                 acc.c_str(), http_status);
                }
                std::this_thread::sleep_for(std::chrono::milliseconds(sleep_ms));
                continue;
            }

            auto records = parse_efetch_response(response);
            if (records.empty()) {
                std::fprintf(stderr, "efetch: no sequence in response for '%s'\n", acc.c_str());
                std::this_thread::sleep_for(std::chrono::milliseconds(sleep_ms));
                continue;
            }

            std::string subseq = std::move(records[0].sequence);
            if (hit.strand == '-') {
                reverse_complement(subseq);
            }

            std::ostringstream header;
            header << hit.accession
                   << " query=" << hit.query_id
                   << " strand=" << hit.strand
                   << " range=" << hr.ext_start << '-' << hr.ext_end
                   << " score=" << hit.chainscore;
            write_fasta_record(out, header.str(), subseq);
            retrieved++;

            // Rate limit sleep
            std::this_thread::sleep_for(std::chrono::milliseconds(sleep_ms));
        }
    }

    return retrieved;
}

#endif // IKAFSSN_ENABLE_REMOTE

} // namespace ikafssn
