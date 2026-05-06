#include "ikafssnhttpd/query_store.hpp"

#include <cerrno>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <utility>

#include <fcntl.h>
#include <unistd.h>

#include "util/zstd_oneshot.hpp"

namespace fs = std::filesystem;

namespace ikafssn {

namespace {

constexpr const char* kSuffix    = ".json.zst";
constexpr const char* kTmpSuffix = ".json.zst.tmp";

} // namespace

QueryStore::QueryStore(std::string root, int level)
    : root_(std::move(root)), level_(level) {}

std::string QueryStore::shard_for(const std::string& job_id) {
    char buf[3] = {'0', '0', '\0'};
    if (job_id.size() >= 1) buf[0] = job_id[0];
    if (job_id.size() >= 2) buf[1] = job_id[1];
    return std::string(buf, 2);
}

std::string QueryStore::path_for(const std::string& job_id) const {
    return root_ + "/" + shard_for(job_id) + "/" + job_id + kSuffix;
}

bool QueryStore::init(std::string& error_msg) {
    std::error_code ec;
    fs::create_directories(root_, ec);
    if (ec) {
        error_msg = "QueryStore::init: mkdir " + root_ + ": " + ec.message();
        return false;
    }
    int n = sweep_tmp(error_msg);
    (void)n;
    return error_msg.empty();
}

bool QueryStore::write_atomic_(const std::string& job_id,
                               const void* compressed, std::size_t n,
                               std::string& error_msg) {
    std::string path = path_for(job_id);
    fs::path p(path);
    std::error_code ec;
    if (p.has_parent_path()) {
        fs::create_directories(p.parent_path(), ec);
        if (ec) {
            error_msg = "QueryStore::write: mkdir "
                      + p.parent_path().string() + ": " + ec.message();
            return false;
        }
    }

    std::string tmp = path + ".tmp";
    int fd = ::open(tmp.c_str(),
                    O_WRONLY | O_CREAT | O_TRUNC | O_CLOEXEC,
                    0644);
    if (fd < 0) {
        error_msg = "QueryStore::write: open " + tmp + ": " + std::strerror(errno);
        return false;
    }

    const uint8_t* p_buf = static_cast<const uint8_t*>(compressed);
    std::size_t remaining = n;
    while (remaining > 0) {
        ssize_t w = ::write(fd, p_buf, remaining);
        if (w < 0) {
            if (errno == EINTR) continue;
            error_msg = "QueryStore::write: write " + tmp + ": "
                      + std::strerror(errno);
            ::close(fd);
            ::unlink(tmp.c_str());
            return false;
        }
        p_buf += w;
        remaining -= static_cast<std::size_t>(w);
    }

#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) \
    || defined(__NetBSD__) || defined(__DragonFly__)
    int sync_rc = ::fsync(fd);
#else
    int sync_rc = ::fdatasync(fd);
#endif
    if (sync_rc != 0 && errno != EINVAL) {
        error_msg = "QueryStore::write: sync " + tmp + ": "
                  + std::strerror(errno);
        ::close(fd);
        ::unlink(tmp.c_str());
        return false;
    }
    if (::close(fd) != 0) {
        error_msg = "QueryStore::write: close " + tmp + ": "
                  + std::strerror(errno);
        ::unlink(tmp.c_str());
        return false;
    }

    if (std::rename(tmp.c_str(), path.c_str()) != 0) {
        error_msg = "QueryStore::write: rename " + tmp + " -> " + path
                  + ": " + std::strerror(errno);
        ::unlink(tmp.c_str());
        return false;
    }
    return true;
}

bool QueryStore::write_compressed_passthrough(const std::string& job_id,
                                              const void* bytes, std::size_t n,
                                              std::string& error_msg) {
    return write_atomic_(job_id, bytes, n, error_msg);
}

bool QueryStore::write_plain(const std::string& job_id,
                             const void* bytes, std::size_t n,
                             std::string& error_msg) {
    std::vector<uint8_t> compressed;
    if (!zstd_compress(bytes, n, compressed, level_, error_msg)) {
        return false;
    }
    return write_atomic_(job_id, compressed.data(), compressed.size(),
                         error_msg);
}

bool QueryStore::read(const std::string& job_id,
                      std::string& out_json,
                      std::string& error_msg) const {
    std::string path = path_for(job_id);
    std::vector<uint8_t> decoded;
    if (!zstd_decompress_file(path, decoded, error_msg)) {
        return false;
    }
    out_json.assign(reinterpret_cast<const char*>(decoded.data()),
                    decoded.size());
    return true;
}

bool QueryStore::exists(const std::string& job_id) const {
    std::error_code ec;
    return fs::exists(path_for(job_id), ec);
}

bool QueryStore::unlink(const std::string& job_id, std::string& error_msg) {
    std::error_code ec;
    fs::remove(path_for(job_id), ec);
    if (ec && ec != std::errc::no_such_file_or_directory) {
        error_msg = "QueryStore::unlink: " + ec.message();
        return false;
    }
    return true;
}

int QueryStore::sweep_tmp(std::string& error_msg) {
    int removed = 0;
    std::error_code ec;
    if (!fs::exists(root_, ec)) return 0;

    fs::recursive_directory_iterator it(root_, ec);
    if (ec) {
        error_msg = "QueryStore::sweep_tmp: open " + root_ + ": "
                  + ec.message();
        return 0;
    }
    const std::string tmp_suffix = kTmpSuffix;
    for (; it != fs::recursive_directory_iterator(); it.increment(ec)) {
        if (ec) {
            error_msg = "QueryStore::sweep_tmp: iterate: " + ec.message();
            return removed;
        }
        if (!it->is_regular_file(ec)) continue;
        std::string name = it->path().filename().string();
        if (name.size() > tmp_suffix.size()
            && name.compare(name.size() - tmp_suffix.size(),
                            tmp_suffix.size(), tmp_suffix) == 0) {
            std::error_code rc;
            fs::remove(it->path(), rc);
            if (!rc) removed++;
        }
    }
    return removed;
}

std::vector<std::string> QueryStore::list_job_ids(std::string& error_msg) const {
    std::vector<std::string> ids;
    std::error_code ec;
    if (!fs::exists(root_, ec)) return ids;

    fs::recursive_directory_iterator it(root_, ec);
    if (ec) {
        error_msg = "QueryStore::list_job_ids: open " + root_ + ": "
                  + ec.message();
        return ids;
    }
    const std::string suffix = kSuffix;
    const std::string tmp_suffix = kTmpSuffix;
    for (; it != fs::recursive_directory_iterator(); it.increment(ec)) {
        if (ec) {
            error_msg = "QueryStore::list_job_ids: iterate: " + ec.message();
            return ids;
        }
        if (!it->is_regular_file(ec)) continue;
        std::string name = it->path().filename().string();
        // Skip *.json.zst.tmp.
        if (name.size() > tmp_suffix.size()
            && name.compare(name.size() - tmp_suffix.size(),
                            tmp_suffix.size(), tmp_suffix) == 0) {
            continue;
        }
        if (name.size() <= suffix.size()) continue;
        if (name.compare(name.size() - suffix.size(),
                         suffix.size(), suffix) != 0) {
            continue;
        }
        ids.emplace_back(name.substr(0, name.size() - suffix.size()));
    }
    return ids;
}

} // namespace ikafssn
