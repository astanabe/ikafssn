#include "ikafssnhttpd/result_store.hpp"

#include <cerrno>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <utility>

#include <fcntl.h>
#include <unistd.h>

#include <zstd.h>

#include "util/zstd_oneshot.hpp"

namespace fs = std::filesystem;

namespace ikafssn {

ResultStore::ResultStore(std::string root, int level)
    : root_(std::move(root)), level_(level) {}

std::string ResultStore::shard_for(const std::string& job_id) {
    // 2-char prefix sharding.  Pad with '0' if the job_id is shorter
    // than two characters (defensive: real UUIDv4 ids are 36 chars).
    char buf[3] = {'0', '0', '\0'};
    if (job_id.size() >= 1) buf[0] = job_id[0];
    if (job_id.size() >= 2) buf[1] = job_id[1];
    return std::string(buf, 2);
}

std::string ResultStore::path_for(const std::string& job_id) const {
    return root_ + "/" + shard_for(job_id) + "/" + job_id + ".bin.zst";
}

bool ResultStore::init(std::string& error_msg) {
    std::error_code ec;
    fs::create_directories(root_, ec);
    if (ec) {
        error_msg = "ResultStore::init: mkdir " + root_ + ": " + ec.message();
        return false;
    }
    int n = sweep_tmp(error_msg);
    (void)n;
    return error_msg.empty();
}

bool ResultStore::write(const std::string& job_id,
                        const void* bytes, std::size_t n,
                        std::string& error_msg) {
    // Compress to a memory buffer first, then write atomically with an
    // explicit fdatasync to harden against crash-after-rename (rename
    // is atomic but the file's data blocks may not yet be on disk).
    std::vector<uint8_t> compressed;
    if (!zstd_compress(bytes, n, compressed, level_, error_msg)) {
        return false;
    }

    std::string path = path_for(job_id);
    fs::path p(path);
    std::error_code ec;
    if (p.has_parent_path()) {
        fs::create_directories(p.parent_path(), ec);
        if (ec) {
            error_msg = "ResultStore::write: mkdir "
                      + p.parent_path().string() + ": " + ec.message();
            return false;
        }
    }

    std::string tmp = path + ".tmp";
    int fd = ::open(tmp.c_str(),
                    O_WRONLY | O_CREAT | O_TRUNC | O_CLOEXEC,
                    0644);
    if (fd < 0) {
        error_msg = "ResultStore::write: open " + tmp + ": " + std::strerror(errno);
        return false;
    }

    const uint8_t* p_buf = compressed.data();
    std::size_t remaining = compressed.size();
    while (remaining > 0) {
        ssize_t w = ::write(fd, p_buf, remaining);
        if (w < 0) {
            if (errno == EINTR) continue;
            error_msg = "ResultStore::write: write " + tmp + ": "
                      + std::strerror(errno);
            ::close(fd);
            ::unlink(tmp.c_str());
            return false;
        }
        p_buf += w;
        remaining -= static_cast<std::size_t>(w);
    }

    if (::fdatasync(fd) != 0 && errno != EINVAL) {
        // EINVAL on filesystems that don't support fdatasync — ignore.
        error_msg = "ResultStore::write: fdatasync " + tmp + ": "
                  + std::strerror(errno);
        ::close(fd);
        ::unlink(tmp.c_str());
        return false;
    }
    if (::close(fd) != 0) {
        error_msg = "ResultStore::write: close " + tmp + ": "
                  + std::strerror(errno);
        ::unlink(tmp.c_str());
        return false;
    }

    if (std::rename(tmp.c_str(), path.c_str()) != 0) {
        error_msg = "ResultStore::write: rename " + tmp + " -> " + path
                  + ": " + std::strerror(errno);
        ::unlink(tmp.c_str());
        return false;
    }
    return true;
}

bool ResultStore::exists(const std::string& job_id) const {
    std::error_code ec;
    return fs::exists(path_for(job_id), ec);
}

bool ResultStore::unlink(const std::string& job_id, std::string& error_msg) {
    std::error_code ec;
    fs::remove(path_for(job_id), ec);
    if (ec && ec != std::errc::no_such_file_or_directory) {
        error_msg = "ResultStore::unlink: " + ec.message();
        return false;
    }
    return true;
}

int ResultStore::sweep_tmp(std::string& error_msg) {
    int removed = 0;
    std::error_code ec;
    if (!fs::exists(root_, ec)) return 0;

    fs::recursive_directory_iterator it(root_, ec);
    if (ec) {
        error_msg = "ResultStore::sweep_tmp: open " + root_ + ": "
                  + ec.message();
        return 0;
    }
    for (; it != fs::recursive_directory_iterator(); it.increment(ec)) {
        if (ec) {
            error_msg = "ResultStore::sweep_tmp: iterate: " + ec.message();
            return removed;
        }
        if (!it->is_regular_file(ec)) continue;
        std::string name = it->path().filename().string();
        if (name.size() > 4 && name.compare(name.size() - 4, 4, ".tmp") == 0
            && name.find(".bin.zst.tmp") != std::string::npos) {
            std::error_code rc;
            fs::remove(it->path(), rc);
            if (!rc) removed++;
        }
    }
    return removed;
}

std::vector<std::string> ResultStore::list_job_ids(std::string& error_msg) const {
    std::vector<std::string> ids;
    std::error_code ec;
    if (!fs::exists(root_, ec)) return ids;

    fs::recursive_directory_iterator it(root_, ec);
    if (ec) {
        error_msg = "ResultStore::list_job_ids: open " + root_ + ": "
                  + ec.message();
        return ids;
    }
    for (; it != fs::recursive_directory_iterator(); it.increment(ec)) {
        if (ec) {
            error_msg = "ResultStore::list_job_ids: iterate: " + ec.message();
            return ids;
        }
        if (!it->is_regular_file(ec)) continue;
        std::string name = it->path().filename().string();
        // Match exactly *.bin.zst (not *.bin.zst.tmp).
        const std::string suffix = ".bin.zst";
        if (name.size() <= suffix.size()) continue;
        if (name.compare(name.size() - suffix.size(), suffix.size(), suffix) != 0)
            continue;
        if (name.size() > 4
            && name.compare(name.size() - 4, 4, ".tmp") == 0)
            continue;
        ids.emplace_back(name.substr(0, name.size() - suffix.size()));
    }
    return ids;
}

} // namespace ikafssn
