#ifndef PARALLEL_SOLVER_HPP
#define PARALLEL_SOLVER_HPP

#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <functional>
#include <queue>
#include <memory>
#include <cmath>

// ==========================================
// Minimal ThreadPool
// ==========================================
class ThreadPool {
public:
    explicit ThreadPool(int n) {
        for (int i = 0; i < n; ++i)
            workers_.emplace_back([this] {
                while (true) {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lk(mtx_);
                        cv_.wait(lk, [this] { return stop_ || !tasks_.empty(); });
                        if (stop_ && tasks_.empty()) return;
                        task = std::move(tasks_.front());
                        tasks_.pop();
                    }
                    task();
                }
            });
    }
    ~ThreadPool() {
        { std::lock_guard<std::mutex> lk(mtx_); stop_ = true; }
        cv_.notify_all();
        for (auto& w : workers_) if (w.joinable()) w.join();
    }
    int size() const { return (int)workers_.size(); }

    template<typename F>
    void parallelFor(int start, int end, F&& fn) {
        int n = end - start;
        if (n <= 0) return;
        int nt = (int)workers_.size();
        int chunk = std::max(1, n / nt);
        std::atomic<int> next{start};
        std::vector<std::thread> threads;
        for (int t = 0; t < nt; ++t) {
            threads.emplace_back([&] {
                while (true) {
                    int i = next.fetch_add(chunk);
                    if (i >= end) break;
                    int e = std::min(i + chunk, end);
                    for (int j = i; j < e; ++j) fn(j);
                }
            });
        }
        for (auto& t : threads) if (t.joinable()) t.join();
    }

private:
    std::vector<std::thread> workers_;
    std::queue<std::function<void()>> tasks_;
    std::mutex mtx_;
    std::condition_variable cv_;
    bool stop_ = false;
};

// ==========================================
// ParallelIncrementalNullspaceSolver
// Same dense RREF but parallel elimination
// ==========================================
template<typename T, typename NullspaceInfo>
class ParallelIncrementalNullspaceSolver {
public:
    explicit ParallelIncrementalNullspaceSolver(size_t nv, int nt = 0)
        : num_vars_(nv) {
        int threads = (nt > 0) ? nt : (int)std::thread::hardware_concurrency();
        if (threads < 1) threads = 1;
        pool_ = std::make_unique<ThreadPool>(threads);
    }
    
    void setNumThreads(int n) {
        int threads = (n > 0) ? n : (int)std::thread::hardware_concurrency();
        if (threads < 1) threads = 1;
        pool_ = std::make_unique<ThreadPool>(threads);
    }
    
    void addRow(const std::vector<T>& row) {
        if (rank_ >= (int)num_vars_) return;
        std::vector<T> r = row;
        
        if (rank_ >= 50 && pool_->size() > 1)
            forwardEliminateParallel(r);
        else
            forwardEliminateSequential(r);
        
        int np = -1;
        for (size_t j = 0; j < num_vars_; ++j) {
            if (!isZero(r[j])) { np = (int)j; break; }
        }
        if (np < 0) return;
        
        T pv = r[np];
        for (size_t j = np; j < num_vars_; ++j) r[j] = r[j] / pv;
        
        if (rank_ >= 50 && pool_->size() > 1)
            backSubstituteParallel(r, np);
        else
            backSubstituteSequential(r, np);
        
        int pos = 0;
        while (pos < rank_ && pivot_cols_[pos] < np) ++pos;
        rref_.insert(rref_.begin() + pos, std::move(r));
        pivot_cols_.insert(pivot_cols_.begin() + pos, np);
        rank_++;
        cache_valid_ = false;
    }
    
    void addRows(const std::vector<std::vector<T>>& rows) {
        for (const auto& r : rows) addRow(r);
    }
    
    int getNullity() const { return (int)num_vars_ - rank_; }
    size_t getNumVars() const { return num_vars_; }
    size_t getNumRows() const { return rref_.size(); }
    
    NullspaceInfo getNullspace() {
        if (!cache_valid_) { cached_info_ = computeNullspace(); cache_valid_ = true; }
        return cached_info_;
    }
    
    void clear() { rref_.clear(); pivot_cols_.clear(); rank_ = 0; cache_valid_ = false; }
    
    T verifyBasis(const std::vector<std::vector<T>>& rows,
                  const NullspaceInfo& info, T tol = T(1e-10)) {
        if (!info.is_valid || info.basis.empty()) return T(0);
        T mr = T(0);
        for (const auto& row : rows) {
            for (const auto& b : info.basis) {
                T dot = T(0);
                for (size_t i = 0; i < num_vars_; ++i) dot += row[i] * b[i];
                T ad = absVal(dot);
                if constexpr (std::is_floating_point_v<T>) {
                    if (ad > mr) mr = ad;
                } else { if (ad != T(0)) return T(1); }
            }
        }
        return mr;
    }

private:
    size_t num_vars_;
    std::vector<std::vector<T>> rref_;
    std::vector<int> pivot_cols_;
    int rank_ = 0;
    NullspaceInfo cached_info_;
    bool cache_valid_ = false;
    std::unique_ptr<ThreadPool> pool_;
    
    static bool isZero(const T& v) {
        if constexpr (std::is_floating_point_v<T>) return std::abs(v) < T(1e-15);
        else return v == T(0);
    }
    static T absVal(const T& v) {
        if constexpr (std::is_floating_point_v<T>) return std::abs(v);
        else return v;
    }
    
    void forwardEliminateSequential(std::vector<T>& r) {
        for (int i = 0; i < rank_; ++i) {
            int pc = pivot_cols_[i];
            if (isZero(r[pc])) continue;
            T f = r[pc];
            for (size_t j = pc; j < num_vars_; ++j)
                r[j] = r[j] - f * rref_[i][j];
        }
    }
    
    void forwardEliminateParallel(std::vector<T>& r) {
        int nt = pool_->size();
        int chunk = std::max(1, rank_ / nt);
        std::vector<std::vector<T>> deltas(nt, std::vector<T>(num_vars_, T(0)));
        std::vector<std::thread> threads;
        for (int t = 0; t < nt; ++t) {
            int start = t * chunk;
            int end = (t == nt - 1) ? rank_ : start + chunk;
            threads.emplace_back([&, t, start, end]() {
                auto& d = deltas[t];
                for (int i = start; i < end; ++i) {
                    int pc = pivot_cols_[i];
                    T factor = r[pc];
                    if (isZero(factor)) continue;
                    const auto& rr = rref_[i];
                    for (size_t j = pc; j < num_vars_; ++j)
                        d[j] = d[j] - factor * rr[j];
                }
            });
        }
        for (auto& th : threads) if (th.joinable()) th.join();
        for (int t = 0; t < nt; ++t) {
            const auto& d = deltas[t];
            for (size_t j = 0; j < num_vars_; ++j)
                if (!isZero(d[j])) r[j] = r[j] + d[j];
        }
    }
    
    void backSubstituteSequential(const std::vector<T>& norm, int np) {
        for (int i = 0; i < rank_; ++i) {
            auto& row = rref_[i];
            if (isZero(row[np])) continue;
            T f = row[np];
            for (size_t j = np; j < num_vars_; ++j)
                row[j] = row[j] - f * norm[j];
        }
    }
    
    void backSubstituteParallel(const std::vector<T>& norm, int np) {
        int nt = pool_->size();
        int chunk = std::max(1, rank_ / nt);
        std::vector<std::thread> threads;
        for (int t = 0; t < nt; ++t) {
            int start = t * chunk;
            int end = (t == nt - 1) ? rank_ : start + chunk;
            threads.emplace_back([&, start, end, np]() {
                for (int i = start; i < end; ++i) {
                    auto& row = rref_[i];
                    if (isZero(row[np])) continue;
                    T f = row[np];
                    for (size_t j = np; j < num_vars_; ++j)
                        row[j] = row[j] - f * norm[j];
                }
            });
        }
        for (auto& th : threads) if (th.joinable()) th.join();
    }
    
    NullspaceInfo computeNullspace() {
        NullspaceInfo info;
        int n = (int)num_vars_;
        info.rank = rank_;
        info.nullity = n - rank_;
        if (info.nullity <= 0) { info.is_valid = (info.nullity == 0); return info; }
        std::vector<bool> ip(n, false);
        for (int pc : pivot_cols_) ip[pc] = true;
        for (int j = 0; j < n; ++j) {
            if (ip[j]) continue;
            std::vector<T> vec(n, T(0));
            vec[j] = T(1);
            for (int r = 0; r < rank_; ++r) {
                int pc = pivot_cols_[r];
                if (!isZero(rref_[r][j]))
                    vec[pc] = vec[pc] - rref_[r][j];
            }
            info.basis.push_back(std::move(vec));
            info.free_cols.push_back(j);
        }
        info.is_valid = true;
        return info;
    }
};

#endif
