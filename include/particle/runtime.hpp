#pragma once
#include <algorithm>
#include <atomic>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>
namespace particle {
inline size_t worker_count(int configured, size_t jobs,
                           size_t automatic_limit = 6) {
  if (!jobs)
    return 0;
  size_t requested =
      configured
          ? size_t(configured)
          : std::min<size_t>(automatic_limit,
                             std::max(1u, std::thread::hardware_concurrency()));
  return std::min(jobs, requested);
}
template <class F> auto parallel_map(size_t jobs, size_t workers, F function) {
  using Result = decltype(function(size_t(0)));
  std::vector<Result> result(jobs);
  if (!jobs)
    return result;
  std::atomic<size_t> next{0};
  std::atomic<bool> failed{false};
  std::exception_ptr error;
  std::mutex mutex;
  auto worker = [&] {
    while (!failed) {
      size_t index = next.fetch_add(1);
      if (index >= jobs)
        return;
      try {
        result[index] = function(index);
      } catch (...) {
        std::lock_guard lock(mutex);
        if (!error)
          error = std::current_exception();
        failed = true;
      }
    }
  };
  {
    std::vector<std::jthread> threads;
    for (size_t i = 1; i < std::min(jobs, workers); ++i)
      threads.emplace_back(worker);
    worker();
  }
  if (error)
    std::rethrow_exception(error);
  return result;
}
} // namespace particle
