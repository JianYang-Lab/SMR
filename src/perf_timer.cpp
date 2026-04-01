#include "perf_timer.hpp"

#include <spdlog/spdlog.h>
#include <sys/resource.h>

/// PerfTimer
PerfTimer::PerfTimer(const std::string& name)
    : name_(name),
      start_time_(std::chrono::high_resolution_clock::now()),
      elapsed_time_(start_time_),
      stopped_(false) {}

PerfTimer::~PerfTimer() {
  if (!stopped_) {
    stop();
  }
}

void PerfTimer::elapsed(const std::string& elapsed_name) {
  if (stopped_) return;  // Avoid double stop
  auto result = elapsed_to(elapsed_time_);
  std::string time_str = result.first;
  elapsed_time_ = result.second;
  spdlog::info("[perf] ===> elapsed {}, cost {}", elapsed_name, time_str);
}

void PerfTimer::stop() {
  if (stopped_) return;  // Avoid double stop
  auto result = elapsed_to(start_time_);
  std::string time_str = result.first;
  spdlog::info("[perf] ===> {}, cost {}", name_, time_str);
  stopped_ = true;
}

std::pair<std::string, std::chrono::high_resolution_clock::time_point> PerfTimer::elapsed_to(
    const std::chrono::high_resolution_clock::time_point& time_point) {
  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = end_time - time_point;
  std::chrono::duration<double, std::milli> milliseconds = duration;
  std::chrono::duration<double> seconds = duration;
  std::chrono::duration<double, std::ratio<60>> minutes = duration;
  std::chrono::duration<double, std::ratio<3600>> hours = duration;

  std::string time_str;
  if (milliseconds.count() < 1000.0) {
    time_str = fmt::format("{}ms", milliseconds.count());
  } else if (seconds.count() < 60.0) {
    time_str = fmt::format("{:02}:{:02}:{:02.3f}s", 0, 0, seconds.count());
  } else if (minutes.count() < 60) {
    auto min = std::trunc(minutes.count());
    auto sec = seconds.count() - 60.0 * min;
    time_str = fmt::format("{:02}:{:02}:{:02.3f}s", 0, min, sec);
  } else {
    auto h = std::trunc(hours.count());
    auto min = std::trunc(minutes.count()) - 60 * h;
    auto sec = seconds.count() - 60 * std::trunc(minutes.count());
    time_str = fmt::format("{:02}:{:02}:{:02.3f}s", h, min, sec);
  }

  return {time_str, end_time};
}

size_t get_memory_usage() {
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  return usage.ru_maxrss * 1024;  // Returns bytes
}
