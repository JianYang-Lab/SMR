#pragma once

#include <chrono>
#include <cstddef>
#include <string>
#include <utility>

class PerfTimer {
 public:
  explicit PerfTimer(std::string name = "");
  ~PerfTimer();

  void elapsed(const std::string& elapsed_name);
  void stop();

 private:
  std::pair<std::string, std::chrono::high_resolution_clock::time_point> elapsed_to(
      const std::chrono::high_resolution_clock::time_point& time_point);

 private:
  std::string name_;
  const std::chrono::high_resolution_clock::time_point start_time_;
  std::chrono::high_resolution_clock::time_point elapsed_time_;
  bool stopped_ = false;
};

size_t get_memory_usage();
