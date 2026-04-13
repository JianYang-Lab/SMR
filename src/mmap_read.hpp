#pragma once

#include <cstddef>
#include <cstdint>
#include <string>

#if defined(__linux__) || defined(__APPLE__)
#include <fcntl.h>     // open
#include <sys/mman.h>  // mmap
#include <sys/stat.h>  // fstat
#include <unistd.h>    // close
#endif

class MappedFile {
 public:
  MappedFile() : addr_(nullptr), size_(0) {}
  explicit MappedFile(std::string filename) : MappedFile() { filename_ = std::move(filename); }

  bool valid() { return addr_ != nullptr && size_ != 0; }

  void* memory_bound() { return (unsigned char*)addr_ + size_; }

  void* addr() { return addr_; }

  template <typename T>
  T offset(size_t n_bytes) {
    return reinterpret_cast<T>(reinterpret_cast<uint8_t*>(addr_) + n_bytes);
  }

  template <typename T>
  T read_from(size_t offset) {
    return *reinterpret_cast<T*>(reinterpret_cast<uint8_t*>(addr_) + offset);
  }

  bool is_end(size_t offset) const { return offset >= size_; }
  size_t remain_size(size_t offset) {
    if (is_end(offset)) return 0;
    return size_ - offset;
  }

  size_t size() { return size_; }
  const std::string& filename() const { return filename_; }

  void unmap() {
    munmap(addr_, size_);
    addr_ = nullptr;
    size_ = 0;
  }

  friend MappedFile mmap_file(const char* filename);

 private:
  void* addr_ = nullptr;
  size_t size_ = 0;
  std::string filename_;
};

MappedFile mmap_file(const char* filename);
