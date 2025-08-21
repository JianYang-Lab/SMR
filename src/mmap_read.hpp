#pragma once

#include <cstddef>
#include <string>

#if defined(__linux__)
#include <fcntl.h>    // open
#include <sys/mman.h> // mmap
#include <sys/stat.h> // fstat
#include <unistd.h>   // close
#endif

class MappedFile {
public:
  MappedFile(): addr_(nullptr), size_(0) {}
  explicit MappedFile(std::string filename) : MappedFile() {
      filename_ = std::move(filename);
  }

  bool valid() { return addr_ != nullptr && size_ != 0; }

  void *memory_bound() { return (unsigned char *)addr_ + size_; }

  void *addr() { return addr_; }
  size_t size() { return size_; }
  const std::string& filename() const { return filename_; }

  void unmap() {
    munmap(addr_, size_);
    addr_ = nullptr;
    size_ = 0;
  }

  friend MappedFile mmap_file(const char* filename);

private:
  void *addr_ = nullptr;
  size_t size_ = 0;
  std::string filename_;
};


MappedFile mmap_file(const char* filename);
