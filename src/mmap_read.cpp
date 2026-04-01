#include "mmap_read.hpp"

#include <cerrno>
#include <cstdio>
#include <cstring>

MappedFile mmap_file(const char* filename) {
  std::string filename_str = filename;
  MappedFile mapped_file(filename_str);
  int fd = open(filename, O_RDONLY);
  if (fd == -1) {
    fprintf(stderr, "open file failed, path=%s", filename);
    return mapped_file;
  }

  struct stat st;
  if (fstat(fd, &st) == -1) {
    close(fd);
    fprintf(stderr, "failed to get file size, path=%s, err=%s", filename, strerror(errno));
    return mapped_file;
  }
  size_t filesize = st.st_size;

  if (filesize % sizeof(float) != 0) {
    fputs("file size is not a multiple of sizeof(float)!", stderr);
    return mapped_file;
  }

  void* mapped_addr = mmap(nullptr, filesize, PROT_READ, MAP_PRIVATE, fd, 0);

  close(fd);

  if (mapped_addr == MAP_FAILED) {
    fprintf(stderr, "failed to mmap file, path=%s", filename);
    return mapped_file;
  }

  madvise(mapped_addr, filesize, MADV_SEQUENTIAL | MADV_WILLNEED);

  mapped_file.addr_ = mapped_addr;
  mapped_file.size_ = filesize;

  return mapped_file;
}
