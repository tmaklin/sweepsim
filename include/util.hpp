#ifndef UTIL_H
#define UTIL_H

#include <fstream>

#include "file.hpp"

template <typename T>
T CountLines(std::istream &fromfile) {
  T lines_count = std::count(std::istreambuf_iterator<char>(fromfile), std::istreambuf_iterator<char>(), '\n');
  lines_count /= 4;
  return(lines_count);
}

template <typename T>
void WriteMetadata(const std::vector<T> &proportions, const std::vector<std::string> &infiles, File::Out &meta, char* argv[]) {
  meta << "#ref" << '\t' << "n_reads" << '\n';
  for (size_t i = 0; i < proportions.size(); ++i) {
    meta << infiles[i] << '\t' << proportions[i] << '\n';
  }
  meta.close();
}

#endif
