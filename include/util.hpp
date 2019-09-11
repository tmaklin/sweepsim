#ifndef UTIL_H
#define UTIL_H

#include <fstream>

template <typename T>
T CountLines(std::istream &fromfile) {
  T lines_count = std::count(std::istreambuf_iterator<char>(fromfile), std::istreambuf_iterator<char>(), '\n');
  lines_count /= 4;
  return(lines_count);
}

template <typename T>
void WriteMetadata(const std::vector<T> &proportions, const std::vector<std::string> &infiles, std::string filename, char* argv[]) {
  filename += "_info.txt";
  std::ofstream meta(filename);
  meta << "#ref" << '\t' << "n_reads" << std::endl;
  for (size_t i = 0; i < proportions.size(); ++i) {
    meta << infiles[i] << '\t' << proportions[i] << std::endl;
  }
  meta.close();
}

#endif
