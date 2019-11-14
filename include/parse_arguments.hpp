#ifndef PARSE_ARGUMENTS_H
#define PARSE_ARGUMENTS_H

#include <vector>
#include <string>
#include <utility>

#include "file.hpp"

struct Arguments {
  std::vector<std::string> infile;
  std::vector<File::In> infiles[2];
  std::string outfile;
  std::vector<double> probs;
  std::pair<File::Out, File::Out> outfiles;
  File::Out meta_file;

  long unsigned total_reads;

  bool randomize = false;
  bool shuffle = false;
  bool compress = false;
  bool gzip = false;
};

void ParseArguments(int argc, char *argv[], Arguments &args);
void PrintHelpMessage();

#endif
