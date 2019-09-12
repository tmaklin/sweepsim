#ifndef PARSE_ARGUMENTS_H
#define PARSE_ARGUMENTS_H

#include <vector>
#include <string>

struct Arguments {
  std::vector<std::string> infiles;
  std::string outfile;
  std::vector<double> probs;

  long unsigned total_reads;

  bool randomize = false;
  bool shuffle = false;
  bool compressed = false;
  bool gzip = false;
};

void ParseArguments(int argc, char *argv[], Arguments &args);
void PrintHelpMessage();

#endif
