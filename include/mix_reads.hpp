#ifndef MIX_READS_H
#define MIX_READS_H

#include <fstream>
#include <vector>

#include "gzstream.h"

void MixReads2(std::ifstream infiles[][2], const std::vector<double> &props, const std::vector<long unsigned> read_counts, const long unsigned &n_reads, std::ofstream &outfile_1, std::ofstream &outfile_2);

void MixReads2(igzstream infiles[][2], const std::vector<double> &props, const std::vector<long unsigned> read_counts, const long unsigned &n_reads, ogzstream &outfile_1, ogzstream &outfile_2);

#endif
