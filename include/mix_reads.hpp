#ifndef MIX_READS_H
#define MIX_READS_H

#include <fstream>
#include <vector>

#include "zstr.hpp"

void MixReads2(std::ifstream infiles[][2], const std::vector<double> &props, const std::vector<long unsigned> read_counts, const long unsigned &n_reads, std::ofstream &outfile_1, std::ofstream &outfile_2);

void MixReads2(zstr::ifstream infiles[][2], const std::vector<double> &props, const std::vector<long unsigned> read_counts, const long unsigned &n_reads, zstr::ofstream &outfile_1, zstr::ofstream &outfile_2);

void MixReads2(std::unique_ptr<std::istream> infiles[][2], const std::vector<double> &props, const std::vector<long unsigned> read_counts, const long unsigned &n_reads, zstr::ofstream &outfile_1, zstr::ofstream &outfile_2);

#endif
