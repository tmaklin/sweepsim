#ifndef MIX_READS_H
#define MIX_READS_H

#include <fstream>
#include <vector>

#include "zstr.hpp"

void MixReads(std::unique_ptr<std::istream> infiles[][2], const std::vector<double> &props, const std::vector<long unsigned> read_counts, const long unsigned &n_reads, std::unique_ptr<std::ostream> outfiles[2]);

#endif
