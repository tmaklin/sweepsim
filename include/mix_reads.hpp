#ifndef MIX_READS_H
#define MIX_READS_H

#include <fstream>
#include <vector>
#include <memory>

#include "cxxio.hpp"

void MixReads(std::vector<cxxio::In> infiles[2], const std::vector<double> &props, const std::vector<long unsigned> read_counts, const long unsigned &n_reads, std::pair<cxxio::Out, cxxio::Out> &outfiles);

#endif
