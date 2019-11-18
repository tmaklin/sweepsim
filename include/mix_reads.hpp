#ifndef MIX_READS_H
#define MIX_READS_H

#include <fstream>
#include <vector>

#include "file.hpp"

void MixReads(std::vector<File::In> infiles[2], const std::vector<double> &props, const std::vector<long unsigned> read_counts, const long unsigned &n_reads, std::pair<File::Out, File::Out> &outfiles);

#endif
