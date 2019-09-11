#include <algorithm>
#include <random>

#include "mix_reads.hpp"
#include "sampling.hpp"

void SampleReads(std::istream &strand, const long unsigned &proportion, std::ostream &outfile, std::vector<long unsigned> read_ids) {
  // Sample reads from a sequence and write to out
  std::string line;
  long unsigned line_nr = 0;
  long unsigned reads_found = 0;

  while (getline(strand, line) && reads_found < proportion) {
    ++line_nr;
    if (read_ids.back() == line_nr) {
      ++reads_found;
      read_ids.pop_back();
      std::string read = line;
      for (unsigned j = 0; j < 3; ++j) {
	getline(strand, line);
	line_nr++;
	read += '\n';
	read +=line;
      }
      outfile << read << '\n';
      while (read_ids.back() == line_nr - 3) {
	outfile << read << '\n';
	read_ids.pop_back();
      }
    }
  }
}

void MixReads2(std::ifstream infiles[][2], const std::vector<double> &props, const std::vector<long unsigned> read_counts, const long unsigned &n_reads, std::ofstream &outfile_1, std::ofstream &outfile_2) {
  // Find how many reads we want from each sequences
  std::vector<long unsigned> proportions(props.size());
  long unsigned total_reads = 0;

  for (size_t i = 0; i < props.size(); ++i) {
    proportions[i] = std::ceil(props[i]*n_reads);
    total_reads += proportions[i];
  }

  // Randomly draw reads from all input samples ***WITH REPLACEMENT***
  // according to the proportions.
  std::vector<std::vector<long unsigned> > read_ids(props.size());
  DrawReadIds(proportions, read_counts, &read_ids);

  // Process the strands
  for (size_t i = 0; i < props.size(); ++i) {
    if (proportions[i] > 0) {
      const std::vector<long unsigned> &read_id = read_ids[i];
      const long unsigned &prop = proportions[i];

      std::ifstream &strand_1 = infiles[i][0];
      SampleReads(strand_1, prop, outfile_1, read_id);

      std::ifstream &strand_2 = infiles[i][1];
      SampleReads(strand_2, prop, outfile_2, read_id);
    }
  }
}

void MixReads2(igzstream infiles[][2], const std::vector<double> &props, const std::vector<long unsigned> read_counts, const long unsigned &n_reads, ogzstream &outfile_1, ogzstream &outfile_2) {
  // Find how many reads we want from each sequences
  std::vector<long unsigned> proportions(props.size());
  long unsigned total_reads = 0;

  for (size_t i = 0; i < props.size(); ++i) {
    proportions[i] = std::ceil(props[i]*n_reads);
    total_reads += proportions[i];
  }

  // Randomly draw reads from all input samples ***WITH REPLACEMENT***
  // according to the proportions.
  std::vector<std::vector<long unsigned> > read_ids(props.size());
  DrawReadIds(proportions, read_counts, &read_ids);

  // Process the strands
  for (size_t i = 0; i < props.size(); ++i) {
    if (proportions[i] > 0) {
      const std::vector<long unsigned> &read_id = read_ids[i];
      const long unsigned &prop = proportions[i];

      igzstream &strand_1 = infiles[i][0];
      SampleReads(strand_1, prop, outfile_1, read_id);

      igzstream &strand_2 = infiles[i][1];
      SampleReads(strand_2, prop, outfile_2, read_id);
    }
  }
}
