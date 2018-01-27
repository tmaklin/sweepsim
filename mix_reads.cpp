#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <algorithm>
#include <string>
#include <unistd.h>
#include <set>
#include <iterator>

std::random_device RD;
std::mt19937 RNG(RD());

void SampleReads(std::ifstream &strand, const long unsigned &proportion, std::ofstream &outfile, std::vector<long unsigned> read_ids) {
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

void DrawReadIds(const std::vector<long unsigned> &how_many, std::vector<long unsigned> read_counts, std::vector<std::vector<long unsigned> > *random_integers) {
  // Draw random ints and return them in reverse order
  for (size_t i = 0; i < read_counts.size(); ++i) {
    std::uniform_int_distribution<int> uni(1, read_counts[i]);
    std::vector<long unsigned> &random = (*random_integers)[i];
    for (unsigned j = 0; j < how_many[i]; ++j) {
      long unsigned nextline = (uni(RNG) - 1)*4 + 1;
      random.push_back(nextline);
    }
    std::sort(random.rbegin(), random.rend());
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

long unsigned CountLines(const std::string &infile) {
  std::ifstream fromfile(infile);
  long unsigned lines_count = std::count(std::istreambuf_iterator<char>(fromfile), std::istreambuf_iterator<char>(), '\n');
  lines_count /= 4;
  return(lines_count);
}

void WriteMetadata(const std::vector<double> &proportions, std::string &filename, char* argv[]) {
  filename += "_info.txt";
  std::ofstream meta(filename);
  meta << "#ref" << '\t' << "n_reads" << std::endl;
  for (size_t i = 0; i < proportions.size(); ++i) {
    meta << argv[i + 1] << '\t' << proportions[i] << std::endl;
  }
  meta.close();
}

void DrawRandomProportions(const int &how_many, std::vector<double> *random_reals) {
  // Draw random proportions
  std::uniform_real_distribution<double> uni(0, 1);
  double total = 0.0;
  for (int i = 0; i < how_many; ++i) {
    double real = uni(RNG);
    total += real;
    random_reals->push_back(real);
  }
  for (int i = 0; i < how_many; ++i) {
    (*random_reals)[i] /= total;
  }
}

int main (int argc, char* argv[]) {
  // Fixed proportions for the input files.
  unsigned n_refs = argc - 3; // Number of samples to mix

  std::cout << "Bootstrapping reads from " << n_refs << (n_refs > 1 ? " samples." : " sample.") << std::endl;
  std::vector<double> props;
  DrawRandomProportions(n_refs, &props);
  
  // Randomly assign the proportions to the input sequences
  std::shuffle(props.begin(), props.end(), RNG);
  
  // Write the proportions to a file
  std::string filename(argv[n_refs + 1]);
  WriteMetadata(props, filename, argv);

  // Prepare the input reads.
  std::ifstream references[n_refs][2];
  std::vector<long unsigned> read_counts(n_refs);
  for (unsigned i = 0; i < n_refs; ++i) {
    std::string strand1(argv[i + 1]);
    std::string strand2(argv[i + 1]);
    strand1 += "_1.fastq";
    strand2 += "_2.fastq";
    read_counts[i] = CountLines(strand1);
    references[i][0].open(strand1);
    references[i][1].open(strand2);
  }

  // Open the outfiles.
  std::string of1(argv[n_refs + 1]);
  std::string of2(argv[n_refs + 1]);
  of1 += "_1.fastq";
  of2 += "_2.fastq";
  std::ofstream outfile_1(of1);
  std::ofstream outfile_2(of2);
  long unsigned total_reads = std::atoi(argv[argc - 1]);
  MixReads2(references, props, read_counts, total_reads, outfile_1, outfile_2);
  return(0);
}
