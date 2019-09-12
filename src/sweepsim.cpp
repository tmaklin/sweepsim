#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <algorithm>
#include <string>
#include <unistd.h>
#include <set>
#include <iterator>

#include "parse_arguments.hpp"
#include "mix_reads.hpp"
#include "sampling.hpp"
#include "util.hpp"
#include "version.h"
#include "zstr.hpp"

int main (int argc, char* argv[]) {
  std::cout << "sweepsim-" << _BUILD_VERSION << '\n' << std::endl;
  Arguments args;
  try {
    ParseArguments(argc, argv, args);
  } catch (std::runtime_error &e) {
    std::cerr << "Error in parsing arguments:\n"
	      << e.what()
	      << "\nexiting" << std::endl;
    return 0;
  } catch (std::invalid_argument &e) {
    std::cerr << e.what() << std::endl;
    PrintHelpMessage();
    return 0;
  }

  // Fixed proportions for the input files.
  short unsigned n_refs = args.infiles.size(); // Number of samples to mix from

  if (args.randomize) {
    std::cout << "\tassigning random proportions to input files" << std::endl;
    DrawRandomProportions(n_refs, &args.probs);
  }
  if (args.shuffle) {
    std::cout << "\tshuffling proportions" << std::endl;
    // Randomly assign the proportions to the input sequences
    Shuffle(args.probs);
  }
  std::cout << "\twriting assigned proportions to a file" << std::endl;
  WriteMetadata<double>(args.probs, args.infiles, args.outfile, argv);

  // Open the outfiles.
  std::cout << "Preparing the output files" << std::endl;
  std::string of1(args.outfile);
  std::string of2(args.outfile);
  of1 += "_1.fastq";
  of2 += "_2.fastq";
  
  // Prepare the input reads.
  std::cout << "Preparing the input files" << std::endl;
  if (!args.compressed) {
    std::ifstream references[n_refs][2];
    std::vector<long unsigned> read_counts(n_refs);
    for (size_t i = 0; i < n_refs; ++i) {
      std::string strand1(args.infiles[i]);
      std::string strand2(args.infiles[i]);
      strand1 += "_1.fastq";
      strand2 += "_2.fastq";
      //      read_counts[i] = CountLines(strand1);
      references[i][0].open(strand1);
      read_counts[i] = CountLines<long unsigned>(references[i][0]);
      references[i][0].close();
      references[i][0].clear();
      references[i][0].open(strand1);
      references[i][1].open(strand2);
    }
    std::ofstream outfile_1(of1);
    std::ofstream outfile_2(of2);

    std::cout << "Bootstrapping " << args.total_reads << " reads from " << n_refs << " input samples" << std::endl;
    MixReads2(references, args.probs, read_counts, args.total_reads, outfile_1, outfile_2);
  } else {
    //    igzstream references[n_refs][2];
    std::unique_ptr<std::istream> references2[n_refs][2];
    //    zstr::ifstream references[n_refs][2];
    std::vector<long unsigned> read_counts(n_refs);
    for (size_t i = 0; i < n_refs; ++i) {
      std::string strand1(args.infiles[i]);
      std::string strand2(args.infiles[i]);
      strand1 += "_1.fastq.gz";
      strand2 += "_2.fastq.gz";
      references2[i][0] = std::unique_ptr<std::istream>(new zstr::ifstream(strand1));
      read_counts[i] = CountLines<long unsigned>(*references2[i][0]);
      std::cout << read_counts[i] << std::endl;

      references2[i][0] = std::unique_ptr<std::istream>(new zstr::ifstream(strand1));
      references2[i][1] = std::unique_ptr<std::istream>(new zstr::ifstream(strand2));
    }
    zstr::ofstream outfile_1(std::string(of1 + ".gz").c_str());
    zstr::ofstream outfile_2(std::string(of2 + ".gz").c_str());

    std::cout << "Bootstrapping " << args.total_reads << " reads from " << n_refs << " input samples" << std::endl;
    MixReads2(references2, args.probs, read_counts, args.total_reads, outfile_1, outfile_2);
  }

  return(0);
}
