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

  // Prepare the input reads.
  std::cout << "Preparing the input files" << std::endl;
  
  std::unique_ptr<std::istream> references[n_refs][2];
  std::vector<long unsigned> read_counts(n_refs);
  for (size_t i = 0; i < n_refs; ++i) {
    std::string strand1(args.infiles[i]);
    std::string strand2(args.infiles[i]);
    strand1 += (args.gzip ? "_1.fastq.gz": "_1.fastq");
    strand2 += (args.gzip ? "_2.fastq.gz": "_2.fastq");
    references[i][0] = std::unique_ptr<std::istream>(new zstr::ifstream(strand1));
    read_counts[i] = CountLines<long unsigned>(*references[i][0]);
    std::cout << read_counts[i] << std::endl;

    references[i][0] = std::unique_ptr<std::istream>(new zstr::ifstream(strand1));
    references[i][1] = std::unique_ptr<std::istream>(new zstr::ifstream(strand2));
  }

  std::unique_ptr<std::ostream> outfiles[2];
  // Open the outfiles.
  std::cout << "Preparing the output files" << std::endl;
  std::string of1(args.outfile);
  std::string of2(args.outfile);
  of1 += (args.compress ? "_1.fastq.gz": "_1.fastq");
  of2 += (args.compress ? "_2.fastq.gz": "_2.fastq");
  
  if (args.compress) {
    outfiles[0] = std::unique_ptr<std::ostream>(new zstr::ofstream(of1));
    outfiles[1] = std::unique_ptr<std::ostream>(new zstr::ofstream(of2));
  } else {
    outfiles[0] = std::unique_ptr<std::ostream>(new std::ofstream(of1));
    outfiles[1] = std::unique_ptr<std::ostream>(new std::ofstream(of2));
  }

  std::cout << "Bootstrapping " << args.total_reads << " reads from " << n_refs << " input samples" << std::endl;
  MixReads(references, args.probs, read_counts, args.total_reads, outfiles);

  return(0);
}
