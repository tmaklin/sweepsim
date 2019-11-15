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
  File::Out log(std::cerr);
  log << "sweepsim-" << SWEEPSIM_BUILD_VERSION << '\n' << '\n';
  Arguments args;
  log << "Parsing arguments" << '\n';
  try {
    ParseArguments(argc, argv, args);
  } catch (std::invalid_argument &e) {
    log << e.what() << '\n';
    PrintHelpMessage(log);
    return 0;
  } catch (std::exception &e) {
    log << "Error in parsing arguments:\n\t"
	      << e.what()
	      << "\nexiting" << '\n';
    return 0;
  }

  // Fixed proportions for the input files.
  short unsigned n_refs = args.infile.size(); // Number of samples to mix from

  if (args.randomize) {
    log << "\tassigning random proportions to input files" << '\n';
    DrawRandomProportions(n_refs, &args.probs);
  }
  if (args.shuffle) {
    log << "\tshuffling proportions" << '\n';
    // Randomly assign the proportions to the input sequences
    Shuffle(args.probs);
  }
  log << "\twriting assigned proportions to a file" << '\n';
  WriteMetadata<double>(args.probs, args.infile, args.meta_file, argv);

  // Prepare the input reads.
  log << "Preparing the input files" << '\n';
  
  std::vector<long unsigned> read_counts(n_refs);
  for (size_t i = 0; i < n_refs; ++i) {
    read_counts[i] = CountLines<long unsigned>(args.infiles[0].at(i).stream());
    args.infiles[0].at(i).rewind();
    log << read_counts[i] << '\n';
  }
  // Open the outfiles.
  log << "Preparing the output files" << '\n';

  log << "Bootstrapping " << args.total_reads << " reads from " << n_refs << " input samples" << '\n';
  MixReads(args.infiles, args.probs, read_counts, args.total_reads, args.outfiles);

  return(0);
}
