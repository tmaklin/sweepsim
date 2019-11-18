#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "zstr/src/zstr.hpp"
#include "cxxargs/include/cxxargs.hpp"

#include "version.h"
#include "mix_reads.hpp"
#include "sampling.hpp"
#include "util.hpp"

void parse_args(int argc, char* argv[], cxxargs::Arguments &args) {
  args.add_long_argument<bool>("help", "Print the help message", false);
  args.add_short_argument<std::string>('o', "outfile");
  args.add_short_argument<std::vector<std::string>>('f', "Prefixes for the fastq-files to mix from.");
  args.add_short_argument<uint32_t>('n', "How many reads to sample in total.");
  args.add_long_argument<bool>("random", "Draw and use random proportions for sampling. (default: false)", false);
  if (!args.value<bool>("random")) {
    args.add_long_argument<std::vector<double>>("props", "Use preset proportions rather than randomly drawn.");
  }
  args.add_long_argument<bool>("shuffle", "Shuffle the proportions before assigning them to input files. (default: false)", false);
  args.add_long_argument<bool>("gzip", "Read from input files compressed with gzip (.gz). (default: false)", false);
  args.add_long_argument<bool>("compress", "Write the output in compressed format (.gz). (default: false)", false);

  args.parse(argc, argv);
}

int main (int argc, char* argv[]) {
  File::Out log(std::cerr);
  cxxargs::Arguments args("sweepsim-" + std::string(SWEEPSIM_BUILD_VERSION), "Usage: mix_reads -o <outputFilePrefix> -f <paired-endInputPrefix> -n <numReads> [OPTIONS]");
  log << args.get_program_name() << '\n';
  log << "Parsing arguments" << '\n';
  std::vector<File::In> infiles[2];
  std::pair<File::Out, File::Out> outfiles;
  File::Out meta_file;
  try {
    parse_args(argc, argv, args);
    if (args.value<bool>("help")) {
      log << args.help() << '\n' << "exiting\n";
      log.close();
      return 0;
    }
    // Open the input and output files and throw errors if they're not accessible.
    for (size_t i = 0; i < args.value<std::vector<std::string>>('f').size(); ++i) {
      std::string name1 = args.value<std::vector<std::string>>('f').at(i) + (args.value<bool>("gzip") ? "_1.fastq.gz": "_1.fastq");
      std::string name2 = args.value<std::vector<std::string>>('f').at(i) + (args.value<bool>("gzip") ? "_2.fastq.gz": "_2.fastq");
      infiles[0].emplace_back(File::In());
      infiles[1].emplace_back(File::In());
      infiles[0].back().open(name1);
      infiles[1].back().open(name2);
    }
    outfiles.first.open(args.value<std::string>('o') + (args.value<bool>("compress") ? "_1.fastq.gz": "_1.fastq"));
    outfiles.second.open(args.value<std::string>('o')+ (args.value<bool>("compress") ? "_2.fastq.gz": "_2.fastq"));
    meta_file.open(args.value<std::string>('o') + "_info.txt");
  } catch (std::exception &e) {
    log << "Error in parsing arguments:\n\t"
	      << e.what()
	      << "\nexiting" << '\n';
    log.close();
    return 0;
  }

  // Fixed proportions for the input files.
  short unsigned n_refs = args.value<std::vector<std::string>>('f').size();

  if (args.value<bool>("random")) {
    log << "\tassigning random proportions to input files" << '\n';
    std::vector<double> probs;
    DrawRandomProportions(n_refs, &probs);
    args.set_val("props", probs);
  }
  if (args.value<bool>("shuffle")) {
    std::vector<double> probs = args.value<std::vector<double>>("props");
    log << "\tshuffling proportions" << '\n';
    // Randomly assign the proportions to the input sequences
    Shuffle(probs);
    args.set_val("props", probs);
  }
  log << "\twriting assigned proportions to a file" << '\n';
  WriteMetadata<double>(args.value<std::vector<double>>("props"),
  			args.value<std::vector<std::string>>('f'),
  			meta_file, argv);

  // Prepare the input reads.
  log << "Preparing the input files" << '\n';
  
  std::vector<long unsigned> read_counts(n_refs);
  for (size_t i = 0; i < n_refs; ++i) {
    read_counts[i] = CountLines<long unsigned>(infiles[0].at(i).stream());
    infiles[0].at(i).rewind();
  }
  // Open the outfiles.
  log << "Preparing the output files" << '\n';

  log << "Bootstrapping " << args.value<uint32_t>('n') << " reads from " << n_refs << " input samples" << '\n';
  MixReads(infiles, args.value<std::vector<double>>("props"), read_counts, args.value<uint32_t>('n'), outfiles);

  return(0);
}
