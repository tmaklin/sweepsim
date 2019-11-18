#include "parse_arguments.hpp"

#include <algorithm>
#include <exception>
#include <sstream>

#include "file.hpp"

void PrintHelpMessage(File::Out &out) {
  out << "Usage: mix_reads -o <outputFilePrefix> -f <paired-endInputPrefix> -n <numReads> [OPTIONS]\n"
	    << "Mixes sequencing reads from multiple files together; randomly or according to some proportions..\n\n"
	    << "Options:\n"
    	    << "\t-o <outputFilePrefix>\n"
	    << "\tOutput file prefix.\n"
	    << "\t-f <paired-endInputPrefix>\n"
	    << "\tprefixes for the fastq-files to mix from.\n"
	    << "\t-n <numReads>\n"
	    << "\thow many reads to sample in total.\n"
    	    << "\t--props <proportions>\n"
	    << "\tUse preset proportions rather (default: randomly drawn from uniform distribution and normalize).\n"
	    << "\t--random <randomizeProportions>\n"
	    << "\tdraw and use random proportions for sampling.\n"
	    << "\t--shuffle <shuffleProportions>\n"
	    << "\tshuffle the proportions before assignign them to input files.\n"
	    << "\t--gzip <gzippedInput>\n"
	    << "\tRead from input files compressed with gzip (.gz).\n"
	    << "\t--compress <compressOutput>\n"
	    << "\tWrite the output in compressed format (.gz).\n"
	    << '\n'
	    << "\t--help"
      << "\tprint this message." << '\n';
}

char* GetCmdOption(char **begin, char **end, const std::string &option) {
  char **it = std::find(begin, end, option);
  return ((it != end && ++it != end) ? *it : 0);
}

bool CmdOptionPresent(char **begin, char **end, const std::string &option) {
  return (std::find(begin, end, option) != end);
}

double ParseDoubleOption(char **begin, char **end, const std::string &option) {
  double opt;
  char* char_opt = GetCmdOption(begin, end, option);
  if (char_opt != 0) {
    opt = std::stod(std::string(char_opt));
  } else {
    throw std::runtime_error(option + " specified but no value given");
  }
  return(opt);
}

void SplitArguments(const std::string &input, std::vector<std::string> *output) {
  std::string object;
  for (std::stringstream sst(input); getline(sst, object, ',');) {
    output->emplace_back(object);
  }
}

std::vector<double> ValidateProbs(const std::vector<std::string> &probs_str) {
  double sum = 0.0;
  std::vector<double> probs;
  for (auto prob_str : probs_str) {
    double prob = std::stod(prob_str);
    sum += prob;
    if (prob < 0 | prob > 1) {
      throw std::runtime_error("given probability: " + std::to_string(prob) + " is not between 0 and 1");
    }
    probs.push_back(prob);
  }
  if (sum < 1.0 - 1e-16 | sum > 1.0 + 1e-16) {
    throw std::runtime_error("probabilities do not sum to 1");
  }
  return probs;
}

void ParseArguments(int argc, char *argv[], Arguments &args) {
  if (CmdOptionPresent(argv, argv+argc, "--help")) {
    throw std::invalid_argument("");
  }

  args.randomize = CmdOptionPresent(argv, argv+argc, "--random");
  args.shuffle = CmdOptionPresent(argv, argv+argc, "--shuffle");
  args.compress = CmdOptionPresent(argv, argv+argc, "--compress");
  args.gzip = CmdOptionPresent(argv, argv+argc, "--gzip");

  if (CmdOptionPresent(argv, argv+argc, "-o")) {
    args.outfile = std::string(GetCmdOption(argv, argv+argc, "-o"));
    args.outfiles.first.open(args.outfile + (args.compress ? "_1.fastq.gz": "_1.fastq"));
    args.outfiles.second.open(args.outfile + (args.compress ? "_2.fastq.gz": "_2.fastq"));
    args.meta_file.open(args.outfile + "_info.txt");
  } else {
    throw std::runtime_error("outfile must be specified");
  }

  if (CmdOptionPresent(argv, argv+argc, "-f")) {
    std::string infile(GetCmdOption(argv, argv+argc, "-f"));
    SplitArguments(infile, &args.infile);
    for (size_t i = 0; i < args.infile.size(); ++i) {
      std::string name1 = args.infile.at(i) + (args.gzip ? "_1.fastq.gz": "_1.fastq");
      std::string name2 = args.infile.at(i) + (args.gzip ? "_2.fastq.gz": "_2.fastq");
      args.infiles[0].emplace_back(File::In());
      args.infiles[1].emplace_back(File::In());
      args.infiles[0].back().open(name1);
      args.infiles[1].back().open(name2);
    }
  } else {
    throw std::runtime_error("no input reads");
  }

  if (CmdOptionPresent(argv, argv+argc, "-n")) {
    args.total_reads = std::stol(std::string(GetCmdOption(argv, argv+argc, "-n")));
  } else {
    throw std::runtime_error("total number of reads to sample is missing");
  }

  if (CmdOptionPresent(argv, argv+argc, "--props") & !args.randomize) {
    std::string probs(GetCmdOption(argv, argv+argc, "--props"));
    std::vector<std::string> probs_str;
    SplitArguments(probs, &probs_str);
    args.probs = ValidateProbs(probs_str);
    if (args.probs.size() != args.infiles[0].size()) {
      throw std::runtime_error("number of input files does not match number of input probabilities");
    }
  } else {
    args.randomize = true;
  }
}
