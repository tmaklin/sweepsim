#include <algorithm>
#include <iostream>
#include <exception>
#include <sstream>
#include "parse_arguments.hpp"

void PrintHelpMessage() {
  std::cerr << "Usage: mix_reads -o <outputFilePrefix> -f <paired-endInputPrefix> -n <numReads> [OPTIONS]\n"
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
	    << '\n'
	    << "\t--help"
	    << "\tprint this message." << std::endl;
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
  std::cerr << "Parsing arguments" << std::endl;

  args.randomize = CmdOptionPresent(argv, argv+argc, "--random");
  args.shuffle = CmdOptionPresent(argv, argv+argc, "--shuffle");

  if (CmdOptionPresent(argv, argv+argc, "-o")) {
    args.outfile = std::string(GetCmdOption(argv, argv+argc, "-o"));
  } else {
    throw std::runtime_error("outfile must be specified");
  }

  if (CmdOptionPresent(argv, argv+argc, "-f")) {
    std::string infiles(GetCmdOption(argv, argv+argc, "-f"));
    SplitArguments(infiles, &args.infiles);
  } else {
    throw std::runtime_error("no input reads");
  }

  if (CmdOptionPresent(argv, argv+argc, "-n")) {
    args.total_reads = std::stol(std::string(GetCmdOption(argv, argv+argc, "-n")));
  } else {
    throw std::runtime_error("total number of reads to sample is missing");
  }

  if (CmdOptionPresent(argv, argv+argc, "--probs") & !args.randomize) {
    std::string probs(GetCmdOption(argv, argv+argc, "--probs"));
    std::vector<std::string> probs_str;
    SplitArguments(probs, &probs_str);
    args.probs = ValidateProbs(probs_str);
    if (args.probs.size() != args.infiles.size()) {
      throw std::runtime_error("number of input files does not match number of input probabilities");
    }
  } else {
    args.randomize = true;
  }
}
