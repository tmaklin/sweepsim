#include <random>

#include "sampling.hpp"

std::random_device RD;
std::mt19937 RNG(RD());

void Shuffle(std::vector<double> &vec) {
  std::shuffle(vec.begin(), vec.end(), RNG);
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
