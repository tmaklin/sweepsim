#ifndef SAMPLING_H
#define SAMPLING_H

#include <vector>

void Shuffle(std::vector<double> &vec);
void DrawReadIds(const std::vector<long unsigned> &how_many, std::vector<long unsigned> read_counts, std::vector<std::vector<long unsigned> > *random_integers);
void DrawRandomProportions(const int &how_many, std::vector<double> *random_reals);

#endif
