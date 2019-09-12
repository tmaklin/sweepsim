# sweepsim
Simulate mixtures of Illumina sequencing reads (paired-end) by
bootstrapping from existing samples, either according to prespecified
or randomized proportions.

# Installation
## Requirements
- C++11 compliant compiler.
- cmake

- Clone repository, enter the directory and run
```
> mkdir build
> cd build
> cmake ..
> make
```
- You can remove the build directory afterwards.

# Input
Input paired-end sample(s) must be in separate files and in the
following format:
```
sample-1_1.fastq
sample-1_2.fastq
```
i.e. the strands are in separate files. The files may be either
uncompressed, or compressed with gzip (.gz format)

# Usage
Bootstrap 100 reads from a single sequencing sample with random proportions
> sweepsim -f sample-1 -o mixture-reads -n 100

... 100 reads from two sequencing sample-s
> sweepsim -f sample-1,sample-2 -o mixture-reads -n 100

... with preset proportions (0.20, 0.80)
> sweepsim -f sample-1,sample-2 --probs 0.20,0.80 -o mixture-reads -n 100

... with preset proportions, but assigned to random inputs
> sweepsim -f sample-1,sample-2 --probs 0.20,0.80 -o mixture-reads -n 100 --shuffle

... from gzipped (.gz) input files
> sweepsim -f sample-1,sample-2 --probs 0.20,0.80 -o mixture-reads -n 100 --shuffle --gzip

... into gzipped (.gz) output files
> sweepsim -f sample-1,sample-2 --probs 0.20,0.80 -o mixture-reads -n 100 --shuffle --gzip --compress

# Output
The output will contain three files:
```
mixture-reads_1.fastq // Reads sampled from the 1st strand
mixture-reads_2.fastq // Reads sampled from the 2nd strand
mixture-reads_info.txt // Info about the input files and their sampling proportions
```

# Help
sweepsim recognizes the following flags:

```
-o <outputFilePrefix>
Output file prefix.
-f <paired-endInputPrefix>
prefixes for the fastq-files to mix from.
-n <numReads>
how many reads to sample in total.
--props <proportions>
Use preset proportions rather (default: randomly drawn from uniform distribution and normalize).
--random <randomizeProportions>
draw and use random proportions for sampling.
--shuffle <shuffleProportions>
shuffle the proportions before assignign them to input files.
--gzip <gzippedInput>
Read from input files compressed with gzip (.gz).
--compress <compressOutput>
Write the output in compressed format (.gz).
--help
print help message.
```
