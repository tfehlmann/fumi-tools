# fumi-tools

This tool is intended to deduplicate UMIs from single-end sequencing data. Reads are considered identical when their UMIs have the same sequence, they have the same length and map at the same position.

## Installation

This code was tested on Ubuntu 16.04, 17.10, 18.04 and Arch Linux with GCC 5, 6, 7 and 9.

### Pre-built binaries

Download the latest binaries from the [release page](https://ccb-gitlab.cs.uni-saarland.de/tobias/fumi_tools/releases). These are self-contained and should work on any 64-bit Linux system.

### Conda

Binaries are also available on our conda channel and can be installed with the following command:

```bash
conda install -c ccb-sb fumi_tools
```

### Build the code

Download the source code either from the [release page](https://ccb-gitlab.cs.uni-saarland.de/tobias/fumi_tools/releases) (Source code /w dependencies) or clone the repository with [git](https://git-scm.com/). Then build the code with [CMake](https://cmake.org/) (if installation_path is omitted, will default to /usr/local). The only dependencies for building the code are a GCC compiler (might work with Clang) and CMake.

```bash
cd {path_to_source_code}
# e.g. cd ~/Downloads/fumi_tools
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX={installation_path} ..
# e.g. to install the binaries in ~/bin
# cmake -DCMAKE_INSTALL_PREFIX=~/ ..
cmake --build . --target install -- -j {number_of_threads_for_compilation}
# e.g. cmake --build . --target install -- -j 4

```



## Usage

```bash
usage: fumi_tools <command> [<args>]
Version: <version>

Available commands are:
    copy_umi  Copy UMI from FASTQ files into their header
    dedup     Deduplicate reads in BAM files
```

### First step - copy UMI into read header

First copy the UMI's into the read name (the sequences remain unchanged):

```bash
usage: fumi_tools copy_umi [-h] -i INPUT -o OUTPUT --umi-length UMI_LENGTH
                           [--threads THREADS] [--version]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input FASTQ file, optionally gzip compressed.
  -o OUTPUT, --output OUTPUT
                        Output FASTQ file, optionally gzip compressed.
  --umi-length UMI_LENGTH
                        Length of the UMI to copy. It is assumed that the UMI
                        starts at the 5' end of the read. (default: None)
  --threads THREADS     Number of threads to use. (default: 1)
  --version             Display version number
```



```bash
# e.g. for dummy.fastq.gz with UMI of length 12
fumi_tools copy_umi --input dummy.fastq.gz --umi-length {umi_length} --output dummy.umi.fastq.gz
```

### Second step - deduplicate alignment file

After you aligned the reads (e.g. with STAR) deduplicate them:

```bash
usage: fumi_tools dedup [-h] -i INPUT -o OUTPUT [--start-only]
                        [--threads THREADS] [--memory MEMORY] [--seed SEED]
                        [--version]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input SAM or BAM file. Needs to be coordinate sorted.
  -o OUTPUT, --output OUTPUT
                        Output SAM or BAM file, sorted by read name. To output
                        SAM on stdout use '-'.
  --start-only          Reads only need the same start position and the same
                        UMI to be considered duplicates. (default: False)
  --threads THREADS     Number of threads to use. (default: 1)
  --memory MEMORY       Maximum memory used for sorting. Units can be K/M/G.
                        (default: 3G)
  --seed SEED           Random number generator seed. (default: 42)
  --version             Display version number.
```



```bash
# e.g. for dummy_aligned.bam using 4 threads and 3 gigabytes of RAM
fumi_tools dedup -i dummy_aligned.bam -o dummy_aligned.dedup.bam --threads 4 --memory 3G
```