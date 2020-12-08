[![pipeline status](https://ccb-gitlab.cs.uni-saarland.de/tobias/fumi_tools/badges/master/pipeline.svg)](https://ccb-gitlab.cs.uni-saarland.de/tobias/fumi_tools/commits/master) [![Anaconda-Server Badge](https://anaconda.org/ccb-sb/fumi_tools/badges/platforms.svg)](https://anaconda.org/ccb-sb/fumi_tools) [![Anaconda-Server Badge](https://anaconda.org/ccb-sb/fumi_tools/badges/version.svg)](https://anaconda.org/ccb-sb/fumi_tools) [![Commitizen friendly](https://img.shields.io/badge/commitizen-friendly-brightgreen.svg)](http://commitizen.github.io/cz-cli/)


# fumi-tools

This tool is intended to deduplicate UMIs from single-end and paired-end sequencing data. Reads are considered identical when their UMIs have the same sequence, they have the same length and map at the same position.

## Installation

This code was tested on Ubuntu 16.04, 17.10, 18.04 and Arch Linux with GCC 5, 6, 7, 9 and 10.

### Pre-built binaries

Download the latest binaries from the [release page](https://ccb-gitlab.cs.uni-saarland.de/tobias/fumi_tools/releases). These are self-contained and should work on any 64-bit Linux system.

### Conda [![Anaconda-Server Badge](https://anaconda.org/ccb-sb/fumi_tools/badges/installer/conda.svg)](https://anaconda.org/ccb-sb/fumi_tools)

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

### Optional dependencies

We optionally use [pigz](https://github.com/madler/pigz) to compress FASTQ files with multiple threads. This tool can usually be installed over the distribution package manager or as alternative over [conda](https://anaconda.org/anaconda/pigz).

#### Ubuntu

```bash
sudo apt install pigz
```

#### Arch Linux

```bash
sudo pacman -S pigz
```

#### Conda

```bash
conda install pigz
```



## Usage

```bash
usage: fumi_tools <command> [<args>]
Version: <version>

Available commands are:
    demultiplex  Demultiplex single FASTQ file
    copy_umi     Copy UMI from FASTQ files into their header
    dedup        Deduplicate reads in BAM files
```

### First step - demultiplexing and/or copying UMI into read header

In case your sequences need to be demultiplexed:

```bash
usage: fumi_tools demultiplex [-h] -i INPUT [-I INPUT_READ2] -s SAMPLE_SHEET -o OUTPUT [-e MAX_ERRORS] [-l LANE [LANE ...]] [--format-umi] [--tag-umi] [--threads THREADS] [--version]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input FASTQ file, optionally gzip compressed.
  -I INPUT_READ2, --input-read2 INPUT_READ2
                        Input paired end R2 FASTQ file, optionally gzip compressed.
  -s SAMPLE_SHEET, --sample-sheet SAMPLE_SHEET
                        Sample Sheet in Illumina format. (default: None)
  -o OUTPUT, --output OUTPUT
                        Output FASTQ file pattern, optionally gzip compressed. Use %i as placeholder for the sample index specified in the sample sheet, %s for the sample name, %l for the lane and optionally %r for the read direction (e.g. demultiplexed_reads/%s_S%i_L%l_R%r.fastq.gz).
  -e MAX_ERRORS, --max-errors MAX_ERRORS
                        Maximum allowed number of errors (mismatches per default). (default: 1)
  -l LANE [LANE ...], --lane LANE [LANE ...]
                        Optionally specify on which lane the samples provided in the sample sheet ran. Can be specified multiple times to pass several lanes. This option takes precedence on the Lane column of the sample sheet. (default: None)
  --format-umi          Add UMI to the end of the FASTQ ID (before the first space in the header), as expected by fumi_tools dedup (default: False)
  --tag-umi             Add UMI to the read ID by adding :FUMI|<UMI_SEQ>| instead of a simple underscore. (default: False)
  --threads THREADS     Number of threads to use. (default: 1)
  --version             Display version number.
```

```bash
# e.g. for dummy_R1.fastq.gz containing multiple samples
fumi_tools demultiplex --input dummy_R1.fastq.gz --sample-sheet sample_sheet.csv --output output_folder/%s_S%i_L%l_R1.fastq.gz
```

The program expects the read header to be formatted as follows (which corresponds to the output of bcl2fastq 2):

```bash
@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<I7-index><UMI>+<I5-index>
```

Below is an example of a minimal sample sheet. The index column contains the i7 index and the index2 column contains the i5 index. Optionally a lane column can be specified.

```bash
Sample_ID,Sample_Name,index,index2
1,Sample234,CGCATGAT,TCAGGCTT
2,Sample745,ACGGAACA,GTTCTCGT
```

All reads not matching any valid index combination will be outputted in a file with Sample_ID 0 and Sample_Name Undetermined.

**ATTENTION: if you are using downstream tools that change the FASTQ header by inserting additional elements using underscores (e.g. bismark), or if you are unsure about it, use the --tag-umi option to copy the UMI into the read header. This will become the default option in the near future.**

Alternatively, if your reads have already been demultiplexed and the UMI sequence is present in the read sequence, copy the UMI into the read header (the sequences remain unchanged):

```bash
usage: fumi_tools copy_umi [-h] -i INPUT [-I INPUT_READ2] -o OUTPUT [-O OUTPUT_READ2] --umi-length UMI_LENGTH [--tag-umi] [--threads THREADS] [--version]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input FASTQ file, optionally gzip, bz2 or xz compressed.
  -I INPUT_READ2, --input-read2 INPUT_READ2
                        Input paired end R2 FASTQ file, optionally gzip, bz2 or xz compressed.
  -o OUTPUT, --output OUTPUT
                        Output FASTQ file, optionally gzip, bz2 or xz compressed.
  -O OUTPUT_READ2, --output-read2 OUTPUT_READ2
                        Output paired end R2 FASTQ file, optionally gzip, bz2 or xz compressed.
  --umi-length UMI_LENGTH
                        Length of the UMI to copy. It is assumed that the UMI starts at the 5\' end of the read. (default: None)
  --tag-umi             Add UMI to the read ID by adding :FUMI|<UMI_SEQ>| instead of a simple underscore. (default: False)
  --threads THREADS     Number of threads to use. (default: 1)
  --version             Display version number.
```

```bash
# e.g. for dummy.fastq.gz with UMI of length 12
fumi_tools copy_umi --input dummy.fastq.gz --umi-length {umi_length} --output dummy.umi.fastq.gz
```

### Second step - deduplicate alignment file

After you aligned the reads (e.g. with STAR) deduplicate them:

```bash
usage: fumi_tools dedup [-h] -i INPUT -o OUTPUT [--paired] [--start-only]
                        [--threads THREADS] [--memory MEMORY]
                        [--seed SEED] [--version]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input SAM or BAM file. Needs to be coordinate sorted.
  -o OUTPUT, --output OUTPUT
                        Output SAM or BAM file, sorted by read name. To output
                        SAM on stdout use '-'.
  --paired              Specify this option if your input contains paired end reads.
  --start-only          Reads only need the same start position and the same
                        UMI to be considered duplicates.
  --chimeric-pairs [{discard,use}]
                        How to handle chimeric read pairs. (default: use)
  --unpaired-reads [{discard,use}]
                        How to handle unpaired reads (e.g. mate did not align) (default: use)
  --sort-adjacent-pairs
                        Keep name sorting, but sort pairs such that the mate always follows the first read.
  --threads THREADS     Number of threads to use. (default: 1)
  --memory MEMORY       Maximum memory used for sorting. Units can be K/M/G. (default: 3G)
  --seed SEED           Random number generator seed. (default: 42)
  --version             Display version number.
```



```bash
# e.g. for dummy_aligned.bam using 4 threads and 3 gigabytes of RAM
fumi_tools dedup -i dummy_aligned.bam -o dummy_aligned.dedup.bam --threads 4 --memory 3G
```
