# MetaCache

[![Linux build status](https://travis-ci.org/muellan/metacache.svg?branch=master)](https://travis-ci.org/muellan/metacache)

MetaCache is a classification system for mapping genomic sequences (short reads, long reads, contigs, ...) from metagenomic samples to their most likely taxon of origin. MetaCache aims to reduce the memory requirement usually associated with k-mer based methods while retaining their speed. It uses locality sensitive hashing to quickly identify candidate regions within one or multiple reference genomes. A read is then classified based on the similarity to those regions.

For an independend comparison to other tools in terms of classification accuracy see the [LEMMI](https://lemmi.ezlab.org) benchmarking site.

**MetaCache's CPU** version classifies around 60 Million reads (of length 100) per minute against all complete bacterial, viral and archaea genomes from NCBI RefSeq Release 97 running with 88 threads on a workstation with 2 Intel(R) Xeon(R) Gold 6238 CPUs.

**MetaCache's [GPU version](docs/gpu_version.md)** classifies around 300 Million reads (of length 100) per minute against all complete bacterial, viral, fungal and archaea genomes from NCBI RefSeq Release 202 running on a workstation with 4 NVIDIA(R) Tesla(R) V100 GPUs (32 GB model).




## Quick Start with NCBI RefSeq
on a Debian/Ubuntu system:

```
sudo apt install -y zlib1g zlib1g-dev
git clone https://github.com/muellan/metacache.git
cd metacache
make
./metacache-build-refseq
```
This will
  * install the zlib library  
  * download the MetaCache source code from GitHub
  * compile MetaCache (without GPU support)
  * download the complete bacterial, viral and archaea genomes from the latest NCBI RefSeq release (this can take some time) 
  * build a classification database

Once the default database is built you can classify reads:
  ```
  ./metacache query refseq myReads.fa -out results.txt
  ./metacache query refseq anyFolderWithFastaOrFastqFiles -out results.txt
  ./metacache query refseq myReads1.fa myReads2.fa -pairfiles -out results.txt
  ./metacache query refseq myPairedReads.fa -pairseq -out results.txt
  ```




## Detailed Installation Instructions

Visit MetaCache's github [repository] to get the latest resources.

* To compile the CPU version: run `make` in the directory containing the Makefile
* To compile the GPU version, follow the instructions provided [here](docs/gpu_version.md).


### CPU Version Requirements

MetaCache itself should compile on any platform for which a C++14 conforming compiler is available. The Makefile is written with g++ or clang++ in mind, but could probably be adapted to (a very recent version of) MSVC or other compilers.

The helper scripts (for downloading genomes, taxonomy etc.) require the Bash shell to run. That means you need a working bash executable as well as some common GNU utilities like "awk" and "wget". On Windows you should use the 'Windows Subsystem for Linux' (which gives you an Ubuntu user mode talking to the Windows Kernel).

MetaCache 2.0.0 was successfully tested on the following platforms (all 64 bit + 64 bit compilers):
- Ubuntu 20.04 with g++ 5.4, g++ 7.4
- Windows 10 20H2 running Ubuntu 20.04 inside WSL2 and g++ 10.3

In order to be able to build the default database (based on NCBI RefSeq Release 97) with default settings your system should have around 64GB of RAM (note that the NCBI RefSeq will still be growing in the near future).
If you don't have enough RAM, you can use [database partitioning](docs/partitioning.md).


### GPU Version Requirements
The GPU version requires a CUDA-capable device of the Pascal generation or newer and either CUDA >= 11 or CUDA 10.2 and a self-provided version of [CUB](https://github.com/NVlabs/cub).

See [here](docs/gpu_version.md) for more.


### Library Requirements (CPU & GPU versions)
MetaCache requires the zlib compression library to be installed on your system in order to be able to process gzipped FASTA/FASTQ files.
On Debian/Ubuntu zlib can be installed with
  ```
  sudo apt install -y zlib1g zlib1g-dev
  ```
If you *don't* have zlib installed or cannot do so you can compile with:
  ```
  make MC_ZLIB=NO
  ```
which will remove the zlib dependency and disables support for gzipped input files.


### Custom Configurations

If you run 'make' without additional parameters MetaCache will be compiled with the default data type settings which support databases with up to 65,535 reference sequences (targets) and k-mer sizes up to 16. This offers a good database space efficiency and is currently sufficient for the complete bacterial, viral and archaea genomes from the NCBI RefSeq.

Using the following compilation options you can compile MetaCache with support for more targets and greater k-mer lengths.

#### number of referece sequences (targets)

* support for up to 65,535 targets (default):
  ```
  make MACROS="-DMC_TARGET_ID_TYPE=uint16_t"
  ```

* support for up to 4,294,967,295 targets (needs more memory):
  ```
  make MACROS="-DMC_TARGET_ID_TYPE=uint32_t"
  ```

* support for more than 4,294,967,295 targets (needs even more memory)
  ```
  make MACROS="-DMC_TARGET_ID_TYPE=uint64_t"
  ```

#### reference sequence lenghts
* support for targets up to a length of 4,294,967,295 windows (default)
  with default settings (window length, k-mer size) no sequence length must exceed 485.3 billion nucleotides
  ```
  make MACROS="-DMC_WINDOW_ID_TYPE=uint32_t"
  ```

* support for targets up to a length of 65,535 windows (needs less memory)
  with default settings (window length, k-mer size) no sequence length must exceed 7.4 million nucleotides
  ```
  make MACROS="-DMC_WINDOW_ID_TYPE=uint16_t"
  ```

#### kmer lengths
* support for kmer lengths up to 16 (default):
  ```
  make MACROS="-DMC_KMER_TYPE=uint32_t"
  ```

* support for kmer lengths up to 32 (needs more memory):
  ```
  make MACROS="-DMC_KMER_TYPE=uint64_t"
  ```

You can of course combine these options (don't forget the surrounding quotes):
  ```
  make MACROS="-DMC_TARGET_ID_TYPE=uint32_t -DMC_WINDOW_ID_TYPE=uint32_t"
  ```

**Note that a database can only be queried with the same variant of MetaCache (regarding data type sizes) that it was built with.**

In rare cases databases built on one platform might not work with MetaCache on other platforms due to bit-endianness and data type width differences. Especially mixing MetaCache executables compiled with 32-bit and 64-bit compilers might be probelematic.


#### disabling zlib support

If you *don't* have the zlib compression library installed and/or want *don't* want gzipped input file support you can compile with:
  ```
  make MC_ZLIB=NO
  ```



## Building Databases


#### Building the Default RefSeq Database
Use the ```metacache-build-refseq``` script to build a MetaCache database based on complete bacterial, viral and archaea genomes from the latest NCBI RefSeq release. Note that the genomes will be downloaded first, which can take some time.

### [Building Custom Databases...](docs/building.md)




## Classification
Once a database (e.g. the standard 'refseq'), is built you can classify reads.
* a single FASTA file containing some reads:
  ```
  ./metacache query refseq my_reads.fa -out results.txt
  ```
* an entire directory containing FASTA/FASTQ files:
  ```
  ./metacache query refseq my_folder -out results.txt
  ```
* paired-end reads in separate files:
  ```
  ./metacache query refseq my_reads1.fa my_reads2.fa -pairfiles -out results.txt
  ```
* paired-end reads in one file (a1,a2,b1,b2,...):
  ```
  ./metacache query refseq my_paired_reads.fa -pairseq -out results.txt
  ```


### [Classification Output Interpretation, Analysis & Formatting Options...](docs/output.md)




## Documentation of Command Line Parameters

* [for mode `build`](docs/mode_build.txt): build database from reference genomes (and write it to disk)
* [for mode `query`](docs/mode_query.txt): query reads against database
* [for mode `build+query`](docs/mode_build_query.txt): build reference database and immediately query reads (mainly recommended for GPU version)
* [for mode `merge`](docs/mode_merge.txt): merge results of independent queries
* [for mode `modify`](docs/mode_modify.txt): add reference genomes to database or update taxonomy
* [for mode `info`](docs/mode_info.txt): obtain information about a database


#### View options documentation from the command line
List available modes:
```
./metacache help
```
or jump directly to a mode's man page with:
```
./metacache help build
./metacache help query
...
```



MetaCache Copyright (C) 2016-2021 [André Müller](https://github.com/muellan) & [Robin Kobus](https://github.com/Funatiq)
This program comes with ABSOLUTELY NO WARRANTY.
This is free software, and you are welcome to redistribute it under certain conditions. See the file 'LICENSE' for details.

[repository]: https://github.com/muellan/metacache
