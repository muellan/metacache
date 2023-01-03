
# MetaCache-GPU


## Example Installation
on an Ubuntu system with NVIDIA Quattro GV100 GPUs and CUDA SDK version 11 installed:
```
sudo apt install -y zlib1g zlib1g-dev
git clone https://github.com/muellan/metacache.git
cd metacache
git submodule update --init --recursive
make gpu CUDA_ARCH=sm_70
```
See below for more details.



## Requirements

### Hardware Requirements

The GPU version of MetaCache requires a CUDA-capable device of the Pascal generation or newer.

In order to be able to build the default database (based on NCBI RefSeq Release 97) with default settings your system will need a total of 120 GB of  GPU memory (e.g. 4x GPUs with 32 GB each).
If you don't have enough GPU memory, you can use [database partitioning](docs/partitioning.md).


### Software Dependencies

* CUDA SDK
  * CUDA >= 11

* Hashtable library [warpcore (github)](https://github.com/sleeepyjack/warpcore) [(HIPC '20 paper)](https://ieeexplore.ieee.org/document/9406635) and sorting library [bb_segsort (github)](https://github.com/Funatiq/bb_segsort). Both repositories are included as submodules and need to be checked out in addition to MetaCache itself. You can do so by calling
  ```git submodule update --init --recursive```

* Support for gzipped FASTA/FASTQ files requires the zlib compression library to be installed on your system.
  On Debian/Ubuntu zlib can be installed with
  `sudo apt install -y zlib1g zlib1g-dev`. If you *don't* have zlib installed or cannot do so you can compile with `make MC_ZLIB=NO`
  which will remove the zlib dependency and disables support for gzipped input files.


## Installation / Compiling

Run `make` in the directory containing the Makefile and set the GPU generation with the `CUDA_ARCH` flag (e.g. `CUDA_ARCH=sm_70` for Quadro GV100):
```
make gpu CUDA_ARCH=sm_70
```

If you don't supply additional parameters MetaCache will be compiled with the default data type settings which support databases with

* up to 4,294,967,295 targets (= reference sequences)
* targets with a length of up to 4,294,967,295 windows (which corresponds to approximately 485.3 billion nucleotides with the default window size of 112)
* kmers with a lengths of up to 16

**A database built by the GPU version can be queried by the corresponding CPU version and vice versa. The only restriction is the available (GPU) memory.**



## Differences to CPU version

MetaCache-GPU allows to **build** distributed databases across multiple GPUs. By default, it will use all available GPUs on the system.
In difference to the [database partitioning](docs/partitioning.md) approach, the reference genomes are automatically distributed across multiple GPUs in a single run. Due to the dynamic distribution scheme and the concurrent execution on the GPUs, two database builds for the same input files will most likely differ. However, this should only have a negligible impact on classification performance.

By default, the **query** mode will use one GPU per database part. It is also possible to utilize more GPUs by replicating the database (see below), which may increase classification speed.

### Build+Query Immediate Mode
Since building databases is significantly faster on the GPU than on the CPU and will often take less than a minute, the [build+query mode](docs/mode_build_query.txt) can be used to build and directly query a database without writing the database to disk.


### Command Line Options

The command line options of the GPU version are similar to the CPU version with a few notable exceptions:

#### mode build

* `-parts <#>` sets the number of GPUs to use (default: all available GPUs).

#### mode query

* `-replicate <#>` enables multiple GPU pipelines (default: 1). Each pipeline occupies one GPU per database part.

#### mode build & mode query

* `-kmerlen` kmer length is limited to 16 (default: 16).
* `-sketchlen` sketch length is limited to 16 (default: 16).
* `-winstride` window stride has to be a multiple of 4 (default: 112).
* `-remove-overpopulated-features` is *not* supported.
* `-remove-ambig-features` is *not* supported.

#### mode info

* submode `locations`is *not* available.
* submode `featurecounts` is *not* available.

#### mode merge

Merging multiple result files will *not* be performed on the GPU and will fall back to the CPU.
