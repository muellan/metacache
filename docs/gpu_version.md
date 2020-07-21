# MetaCache-GPU

## Installation Instructions

#### Requirements

The GPU version of MetaCache requires a CUDA-capable device of the Pascal generation or newer and either:

* CUDA >= 11
* CUDA 10.2 and a self-provided version of [CUB](https://github.com/NVlabs/cub)

Make sure to adjust the Makefile to the GPU generation you want to use by setting the `-arch` flag (e.g. `-arch=sm_70` for Quadro GV100). You also have to set the include path for CUB if your CUDA version is below CUDA 11.

MetaCache-GPU depends on the hashtable implementation of [warpcore](https://github.com/sleeepyjack/warpcore) and the sorting algorithm [bb_segsort](https://github.com/Funatiq/bb_segsort). Both repositories are included as submodules and need to be checked out in addition to MetaCache itself. You can do so be calling

```git submodule update --init --recursive```

In order to be able to build the default database (based on NCBI RefSeq Release 97) with default settings your system will need a total of 120 GB of GPU memory (e.g. 4x GPUs with 32 GB each).
If you don't have enough GPU memory, you can use [database partitioning](docs/partitioning.md).

#### Compile
Run '`make gpu_release`' in the directory containing the Makefile.
This will compile MetaCache-GPU with support for:

* up to 4,294,967,295 reference sequences
* targets up to a length of 4,294,967,295 windows
* kmer lengths up to 16

This corresponds to the CPU version compiled with `make MACROS="-DMC_TARGET_ID_TYPE=uint32_t"`

**Note that a database build by the GPU version can be queried by the corresponding CPU version and vice versa. The only restriction is the available (GPU) memory.**


## Differences to CPU version

MetaCache-GPU allows to **build** distributed databases across multiple GPUs.
In difference to the [database partitioning](docs/partitioning.md) approach, the program distributes the reference genomes automatically across the GPUs in a single run. Due to the dynamic distribution scheme and the concurrent execution on the GPUs, two database builds for the same input files will most likely differ. However, this should have only a small impact on classification performance.

In order to **query** a multi-GPU database make sure to set the same number of GPUs when using the query mode. Note, that only a small number of threads is needed to saturate the GPU query pipeline.

#### Command Line Options

The command line options of the GPU version are similar to the CPU version with a few notable exceptions:

##### mode build & mode query

* `-gpus <#>` sets the maximum number of GPUs to use (default: all available GPUs).
* `-winlen` window length is limited to 127 (default: 127).
* `-winstride` window stride has to be multiple of 4 (default: 112).
* `-remove-overpopulated-features` is not supported.
* `-remove-ambig-features` is not supported.

##### mode info

* feature map is not available.
* feature counts are not available.

##### mode merge

* merging on GPU is not available and will fall back to CPU version.
