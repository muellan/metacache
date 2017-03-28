# MetaCache

MetaCache is a taxnomomic classification system for metagenomic read mapping.


## Quick Installation & Usage
This will download MetaCache, compile it, download the complete bacterial, viral and archaea genomes from the latest NCBI RefSeq release (this can take some time) and build a classification database from them:

```
git clone https://github.com/muellan/metacache.git 
cd metacache
make
./metacache-build-refseq
```

Once the default database is built you can classify reads:
  ```
  ./metacache query refseq myReads.fa -out results.txt
  ./metacache query refseq anyFolderWithFastaOrFastqFiles -out results.txt
  ./metacache query refseq -pair_files myReads1.fa myReads2.fa -out results.txt
  ./metacache query refseq -pair_sequences myPairedReads.fa -out results.txt
  ```


## Pre-built Databases
Coming soon...


## Detailed Installation Instructions

#### Requirements
MetaCache itself should compile on any platform for which a C++11 conforming compiler is available.

The helper scripts (for downloading genomes, taxonomy etc.) however require the Bash shell to run. So on Windows you need a working bash executable as well as some common GNU utilities like "awk" and "wget". There is for example "Git Bash" which comes with git for Windows or the 'Windows Subsystem for Linux'.

There are no dependencies on third party libraries.
MetaCache was successfully tested on the following platforms (all 64 bit + 64 bit compilers):
- Ubuntu 14.04 with g++ 4.9, g++ 5.2 or g++ 5.4
- Ubuntu 16.04 with g++ 4.9, g++ 5.3 or g++ 5.4
- Windows 10.1511 64bit with MinGW-w64 g++ 5.1 and MinGW-w64 g++ 5.3
- Windows 10.1607 64bit with MinGW-w64 g++ 5.3
- Windows 10.1607 64bit running Ubuntu 14.04 inside WSL and g++ 5.4.1 

In order to be able to build the default database with default settings your system should have around 32GB of RAM (note that the NCBI RefSeq will still be growing in the near future).


#### Get The Latest Sources
Visit MetaCache's github [repository].


#### Compile
Run 'make' in the directory containing the Makefile. 
This will compile MetaCache with the default data type settings which support databases with up to 65,535 reference sequences (targets) and k-mer sizes up to 16. This offers a good database space efficiency and is currently sufficient for the complete bacterial, viral and archaea genomes from the NCBI RefSeq.

Using the following compilation options you can compile MetaCache with support for more reference sequences and greater k-mer lengths.

##### number of referece sequences (targets)

* support for up to 65,535 reference sequences (default):
  ```
  make MACROS="-DMC_TARGET_ID_TYPE=uint16_t"
  ```

* support for up to 4,294,967,295 reference sequences (needs more memory):
  ```
  make MACROS="-DMC_TARGET_ID_TYPE=uint32_t"
  ```

* support for more than 4,294,967,295 reference sequences (needs even more memory)
  ```
  make MACROS="-DMC_TARGET_ID_TYPE=uint64_t"
  ```

##### reference sequence lenghts
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


##### kmer lengths
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


## Usage
   
#### Build a Reference Database

##### Building the Default RefSeq Database
Use the ```metacache-build-refseq``` script to build a MetaCache database based on complete bacterial, viral and archaea genomes from the latest NCBI RefSeq release. Note that the genomes will be downloaded first, which can take some time. 

##### Custom Builds
Metacache has different modes. The 'build' mode is used for creating databases. See its documentation for more information:
```
./metacache help build
```

MetaCache also comes with these helper scripts:
* ```download-ncbi-genomes``` downloads NCBI reference genomes.
* ```download-ncbi-taxonomy``` downloads NCBI taxonomy.
* ```download-ncbi-taxmaps``` downloads NCBI accession to taxon ID maps.
     Note that these maps are not needed for the latest NCBI RefSeq releases.

Note: In rare cases databases built on one platform might not work with MetaCache on other platforms due to bit-endianness and data type width differences. Especially mixing MetaCache executables compiled with 32-bit and 64-bit compilers might be probelematic.


#### Classification 
Metacache has different modes, one of them is the 'query' mode. Once a database (e.g. the standard 'refseq'), is built you can classify reads.
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
  ./metacache query refseq -pair_files my_reads1.fa my_reads2.fa -out results.txt
  ```
* paired-end reads in one file:
  ```
  ./metacache query refseq -pair_sequences my_paired_reads.fa -out results.txt
  ```

#### View Documentation
The operating manual consists of several text files (one for each mode) located in the 'docs' directory.
Once MetaCache is installed you can also view the documentation with 
```
./metacache help
```
or jump directly to specific topics with
```
./metacache help build
./metacache help query
...
```

MetaCache  Copyright (C) 2016  André Müller
This program comes with ABSOLUTELY NO WARRANTY.
This is free software, and you are welcome to redistribute it under certain
conditions. See the file 'LICENSE' for details.

[repository]: https://github.com/muellan/metacache
