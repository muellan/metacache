# MetaCache taxonomic classification system

## Installation Instructions

#### Requirements
MetaCache itself should compile on any platform for which a C++11
conforming compiler is available.

The helper scripts (for downloading genomes, taxonomy etc.) however require
the Bash shell to run. So on Windows you need a working bash executable,
for example "Git Bash" which comes with git for Windows as well as some 
common GNU utilities like 'awk' and 'wget'.

There are no dependencies on third party libraries.
Compilation was successfully tested on the following platforms:
- Ubuntu 14.04 with g++ 4.9
- Ubuntu 14.04 with g++ 5.2
- Ubuntu 16.04 with g++ 4.9
- Ubuntu 16.04 with g++ 5.3
- Windows 10.1511 with MinGW-w64 g++ 5.1
- Windows 10.1511 with MinGW-w64 g++ 5.1
- Windows 10.1511 with MinGW-w64 g++ 5.3

#### Get The Latest Sources
Visit MetaCache's github repository [repo].


#### Compile
Run 'make' in the directory containing the Makefile.


## Usage

#### Build A Reference Database
Use the ```metacache-build-refseq``` script to build a MetaCache database based
on complete genomes from the latest NCBI RefSeq or Genbank. Note that the
genomes will be downloaded first, which can take some time.
Currently you can choose between three default settings: standard, big and small:
```
./metacache-build-refseq standard
./metacache-build-refseq big
./metacache-build-refseq small
```

If you want full control over the individual steps, there are several
helper scripts:
- ```download-ncbi-genomes``` downloads NCBI reference genomes.
- ```download-ncbi-taxonomy``` downloads NCBI taxonomy.
- ```download-ncbi-taxmap``` downloads NCBI accession to taxon ID maps.
  These maps are not needed for the latest NCBI RefSeq releases.

#### Classification TL;DR 
Metacache has different modes, one of them is the 'query' mode for classification.
Once the standard 'refseq' database is built you can classify reads:
```
./metacache query refseq my_reads.fa -out results.txt
```
```
./metacache query refseq directory_with_read_files -out results.txt
```
paired-end reads in files:
```
./metacache query refseq -pair_files my_reads1.fa my_reads2.fa -out results.txt
```
paired-end reads in one file:
```
./metacache query refseq -pair_sequences my_paired_reads.fa -out results.txt
```

#### View Documentation
The operating manual consists of several text files (one for each mode)
located in the 'docs' directory.
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
This is free software, and you are welcome to redistribute it
under certain conditions. See the file 'LICENSE' for details.

[repo]: https://github.com/muellan/metacache
