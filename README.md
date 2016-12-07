# MetaCache taxonomic classification system

## Installation Instructions

#### Requirements
MetaCache itself should compile on any platform with a C++11 conforming compiler.

The helper scripts (for downloading genomes, taxonomy etc.) however require the Bash shell to run. So on Windows you need a working bash executable, for example "Git Bash" which comes with git for Windows as well as some common GNU utilities like 'awk' and 'wget'.

There are no dependencies on third party libraries.
Compilation was successfully tested on the following platforms:
- Ubuntu 14.04 with g++ 4.9 and g++ 5.2
- Ubuntu 16.04 with g++ 4.9 and g++ 5.3
- Windows 10.1511 with MinGW-w64 g++ 5.1 and MinGW-w64 g++ 5.3
- Windows 10.1607 with MinGW-w64 g++ 5.3


#### Get The Latest Sources
Visit MetaCache's github [repository].


#### Compile
Run 'make' in the directory containing the Makefile.

The default supports databases with up to 65535 reference sequences (genomes) with a length
of 8,388,480 nucleotides each and k-mer sizes up to 16. This offers a good database space efficiency and should be sufficient for bacterial and viral genomes.
Using the following compilation options you can compile MetaCache with support for more reference sequences of longer size.

##### number of referece sequences (genomes)

* support for up to 255 reference sequences (needs less memory):
  ```
  make MACROS="-DMC_GENOME_ID_TYPE=uint8_t"
  ```

* support for up to 65535 reference sequences (default):
  ```
  make MACROS="-DMC_GENOME_ID_TYPE=uint16_t"
  ```

* support for up to 4,294,967,295 reference sequences (needs more memory):
  ```
  make MACROS="-DMC_GENOME_ID_TYPE=uint32_t"
  ```

* support for more than 4,294,967,295 reference sequences (needs even more memory)
  ```
  make MACROS="-DMC_GENOME_ID_TYPE=uint64_t"
  ```

##### reference sequence lenghts

* support for genomes up to a length of 65536 windows (default)
  if window size is 128 (default) then max. genome length supported is 8,388,480 nucleotides
  ```
  make MACROS="-DMC_WINDOW_ID_TYPE=uint16_t"
  ```

* support for genomes up to a length of 4,294,967,295 windows (needs more memory)
  if window size is 128 (default) then max. genome length supported is 549,755,813,760 nucleotides
  ```
  make MACROS="-DMC_WINDOW_ID_TYPE=uint32_t"
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
  make MACROS="-DMC_GENOME_ID_TYPE=uint32_t -DMC_WINDOW_ID_TYPE=uint32_t"
  ```

**Note that a database can only be queried with the variant of MetaCache that it was built with.**


## Usage
   
#### Build a Reference Database

##### Building the Default RefSeq Database
Use the ```metacache-build-refseq``` script to build a MetaCache database based on complete bacterial, viral and archaea genomes from the latest NCBI RefSeq release. Note that the genomes will be downloaded first, which can take some time. Currently you can choose between three default settings: standard, big and small.
```
./metacache-build-refseq standard
./metacache-build-refseq big
./metacache-build-refseq small
```

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


#### Classification TL;DR 
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
