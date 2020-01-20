# AFS-MetaCache: Food Ingredient Detection & Abundance Analysis

MetaCache is a classification system for mapping (short or long) reads from metagenomic samples to their most likely taxon of origin. It uses locality sensitive hashing to quickly identify candidate regions within one or multiple reference genomes. A read is then classified based on the similarity to those regions. 


* [MetaCache Github Repository](https://github.com/muellan/metacache)
* [**MetaCacheSpark**: Apache Spark&trade; implementation of MetaCache for big data clusters](https://github.com/jmabuin/MetaCacheSpark)



## Installation

Clone MetaCache's Github repository and compile the executable:

```
git clone https://github.com/muellan/metacache.git 
cd metacache
make MACROS="-DMC_TARGET_ID_TYPE=uint32_t"
```

The data type configuration in the last line makes sure that MetaCache can handle up to 4,294,967,295 reference sequences per database.
This is necessary because many eukaryotic reference genome files contain _a lot_ of genomic sequences.
If you know that your database will include less than 65,535 sequences you can safe memory by just calling `make` without the `MACROS=` part. 



## Building Reference Databases

#### Example

The following snippet downloads the NCBI taxonomic hierarchy and builds a database `afs.db` from all genomes in `genomes_folder`. The reference sequence to taxon id mapping is supplied in the file `taxa.accession2taxid` which
conforms to the NCBI's `accession2taxid` format.
```
download-ncbi-taxonomy ncbi_tax

metacache build afs genomes_folder \
          -taxonomy ncbi_tax -taxpostmap taxa.accession2taxid \
          -remove-overpopulated-features 
```

It is important that you supply the option `-remove-overpopulated-features` if you add large eukaryotic genomes. This will improve classification accuracy and runtime performance.


#### For more information see
* [Building custom databases](building.md)
* [Using partitioned databases](partitioning.md)



## Querying & Abundance Estimation

#### Example

This queries all reads in all FASTA and FASTQ files in `myreads_folder` against database `afs.db`. The individual read mappings are written to `read_mappings.txt` and the abundance analysis to `abund.txt`.
```
metacache query afs myreads_folder -pairfiles \
          -max-cand 4 -hitmin 4 -hitdiff 80 \
          -out read_mappings.txt -splitout -abundances abund.txt -abundance-per species
```

The option `-pairfiles` means that each mate of a read pair will be read from a separate file (e.g., `fileA_1.fa` + `fileA_2.fa`, `fileB_1.fa` + `fileB_2.fa`, ...). If both mates of your paired-end reads are in one file (1+2, 3+4, ...) then use the option `-pairseq` instead.
If the option `-splitout` is given, mapping and abundance results will be written to a separate file for each input file(-pair).

**Important:** For abundance analysis tasks involving large eukaryotic genomes, make sure to supply options 
```
-max-cand 4 -hitmin 4 -hitdiff 80
```
for optimal results. This tells MetaCache to consider the 4 best matching candidates per read (the default is 2, which is fine for bacteria). It also makes sure that the best matching candidate wins over the lowest common ancestor of all candidates if the input read has at least 80% features more in common with this best candidate than it has with the 2nd best.



## More Documentation

* [Output Formatting](output.md)
* [Building custom databases](building.md)
* [Using partitioned databases](partitioning.md)


### All Command Line Options

* [for mode `build`](mode_build.txt): build database(s) from reference genomes
* [for mode `query`](mode_query.txt): query reads against databases
* [for mode `merge`](mode_merge.txt): merge results of independent queries
* [for mode `modify`](mode_modify.txt): add reference genomes to database or update taxonomy
* [for mode `info`](mode_info.txt): get information about a database

