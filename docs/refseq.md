# Building a Database from NCBI RefSeq Genomes

Use the ```metacache-build-refseq``` script to build a MetaCache database based on complete bacterial, viral and archaea genomes from the latest NCBI RefSeq release:

```
git clone https://github.com/muellan/metacache.git 
cd metacache
make
./metacache-build-refseq
```

Note that the genomes will be downloaded first, which can take some time.

Once the default database `refseq.db` is built you can classify reads:
```
./metacache query refseq myReads.fa -out results.txt
./metacache query refseq anyFolderWithFastaOrFastqFiles -out results.txt
./metacache query refseq myReads1.fa myReads2.fa -pairfiles -out results.txt
./metacache query refseq myPairedReads.fa -pairseq -out results.txt
```


## Custom RefSeq Database Builds

MetaCache comes with the following helper scripts:
* `download-ncbi-genomes` downloads NCBI reference genomes
* `download-ncbi-taxonomy` downloads NCBI taxonomy
* `download-ncbi-taxmaps` downloads NCBI accession to taxon ID maps (not needed for the latest NCBI RefSeq releases)


### Example
Build a database with bacteria, viruses, plant and fungal genomes from RefSeq:
```
./download-ncbi-taxonomy tax_folder 

./download-ncbi-genomes refseq/bacteria genomes_folder
./download-ncbi-genomes refseq/viral    genomes_folder
./download-ncbi-genomes refseq/plant    genomes_folder
./download-ncbi-genomes refseq/fungi    genomes_folder

./metacache build myrefseq genomes_folder
```



## See also

* [Output Formatting](output.md)
* [Building Custom Databases](building.md)
* [Using Partitioned Databases](partitioning.md)

### Command Line Options

* [for mode `build`](mode_build.txt)
* [for mode `query`](mode_query.txt)
* [for mode `merge`](mode_merge.txt)


