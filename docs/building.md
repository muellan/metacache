# Building a Custom Reference Genome Database

Make sure that you compile MetaCache with the right data type sizes so that you 
* don't waste memory 
* have support for enough reference sequences per database
* and/or long enough reference sequences.

See the [installation instructions](https://github.com/muellan/metacache#detailed-installation-instructions)
for more details.

------------------------

MetaCache's [build mode](build.txt) is used for creating databases.

You need
* taxonomic hierarchy files
* reference genome files 
* target-to-taxon mappings that specify each reference genome's taxon


If your machine doesn't have enough RAM to fit an entire database you can use [database partitioning](partitioning.md) to split up the reference genomes into several partitions.


## Taxonomic Hierarchy
You can download the NCBI's taxonomy with an included helper script:
```
download-ncbi-taxonomy ncbi_taxonomy
```
This will download the taxonomy and put it in a folder called `ncbi_taxonomy`



## Reference Genome Files
Must be in (uncompressed) FASTA or FASTQ format.

You can either specify all input genome files separately:
```
metacache build mydatabase genome1.fna genome2.fnq genome3.fa -taxonomy ncbi_taxonomy
```

or put them all in one folder and then use that folder as input source:
```
metacache build mydatabase genomes_folder -taxonomy ncbi_taxonomy
```



## Target To Taxon Mapping
There are three ways for MetaCache to obtain the taxon of a reference genome.


### 1. NCBI-style `assembly_summary.txt` files

MetaCache looks for a file named `assembly_summary.txt` in each refenrece genome input folder.
Such files map reference genome filenames to taxon ids, so that all sequences in a reference genome file will get the same taxid.

A proper `assembly_summary.txt` file must contain at least:
* first line that starts with `#`
* second line that starts with `#` and all column names separated by tabs
* the first column must contain the file names
* some other column must be named `taxid` (and contain the taxon ids)

Example:
```tsv
# first comment line needs to be there
# assembly_accession	taxid
filename1	438753
filename2	198804
...
```



### 2. Specify the taxon id for each sequence directly in the reference genome files:

```FASTA
>NC_017430.1 | taxid|813 
GCGGCCGCCCGGGAAATTGCTAAAAGATGGGAGCAAAGAGTTAGAGATCTACAAGATAAAGGTGCTGCACGAAAATTATT
...
```

MetaCache looks for such annotations by default, so this will automatically work.




### 3. NCBI-style `accession2taxid` tab-separated files that look like this:

#### Use the NCBI' mapping files
The NCBI's default accession to taxon mappings can be downloaded with the included helper script `download-ncbi-taxmaps`.
You should put them into the same folder as the taxonomy: 
```
download-ncbi-taxmaps taxonomy_folder
```
MetaCache will automatically look in the taxonomy folder (determined by build option `-taxonomy`) for them. They will be used at the end of the database build phase for all reference sequences that don't have a taxid assigned yet (by any of the other methods).

#### Use custom mapping files
```tsv
  accession	accession	version	taxid	gid
  A00002	A00002.1	9913	2
  A00003	A00003.1	9913	3
  X52700	X52700.1	9771	16
  ...
```

The mappings file name(s) go after the parameter `-taxpostmap`:
```
metacache build mydb mygenomes_folder -taxonomy ncbi_taxonomy -taxpostmap mymap.accession2taxid
```

In order for this to work the `accession` or `accession.version` must be contained as sequence id in the
reference genome FASTA or FASTQ files:
```FASTA
>NZ_AAGY02000218.1 Streptococcus pneumoniae TIGR4 ctg00649, whole genome shotgun sequence
TTTGAGCCACTTCGTCTTTAACGGCTTTATTCATAAGCTCTTGTAATTTTTCTTTACTATCAATTACTTCTGATTTTCCG
...
```

