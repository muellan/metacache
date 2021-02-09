# Building a Custom Reference Genome Database

Make sure that you compile MetaCache with the right data type sizes so that you 
* don't waste memory 
* have support for enough reference sequences per database
* and/or support for long enough reference sequences.

See the [installation instructions](https://github.com/muellan/metacache#detailed-installation-instructions)
for more details.

------------------------

MetaCache's [build mode](mode_build.txt) is used for creating databases.

You need
* taxonomic hierarchy files
* reference genome/scaffold/contig files 
* target-to-taxon mappings that specify each reference genome's taxon


If your machine doesn't have enough RAM to fit an entire database you can use [database partitioning](partitioning.md) to split up the reference genomes into several partitions.


## Taxonomic Hierarchy
You can download the NCBI's taxonomy with an included helper script:
```
download-ncbi-taxonomy ncbi_taxonomy
```
This downloads the taxonomy and puts it in a folder called `ncbi_taxonomy`



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

MetaCache looks for a file named `assembly_summary.txt` (as used by the NCBI) in each refenrece genome input folder.
Such files map reference genome filenames to taxon ids, so that all sequences in a reference genome file get the same taxid.

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

MetaCache looks for such annotations by default, so this does just work.




### 3. NCBI-style `accession2taxid` tab-separated files

#### Use the NCBI's bulk mapping files
The NCBI's accession to taxon mapping files can be downloaded with the included helper script `download-ncbi-taxmaps`.
You should put them into the same folder as the taxonomy: 
```
download-ncbi-taxmaps taxonomy_folder
```
MetaCache looks in the taxonomy folder (determined by build option `-taxonomy`) for them. They are used at the end of the database build phase for all reference sequences that don't have a taxid assigned yet (by any of the other methods).

#### Use custom mapping files
`accession2taxid` files contain 4 tab-separated columns:
```tsv
  accession	accession_version	taxid	gid
  A00002	A00002.1	9913	0
  A00003	A00003.1	9913	0 
  X52700	X52700.1	9771	0  
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

If your genomes are already annotated with taxon ids (see section 2) and/or you also have `assembly_summary.txt` with taxon ids in the genomes folder, than these taxon ids are used with higher priority. 
In case you want the taxon ids from accession2taxid files to supersede those from other sources you can supply the option `-reset-taxa`.
```
metacache build mydb mygenomes_folder -taxonomy ncbi_taxonomy -reset-taxa -taxpostmap mymap.accession2taxid
```
