# Using Partitioned Databases


## Building

For workstations with limited RAM we provide a script to split the set of all reference genomes into multiple smaller partitions. 
The partitions created by the script are directories with symlinks to the original genome files.
With default settings the resulting MetaCache database built from a partition will be about the same size as the combined file sizes within the partition.

To split the set of all reference genomes into partitions with a given maximum size run

```
metacache-partition-genomes <path to genomes> <partition size in MB>
```

Then run `metacache build` for each of the partition directories.

### Example
Create 16 GB partitions from 40 GB of reference genomes (which results in 3 partitions) and build databases from them:
```
metacache-partition-genomes path/to/mygenomes 16000

metacache build mydb_1 path/to/mygenomes_1 ...
metacache build mydb_2 path/to/mygenomes_2 ...
metacache build mydb_3 path/to/mygenomes_3 ...
```

### See also
* [Building Custom Databases](building.md)
* [build mode man page](build.txt)



## Querying

After database construction

1. query your read dataset(s) against each database partition with `metacache query`
2. merge the individual query results with `metacache merge` to obtain a global result based on all reference genomes


### Example
```
  metacache query mydb_1 myreads.fa -tophits -queryids -lowest species -out res1.txt 
  metacache query mydb_2 myreads.fa -tophits -queryids -lowest species -out res2.txt
  metacache query mydb_3 myreads.fa -tophits -queryids -lowest species -out res3.txt

  metacache merge res1.txt res2.txt res3.txt -taxonomy ncbi_taxonomy
```

### Important!
In order to be mergable, queries must be run with command line options
```
-tophits -queryids -lowest species
```

and must <strong>NOT</strong> be run with options that suppress or alter the default output
like, e.g.: `-nomap`, `-no-summary`, `-separator`, etc.



### See also
* [all query mode command line options](query.txt)
* [all merge mode command line options](merge.txt)
* [output formatting](output.md)


