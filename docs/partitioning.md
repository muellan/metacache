# Using Partitioned Databases

MetaCache offers two ways of splitting larger databases into multiple smaller partitions.
This can by useful either for workstations with a limited amount of RAM or in case you want to build very large databases in the range of 100s of GBs or even TBs.

MetaCache can either automatically partition databases during the build phase or the user can
manually build separate databases. The output of querying several or all partitions can be
merged afterwards to obtain final classication results.

Automatic partitioning is the default for the GPU version of MetaCache where each database 
partition resides on a different GPU and the partition sizes are determined/limited
by the GPUs' VRAM sizes.

The CPU version builds a single, monolithic database by default, but
automatic partitioning can be chosen by suppying options `-parts <#>`
and/or `-max-part-size <GB>`.



## Building

### Automatic Partitioning

Automatic partitioning produces databases that consist of one `.meta` file (metadata)
and multiple `.cache#` files (the actual partitions) where `#` is a number from 0 to
one less than the number of partitions.

Maximum partition size values that are much smaller than the reference sequence sizes
might produce partitions that are larger than the size threshold, since any partition must
contain at least one reference. 

#### GPU-Version
- The number of parts should be equal or smaller to the number of
available GPUs. By default all available GPUs in a system are automatically used
and each GPU will hold one partition.
- The `-max-part-size` option is not available for the GPU version of MetaCache.

#### CPU-Version
Reference sequences will be processed in parallel followed by a consolidation
phase which generates the database parts in a way that minimizes feature collisions
(hash table bucket overflows) if multiple partitions are requested.

**Important!**
Automatic partitioning can only be used, if there's enough memory to hold
all database partitions at once (similar to a monolothic single part database).
If you only have enough memory to fit a single partition, you need to manually
partition the reference genomes (see below).

Querying can always be done one partition at a time on a system that can only fit a single partition
into memory; to obtain the final classification the individual query results of all partitions 
have to merged afterwards (see below). 


#### Examples

In the following examples all reference sequence files are in directory `genomes_folder` 
and its subdirectories and the NCBI taxonomy files are in directory `ncbi_taxonomy`.

* Build a database with exactly 4 partitions (no limit on sizes of each partition):
  ```
  metacache build mydatabase genomes_folder -taxonomy ncbi_taxonomy -parts 4
  ```

* Build a database with a maximum size of 30GB (no limit on number of partitions):
  ```
  metacache build mydatabase genomes_folder -taxonomy ncbi_taxonomy -max-part-size 30
  ```

* Build a database with a minimum of 4 partitions and a maximum size of 30GB for each partition:
  ```
  metacache build mydatabase genomes_folder -taxonomy ncbi_taxonomy -parts 4 -max-part-size 30
  ```
  If 4 partitions are not enough to fit all reference sequences while keeping
  the 30GB size limit, one or more additional partitions will be created.


### Manual Partitioning

The script `metacache-partition-genomes` can be used to split the set of all reference genomes 
into multiple smaller partitions based on the sequence files' sizes. 
The partitions created by the script are subdirectories with symlinks to the original genome files.
With default settings the resulting MetaCache database built from a partition will be about the same size as the combined uncompressed FASTA file sizes within the partition.

To split the set of all reference genomes into partitions with a given maximum size, run:

```
metacache-partition-genomes <path to genomes> <partition size in MB>
```

Then run `metacache build` for each of the partition directories.


#### Example
Create 16 GB partitions from 40 GB of reference genomes (which results in 3 partitions) and build databases from them:
```
metacache-partition-genomes path/to/mygenomes 16000

metacache build mydb_1 path/to/mygenomes_1 -taxonomy ncbi_taxonomy
metacache build mydb_2 path/to/mygenomes_2 -taxonomy ncbi_taxonomy
metacache build mydb_3 path/to/mygenomes_3 -taxonomy ncbi_taxonomy
```

(This assumes that the NCBI taxonomy files are in directory `ncbi_taxonomi`.)


### See also
* [Building Custom Databases](building.md)
* [build mode man page](mode_build.txt)




## Querying


### Sequential Querying (one partition at a time) in order to save memory

1. **query** your read dataset(s) against each database partition with `metacache query`
2. **merge** the individual query results with `metacache merge` to obtain a global result based on all reference genomes


#### Example: Sequential querying of 3 individual databases (as produced by manually partitioning)
```
  metacache query mydb_1 myreads.fa -tophits -queryids -lowest species -out query1.txt 
  metacache query mydb_2 myreads.fa -tophits -queryids -lowest species -out query2.txt
  metacache query mydb_3 myreads.fa -tophits -queryids -lowest species -out query3.txt

  metacache merge query1.txt query2.txt query3.txt -taxonomy ncbi_taxonomy -out result.txt
```

#### Example: Sequential querying of an automatically partitioned database with 3 parts
If the name of a partition file is given explicitly with the `.cache#` extension, only that
file will be loaded and queried. (If the extension was left out, e.g. just `mydb` was given,
in the example below, then all partitions would be loaded at the same time.)
```
  metacache query mydb.cache0 myreads.fa -tophits -queryids -lowest species -out query1.txt 
  metacache query mydb.cache1 myreads.fa -tophits -queryids -lowest species -out query2.txt
  metacache query mydb.cache2 myreads.fa -tophits -queryids -lowest species -out query3.txt

  metacache merge query1.txt query2.txt query3.txt -taxonomy ncbi_taxonomy -out result.txt
```


#### Important!
In order to be mergable, queries must be run with command line options
```
-tophits -queryids -lowest species
```

and must <strong>NOT</strong> be run with options that suppress or alter the default output content or
output formatting like, e.g.: `-no-map`, `-no-summary`, `-separator`, etc.



### Parallel Querying

of all database parts at once.
This works only for automatically partitioned databases.


#### Example: Query automatically partitioned database in one go
All parts will be loaded into memory and queried just like a monolithic single-partition database.
```
  metacache query mydb myreads.fa -out result.txt 
```



### See also
* [all query mode command line options](mode_query.txt)
* [all merge mode command line options](mode_merge.txt)
* [output formatting](output.md)


