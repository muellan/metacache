Configuration
---------------------------------------------------------------------
```
config.h                     all basic datatype definitions

version.h                    database and project version constants;
                             needed for database compatibility check
```

Mode Starting
---------------------------------------------------------------------
```
main.cpp                     main entry point; selects modes

modes.h                      declares mode starting functions
mode_build.cpp               build database and write to disk
mode_build_query.cpp         build database and directly query it
mode_help.cpp                shows help files from /docs
mode_info.cpp                database property queries
mode_merge.cpp               merge query results from different databases
mode_query.cpp               load database and start classification
```

Classification ("query" mode)
---------------------------------------------------------------------
```
querying.h                   process read files
classification.h             classification starting function
classification.cpp           query database and classify reads
candidates_generation.h      determine top hits in contiguous windows
candidates_structs.h         data structures used for candidate generation
gpu_result_processing.cuh    compact and sort query results, generate candidates on GPU
query_batch.cuh/cu           gpu data handling for query
matches_per_target.h         target->matches map construction

database_query.h             multi-threaded, batched database query functions;
                             also paired-end read handling


printing.h/cpp               classification and analysis output functions
```

Database Operations / Sketching
---------------------------------------------------------------------
```
building.h                   add targets to database
database.h/cpp               wrapper for hashmap + sketchers + taxonomy + auxiliary data

host_hashmap.h               specialized interface to hash_multimap
hash_multimap.h              multi-value hashmap
chunk_allocator.h            default allocator used by hash_multimap

gpu_hashmap.cuh/cu           specialized interface to GPU hashmap
gpu_hashmap_operations.cuh   query and insert kernels for GPU hashmap
dep/warpcore                 submodule containing GPU hashmap
sequence_batch.cuh/cu        gpu data handling for insert

taxonomy.h                   taxonomic tree
taxonomy_io.h/cpp            NCBI taxonomic files -> taxonomic tree

bitmanip.h                   bit-level manipulation functions
dna_encoding.h               ASCII DNA strings -> 2-bit encoded kmers
hash_dna.h                   sketcher classes
hash_int.h                   integer hash functions

sequence_io.h/cpp            reference sequence readers for FASTA/FASTQ files
sequence_view.h              non-owning string view class
```

Analysis
---------------------------------------------------------------------
```
classification_statistics.h  rank-based classification statistics

alignment.h                  construct (semi-global) alignments

stat_confusion.h             thread-safe (binary) confusion statistics
stat_moments.h               accumulators for statistical moments
                             (mean, variance, skewness) and extrema (min/max)
stat_combined.h              combined statistical accumulators
stat_combined.cuh/cu         combined statistical accumulators (GPU version)
```

Utilities
---------------------------------------------------------------------
```
options.h                    default settings for all modes
options.cpp                  command line args parsing for all modes

filesys_utility.h/cpp
cmdline_utility.h/cpp

io_error.h                   I/O exception definitions
io_options.h                 I/O related option types
io_serialize.h               fundamental datatype (de-)serialization functions

batch_processing.h           concurrent batch processing

timer.h                      simple std::chrono based timer class
string_utils.h               string processing utilities (trim)
typename.h                   name demangling functions
```
