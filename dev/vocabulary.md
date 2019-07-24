MetaCache Vocabulary
=====================
```
 term        explanation

 taxon       node in the tree of life
 rank        named level in the tree of life (species, genus, family , ...)

 query       read or pair of reads to be classified
 target      (genomic) reference sequence

 window      sequence of w consecutive characters used for constructing a sketch
 kmer        sequence of k < w consecutive characters within a window
 feature     hash of a kmer
 sketch      set of s <= w-k+1 features
 sketcher    takes a window and returns a sketch

 location    {target id , window index within target}

 database    feature->location hash table + taxonomy + some auxiliary data
```
