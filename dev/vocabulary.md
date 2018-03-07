MetaCache Vocabulary
=====================
```
 term        explanation                                                           

 taxon       node in the tree of life                                       
 rank        named level in the tree of life (species, genus, family , ...) 

 target      (genomic) reference sequence                       
 query       short read or pair of short reads to be classified 

 window      sequence of 'w' consecutive characters used for constructing a sketch 
 kmer        sequence of 'k' consecutive characters                                
 feature     hash of a kmer                                                        
 sketch      set of 's' features                                                   
 sketcher    takes a window and returns a sketch                                   

 location    {target id , window index within target} 

 database    feature -> location hash table + taxonomy + some auxiliary data 
```
