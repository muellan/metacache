# Output Formatting

MetaCache's default read mapping output format is: 
```read_header | taxon_rank:taxon_name```

This will not be changed in the future to avoid breaking anyone's pipelines. Command line options won't change in the near future for the same reason. 

The following tables show some of the possible mapping layouts with their associated command line arguments.

| mapping layout                                   | command line arguments                   |
| -----------------------------------------------  | ---------------------------------------- |
| ```read_header │ rank:taxon_name```              |                                          |
| ```read_header │ rank:taxon_name(taxon_id)```    | ```-taxids```                            |
| ```read_header │ rank │ taxon_id```              | ```-taxids-only -separate-cols```        |
| ```read_header │ rank │ taxon_name```            | ```-separate-cols```                     |
| ```read_header │ rank │ taxon_name │ taxon_id``` | ```-taxids -separate-cols```             |
| ```read_header │ taxon_id```                     | ```-taxids-only -omit-ranks```           |
| ```read_header │ taxon_name```                   | ```-omit-ranks```                        |
| ```read_header │ taxon_name(taxon_id)```         | ```-taxids -omit-ranks```                |
| ```read_header │ taxon_name │ taxon_id```        | ```-taxids -omit-ranks -separate-cols``` |
| ```read_header │ rank:taxon_id```                | ```-taxids-only```                       |


Note that the default lowest taxon rank is "sequence". Sequence-level taxon ids have negative numbers in order to not interfere with NCBI taxon ids. Each target sequence (reference genome) is added as its own taxon below the lowest known NCBI taxon for that sequence. If you do not want to classify at sequence-level, you can set a higher rank as lowest classification rank with the ```-lowest``` command line option: ```-lowest species``` or ```-lowest subspecies``` or ```-lowest genus```, etc.


### [All query mode command line options...](mode_query.txt)
