# Classification Output 

- [Read Mappings](#read-mappings)
- [Read Mapping Output Formatting](#read-mapping-output-formatting-options)
- [Read Mapping Summary](#read-mapping-summary)
- [Abundance Summary](#abundance-summary)
- [Abundance Estimation On A Given Taxonomic Rank](#abundance-estimation-on-taxonomic-rank)
- [Overview Of Mapped-To Reference Genome Locations](#reference-hit-statistics)
- [Ground Truth Based Evaluation](#ground-truth-based-evaluation)




## Read Mappings

MetaCache's default read mapping output format is: 
```read_header | taxon_rank:taxon_name```

This will not be changed in the future to avoid breaking anyone's pipelines. Command line options won't change in the near future for the same reason. 

#### Example Output

All lines that don't contain read mappings start with `#`.

The first column contains the read ids obtained from the input files, the second column contains the classification result.
The first read in the example has id `A_hydrophila_HiSeq.20480` and was mapped to species Aeromonas hydrophila.

```
# Reporting per-read mappings (non-mapping lines start with '# ').
# Only the lowest matching rank will be reported.
# Classification will be constrained to ranks from 'sequence' to 'domain'.
# Classification hit threshold is 5 per query
# At maximum 2 classification candidates will be considered per query.
# Using 64 threads
# TABLE_LAYOUT: query_header    |   rank:taxname
# /mnt/ssd/metagenomics/HiSeq_timing.fna
A\_hydrophila\_HiSeq.20480    |   species:Aeromonas hydrophila
A\_hydrophila\_HiSeq.20481    |   species:Aeromonas hydrophila
A\_hydrophila\_HiSeq.20482    |   genus:Aeromonas
A\_hydrophila\_HiSeq.20483    |   sequence:NZ\_CP007518.2
A\_hydrophila\_HiSeq.20484    |   sequence:NZ\_CP007518.2
A\_hydrophila\_HiSeq.20485    |   genus:Aeromonas
A\_hydrophila\_HiSeq.20486    |   sequence:NZ\_CP007518.2
A\_hydrophila\_HiSeq.20487    |   genus:Aeromonas
A\_hydrophila\_HiSeq.20488    |   species:Aeromonas hydrophila
A\_hydrophila\_HiSeq.20489    |   sequence:NZ\_CP007518.2
A\_hydrophila\_HiSeq.20490    |   --
...
```




## Read Mapping Output Formatting Options

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




## Read Mapping Summary

#### Example Output

```
# queries: 10000
# time:    14 ms
# speed:   4.24286e+07 queries/min
# unclassified: 19.02% (1902)
# classified:
#   sequence   9.54% (954)
#   subspecies 12.72% (1272)
#   species    65.95% (6595)
#   genus      78.11% (7811)
#   family     78.87% (7887)
#   order      78.93% (7893)
#   class      79.17% (7917)
#   phylum     79.41% (7941)
#   kingdom    79.41% (7941)
#   domain     80.98% (8098)
#   root       80.98% (8098)
```

#### Related Command Line Options

| command line argument | purpose / effect                                              |
| --------------------  | ------------------------------------------------------------- |
| ```-no-map```         | suppresses detailed, per-read mappings                        |
| ```-no-summary```     | suppresses the mapping summary                                |
| ```-lowest <rank>```  | lowest taxonomic rank for classification, default: 'sequence' |




## Abundance Summary
with command line option `-abundances [<filename>]`. If no filename is given the list is appended to the mappings output.

```
# query summary: number of queries mapped per taxon
# rank:name                     | taxid   | reads | abundance
domain:Bacteria                 | 2       | 30    |  0.3%
phylum:Proteobacteria           | 1224    | 16    |  0.16%
phylum:Firmicutes               | 1239    | 1     |  0.01%
class:Gammaproteobacteria       | 1236    | 28    |  0.28%
order:Bacillales                | 1385    | 4     |  0.04%
order:Corynebacteriales         | 85007   | 1     |  0.01%
order:Enterobacterales          | 91347   | 3     |  0.03%
family:Enterobacteriaceae       | 543     | 6     |  0.06%
family:Streptococcaceae         | 1300    | 1     |  0.01%
genus:Xanthomonas               | 338     | 233   |  2.33%
genus:Streptococcus             | 1301    | 2     |  0.02%
genus:Bacillus                  | 1386    | 5     |  0.05%
subgenus:Bacillus cereus group  | 86661   | 178   |  1.78%
species:Xanthomonas citri       | 346     | 3     |  0.03%
species:Vibrio cholerae         | 666     | 4     | 32.04%
species:Streptococcus pneumoniae| 1313    | 768   | 17.68%
species:Streptococcus pyogenes  | 1314    | 2     |  9.02%
species:Bacillus cereus         | 1396    | 70    |  8.7%
sequence:NC_005707.1            | 222523  | 1     |  0.03%
sequence:NC_003909.8            | 222523  | 9     |  0.09%
sequence:NC_002944.2            | 262316  | 1     |  0.01%
sequence:NC_003997.3            | 198094  | 7     |  0.07%
sequence:NC_004722.1            | 226900  | 7     |  0.07%
...
```



## Abundance Estimation On Taxonomic Rank

with command line option `-abundance-per <taxonomic rank>`. 

If the raw abundance list (see above) is redirected to a dedicated file with `-abundances <filename>` then the estimation list is also written to the same file.

#### Example Output: `abundance-per species`

```
# estimated abundance (number of queries) per species
# rank:name                     | taxid | queries | abundance
family:Pasteurellaceae          | 712   | 1.04485 |  0.0104485%
species:Xanthomonas campestris  | 339   | 2.51262 |  0.0251262%
species:Xanthomonas citri       | 346   | 412.069 | 34.12069%
species:Xanthomonas oryzae      | 347   | 7.53785 |  0.0753785%
species:Sinorhizobium meliloti  | 382   | 1.02339 |  0.0102339%
species:Escherichia coli        | 562   | 3.32452 |  0.0332452%
species:Shigella flexneri       | 623   | 3.32452 |  0.0332452%
species:Yersinia pestis         | 632   | 1.32981 |  0.0132981%
species:Vibrio cholerae         | 666   | 957.083 | 28.57083%
species:Streptococcus pyogenes  | 1314  | 2.03794 |  0.0203794%
species:Enterococcus faecalis   | 1351  | 1.01568 |  0.0101568%
species:Bacillus anthracis      | 1392  | 21.2461 |  0.212461%
species:Bacillus cereus         | 1396  | 264.059 | 12.64059%
species:Mycobacterium avium     | 1764  | 2.02304 |  0.0202304%
species:Salmonella enterica     | 28901 | 6.64905 |  0.0664905%
...

```

#### Example Output: `abundance-per phylum`
```
# estimated abundance (number of queries) per phylum
# rank:name           | taxid  | queries | abundance
phylum:Proteobacteria | 1224   | 1395.9  | 73.959%
phylum:Firmicutes     | 1239   | 1236.08 | 12.3608%
phylum:Actinobacteria | 201174 | 2.02304 | 0.0202304%
...
```



## Reference Hit Statistics

The command line option `-hits-per-ref [<filename>]` generates an overview of all locations in the reference genomes that the query sequences (reads/read pairs) were mapped to. If no filename is given the list is appended to the mappings output.
Note that the generated output can get **very** large!

#### Example Output

```
# --- list of hits for each reference sequence ---
# window start position within sequence = window_index * window_stride(=113)
# TABLE_LAYOUT:  sequence       |        windows_in_sequence    |       queryid/window_index:hits/window_index:hits/...,queryid/...
sequence:NC_003098.1    |       18041   |       7993/8:12,7008/11:6/12:10,7467/24:10
sequence:NC_005004.1    |       216     |       6957/60:5
sequence:NC_003047.1    |       32338   |       5295/13194:6
sequence:NC_004557.1    |       24773   |       6114/1562:5
sequence:NC_004070.1    |       16819   |       7387/179:4/180:12
sequence:NC_002932.3    |       19071   |       2371/1267:2/1268:4,8667/17990:5
sequence:NC_002162.1    |       6653    |       8667/3076:5
sequence:NC_003922.1    |       575     |       9474/120:10,9126/143:8/144:8,9247/144:6/145:7,9410/144:5/145:6,9959/146:3/147:3,9053/153:2/154:4,9195/216:7,9460/316:3/317:2,9465/354:7,9800/363:12,9652/366:8,9397/379:6
sequence:NC_000913.3    |       41077   |       8350/2013:5/2014:5,5994/2020:6,25/7580:6,7286/13390:6/13391:6,7162/17432:1/17433:13,494/24155:2/24156:5,2462/26280:7/26281:4,7304/30132:1/30133:6,3763/37011:10
sequence:NC_003197.2    |       42987   |       8647/24789:8,52/31603:3/31604:13,46/38913:8,563/38913:5
sequence:NC_002662.1    |       20935   |       7127/17564:13
sequence:NC_002940.2    |       15035   |       8752/96:5/97:7
...
```

This means that the reference sequence 'NC_002932.3' is divided into 19071 windows and that the input read with id '2371' produced 2 hits in window number '1267' and 4 hits in windows 1268 and the input read with id '8667' produced 5 hits in window number 17990 of that reference sequence.

The positions in the genome can be calculated by multiplying a window index with the window stride, which is shown in the header of the table (in this case the stride is set to MetaCache's default of 113).

This table can be used for, e.g., coverage analyses.




## Ground Truth Based Evaluation

can be turned on with the option `-precision`, which has the effect that precision, sensitivity, etc. will be reported in the mapping summary.
The ground truth taxon can also be shown in the read mappings (next to the classification) with option `-ground-truth`.
Note that in order to obtain the ground truth for an input read, each read (pair) needs to be annotated in its header with either an `taxid|<number>` entry or a sequence id (e.g., `NC_002505.1`) of a reference sequence that is present in the database.


#### Example Read Mappings Output

run with options `-precision -ground-truth -omit-ranks -taxids`.
Note that negative reference sequences in the database have nagative taxon ids in order to not interfere with the NCBI's taxon ids.

```
# Reporting per-read mappings (non-mapping lines start with '# ').
# Only the lowest matching rank will be reported.
# Classification will be constrained to ranks from 'sequence' to 'domain'.
# Classification hit threshold is 5 per query
# At maximum 2 classification candidates will be considered per query.
# Using 4 threads
# TABLE_LAYOUT: query_header    |   truth_taxname(truth_taxid)  |   taxname(taxid)
# analysis/HiSeq.fna
V_cholerae_HiSeq.195766 |   Vibrio cholerae(666)    |   NC_002505.1(-15)
V_cholerae_HiSeq.196272 |   Vibrio cholerae(666)    |   NC_002506.1(-16)
V_cholerae_HiSeq.196439 |   Vibrio cholerae(666)    |   Proteobacteria(1224)
V_cholerae_HiSeq.196534 |   Vibrio cholerae(666)    |   NC_002505.1(-15)
V_cholerae_HiSeq.196565 |   Vibrio cholerae(666)    |   Vibrio cholerae(666)
V_cholerae_HiSeq.198656 |   Vibrio cholerae(666)    |   --
V_cholerae_HiSeq.200147 |   Vibrio cholerae(666)    |   Vibrio cholerae(666)
V_cholerae_HiSeq.200184 |   Vibrio cholerae(666)    |   Gammaproteobacteria(1236)
V_cholerae_HiSeq.201960 |   Vibrio cholerae(666)    |   NC_002505.1(-15)
V_cholerae_HiSeq.202496 |   Vibrio cholerae(666)    |   NC_002506.1(-16)
...
```

Same with options `-precision -ground-truth -omit-ranks -taxids-only`:

```
...
# TABLE_LAYOUT: query_header    |   truth_taxid  |   taxid
# analysis/HiSeq.fna
V_cholerae_HiSeq.195766 |   666    |   -15
V_cholerae_HiSeq.196272 |   666    |   -16
V_cholerae_HiSeq.196439 |   666    |   1224
V_cholerae_HiSeq.196534 |   666    |   -15
V_cholerae_HiSeq.196565 |   666    |   666
V_cholerae_HiSeq.198654 |   666    |   0
V_cholerae_HiSeq.200147 |   666    |   666
V_cholerae_HiSeq.200184 |   666    |   1236
V_cholerae_HiSeq.201960 |   666    |   -15
V_cholerae_HiSeq.202496 |   666    |   -16
...
```


#### Example Summary Output
```
# queries: 10000
# time:    14 ms
# speed:   4.24286e+07 queries/min
# unclassified: 19.02% (1902)
# classified:
#   sequence   9.54% (954)
#   subspecies 12.72% (1272)
#   species    65.95% (6595)
#   genus      78.11% (7811)
#   family     78.87% (7887)
#   order      78.93% (7893)
#   class      79.17% (7917)
#   phylum     79.41% (7941)
#   kingdom    79.41% (7941)
#   domain     80.98% (8098)
#   root       80.98% (8098)
# ground truth known:
#   sequence   0% (0)
#   subspecies 20% (2000)
#   species    100% (10000)
#   genus      100% (10000)
#   family     100% (10000)
#   order      100% (10000)
#   class      100% (10000)
#   phylum     100% (10000)
#   kingdom    100% (10000)
#   domain     100% (10000)
#   root       100% (10000)
# correctly classified:
#   sequence   0
#   subspecies 286
#   species    5032
#   genus      7767
#   family     7874
#   order      7884
#   class      7904
#   phylum     7935
#   kingdom    7935
#   domain     8098
#   root       8098
# precision (correctly classified / classified) if ground truth known:
#   sequence   0%
#   subspecies 15.3351%
#   species    76.1271%
#   genus      99.3731%
#   family     99.7846%
#   order      99.8354%
#   class      99.8358%
#   phylum     99.9244%
#   kingdom    99.9244%
#   domain     100%
#   root       100%
# sensitivity (correctly classified / all) if ground truth known:
#   sequence   0%
#   subspecies 14.3%
#   species    50.32%
#   genus      77.67%
#   family     78.74%
#   order      78.84%
#   class      79.04%
#   phylum     79.35%
#   kingdom    79.35%
#   domain     80.98%
#   root       80.98%
```



## [All query mode command line options...](mode_query.txt)
