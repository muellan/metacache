# Generating Krona Diagrams

The script `krona-from-abundances.py` can generate [Krona](https://github.com/marbl/Krona/wiki/KronaTools) plots from MetaCache abundance results.


1. Query MetaCache datatabase with options `-abundances` and `-abundance-per`:

```
metacache query mydatabase myreads -out classification.txt -abundances abund.txt -abundance-per species
```


2. Then run the script `krona-from-abundances.py` on `abund.txt`:

```
krona-from-abundances abund.txt
```

This script 
  - splits the MetaCache abundance file into two separate files for the original abundance counts and the estimated abundances (here `abund_orig.txt` and `abund_est.txt`)
  - creates a Krona plot HTML file (here `abund_.krona.html`) that can be viewed in a web browser

