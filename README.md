# Map reads with suffix-array and Burrows-Wheeler Transformation (BWT) algorithm

This program is **not intended for practical purposes** but to understand how NGS reads are mapped to the reference genome. Efficiency of suffix-array-based and BWT-based can be compared. This script will simulate both NGS reads and the refrence sequence. It is adapted from https://github.com/rayguo233/read-mapping/.

## Useage
for details: 
```{bash}
map_reads.py -h
```

```{bash}
map_reads.py -n [number of NGS reads] -l [length of NGS reads (e.g. 150 bp)] -g [length of simulated genome] -a [algorithm to map reads (0|1|2)]
```
