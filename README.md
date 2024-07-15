# Map reads with suffix-array Burrows-Wheeler Transformation (BWT) algorithm

This program is not intended for practical purposes but to understand how NGS reads are mapped to the reference genome. Code is adapted from https://github.com/rayguo233/read-mapping/.

## Useage
for help ```map_reads.py -h```

```{bash}
map_reads.py -n [number of NGS reads] -l [length of NGS reads (e.g. 150 bp)] -g [length of simulated genome] -a [algorithm to map reads (0|1|2)]
```
