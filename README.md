# Welcome to collapse RepeatMasker annotation page

### About
The script `collapse_RM_annotation.py` collapses TE annotations from RepeatMasker (i.e. genome.fa.out file). The TE annotations from RepeatMasker are frequently overlap each other. If so, when counting NGS reads mapping to TEs, reads mapping to two or more TE annotations may not be counted. In such case, it is preferable to use non-overlapping TE annotations. This script makes non-overlapping TE annotations by removing TEs with lower bit score. Processing of one file (e.g. hg38, mm10) requires ~3 min and ~100 MB RAM. This script only uses single thread.  
  
### Requirement
- Python 3.6 or later.
- Python modules (all should be built-in): os, sys, gzip, datetime, collections, argparse, errno
  
### How to use
If you do not have `genome.fa.out`, you need to generate it by RepeatMasker.
The command below will output `./RM_out/GRCm38.p6.genome.fa.out` file, which contains repeat annotations.
```
# this is the example of GRCm38
RepeatMasker \
GRCm38.p6.genome.fa \
-species 10090 \
-s -no_is \
-dir ./RM_out \
-pa 8
```
  
The script `collapse_RM_annotation.py` will remove overlapping repeat annotations based on bit score.
The command below will generate two files: `GRCm38.p6.genome.fa.out.collapsed.gtf.gz` and `GRCm38.p6.genome.fa.out.collapsed.bed.gz`.
The `.gtf.gz` file can be used for counting reads mapping to repeat anotations (e.g., STAR, featureCount).  
  
The two flags below are required. Please specify an input `.fa.out` file with the `-i` flag, and an output file basename with the `-o` flag.
  
- `-i [input.fa.out]`  
- `-o [output basename]`
  
```
python collapse_RM_annotation.py \
-i ./RM_out/GRCm38.p6.genome.fa.out \
-o GRCm38.p6.genome.fa.out.collapsed
```
  
### Output files
It will generate two files: `.gtf.gz` and `.bed.gz`.  

- `.gtf.gz`  
Each annotation in this GTF file will be composed of one gene, one transcript, and one exon.  
`gene_id` in the attribution field will be `RM_n`.`repeat`.`repeat_class`.`original_id`, where `n` is the serial numbering, `repeat` is the 10th column of the input `.fa.out` file, `repeat_class` is the 11 th column of the input `.fa.out` file, and `original_id` is the 15th column of the input `.fa.out` file.  

- `.bed.gz`  
The 4th column (name field) in the bed file will be `RM_n`.`repeat`.`repeat_class`.`original_id`, which is the same naming convention as the `gene_id` in `.gtf.gz` file.  
The 5th column (score filed) will have bit score (1st column of the input `.fa.out` file).  
  
### Options that affect output results
- `-keep_simple_repeat`  
By default, the script will remove "Simple_repeat" and "Low_complexity" from output files.
If you want to include such annotations, please add the `-keep_simple_repeat` option.  
```
python collapse_RM_annotation.py \
-i ./RM_out/GRCm38.p6.genome.fa.out \
-o GRCm38.p6.genome.fa.out.collapsed \
-keep_simple_repeat
```
  
- `-gap [int]` (default = 0)  
Sometimes, one TE annotation will be split into several fragments in the `.fa.out` file.
In such case, This script can connect those fragments if the fragments were close each other.
By default, if the same two or more repeats are nest to each other without gap, it connects those and report one repeat.
On the other hand, if you specify the `-gap 5` option, it will connect the same two or more repeats with gap distance 5 bases or less.  
```
python collapse_RM_annotation.py \
-i ./RM_out/GRCm38.p6.genome.fa.out \
-o GRCm38.p6.genome.fa.out.collapsed \
-gap 5
```
  
### Further options that does not affect results
- `-quiet`  
- `-version`  
- `-h`  
  
Please see further information by the `-h` option.  
```
python collapse_RM_annotation.py -h
```


