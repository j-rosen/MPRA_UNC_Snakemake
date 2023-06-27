#!/bin/bash

fastq=$1
sample=$2
fasta=$3

# Find exact barcode matches in fastq files
export LC_ALL=C # Speed up sort
awk 'NR==FNR {a[$2] = $1; next} {if (FNR%4 == 2 && substr($1, 1, 20) in a) print a[substr($1, 1, 20)]"\t"substr($1, 1, 20)}' \
    <( sed 's/^>//g' ${fasta} | paste - - ) \
    <( zcat ${fastq} ) | \
    sort | uniq -c \
    > ${sample}.reads.txt

