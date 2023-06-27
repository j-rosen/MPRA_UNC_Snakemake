#!/bin/bash

# Assign barcodes to designed sequences

design_file=$1
bc_file=$2
oligo_file=$3
sample=$4

# Combine barcodes and oligos
paste ${bc_file} ${oligo_file} | \
awk 'NR%4==2 && NF==2 && length($1)==20' \
> ${sample}.bc.oligo.temp

# Format design file
sed 's/^>//g' ${design_file} | \
paste - - \
> ${sample}.design.temp

# Look up exact matches in design file
export LC_ALL=C # speed up sort

awk 'NR==FNR {a[$2] = $1; next} {if ($2 in a) print a[$2]"\t"$1; else print "NA\t"$1}' \
  ${sample}.design.temp ${sample}.bc.oligo.temp | \
sort -S 7G --parallel=24 | \
uniq -c > ${sample}.sequences.barcodes.txt

rm ${sample}.bc.oligo.temp ${sample}.design.temp



