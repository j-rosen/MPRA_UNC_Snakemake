#!/bin/bash

# Count number of arguments
n=$( echo $# )

# Throw error if less than two arguments
if [[ $n -lt 3 ]]; then
  echo "At least two arguments required!"
  exit 1
fi

# Join first two files
# -e sets the unmatched output
# -a1 and -a2 necessary for outputting unmatched fields
# -o auto necessary to print the 0 from the -e flag

join -e 0 -a1 -a2 -o auto \
  <( sort -k2,3 $1 | awk '{print $2","$3"\t"$1}' ) \
  <( sort -k2,3 $2 | awk '{print $2","$3"\t"$1}' ) \
  > ${!n}


# If more than two files, join them iteratively
if [[ $n -gt 3 ]]; then
  for i in `seq 3 $(( $n - 1 ))`
  do
    join -e 0 -a1 -a2 -o auto \
      ${!n} \
      <( sort -k2,3 ${!i} | awk '{print $2","$3"\t"$1}' ) \
      > ${!n}.temp

    mv ${!n}.temp ${!n}
  done
fi

# Process results
sed 's/,/ /g' ${!n} > ${!n}.temp
mv ${!n}.temp ${!n}


