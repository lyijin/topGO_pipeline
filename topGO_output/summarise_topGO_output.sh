#!/bin/bash

# remove previous summary files
rm -f summary_*.txt

# generate new summary files
# two criteria to pass:
# 1. At least 5 terms in the universe ($4)
# 2. P value < 0.05 OR if it contains the character '<' (because R outputs stuff like '< 1e-30') ($7)
for a in bp*.txt; do
    for b in $a ${a/bp/cc} ${a/bp/mf}; do
        echo -- $b -- >> ${a/bp_/summary_}
        awk -F $'\t' '{if ($4 >= 5 && ($7 < 0.05 || $7 ~ /^</)) print}' $b >> ${a/bp_/summary_}
    done
done

# prettify output
for a in summary_*.txt; do
    sed 's/GO\.ID/\tGO.ID/' $a | sed 's/-- /\n-- /' | sed 1d > tmp && mv -f tmp $a
done