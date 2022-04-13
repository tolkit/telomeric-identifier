#!/usr/bin/env bash

# parse this data to Rust
python3 build_clades.py > clades.txt

# remove in place the variables
sed '/\/\/ automated input start/,/\/\/ automated input end/{//!d}' ../src/clades.rs > ./temp.rs

# insert the new ones
# get line first
LINE=$(grep -n '\/\/ automated input start' ./temp.rs | cut -d ":" -f 1)
# make new file
{ head -n $LINE ./temp.rs; cat clades.txt; tail -n +$(($LINE+1)) ./temp.rs; } > temp2.rs

rm ../src/clades.rs
mv temp2.rs ../src/clades.rs

# need to format clades too.
rustfmt ../src/clades.rs