#!/bin/bash
# Accuracy simulation from the tidk paper (Supplementary Table S3), run
# against one or more tidk binaries on the same simulated fastas.
#
# usage: ./run.bash name1=/path/to/tidk [name2=/path/to/other/tidk ...] > results.tsv
# e.g.   ./run.bash old=../old/tidk new=../../target/release/tidk
#
# Simulated fastas and raw explore output go in sim/ and out/ (gitignored).
set -euo pipefail

if [ $# -eq 0 ]; then
  echo "usage: $0 name=/path/to/tidk [name=/path/to/tidk ...]" >&2
  exit 1
fi

cd "$(dirname "$0")"
mkdir -p sim out

sequence_lengths=(600 12000 30000)
error_rates=(0.0 0.001 0.01 0.015 0.02 0.05 0.1)
num_sequences=1000
motif="TTAGGG"

echo -e "length\terror\tversion\tmost_abundant\tmost_abundant_count\tmost_abundant_6nt\tAACCCT_count\tAACCCT_rank\ttotal_units"
for seq_len in "${sequence_lengths[@]}"; do
  for err_rate in "${error_rates[@]}"; do
    fa=sim/len${seq_len}_err${err_rate}.fasta
    python3 random_telomere_generator.py $motif $err_rate $seq_len $num_sequences $fa

    if (( seq_len < 1000 )); then threshold=10; else threshold=100; fi

    for spec in "$@"; do
      name=${spec%%=*}
      bin=${spec#*=}
      tsv=out/len${seq_len}_err${err_rate}_${name}.tsv
      $bin explore -t $threshold --minimum 5 --maximum 12 --distance 0.5 $fa > $tsv 2>/dev/null

      body=$(tail -n +2 $tsv)
      top=$(echo "$body" | head -1 | cut -f1)
      top_count=$(echo "$body" | head -1 | cut -f2)
      six=$(echo "$body" | grep -m1 -E '^\S{6}\s' | cut -f1 || true)
      aaccct=$(echo "$body" | grep -E '^AACCCT\s' | cut -f2 || true)
      rank=$(echo "$body" | grep -n -E '^AACCCT\s' | cut -d: -f1 || true)
      n=$(echo "$body" | grep -c . || true)
      echo -e "$seq_len\t$err_rate\t$name\t${top:-N}\t${top_count:-0}\t${six:-N}\t${aaccct:-0}\t${rank:-NA}\t$n"
    done
  done
done
