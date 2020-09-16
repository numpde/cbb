#!/usr/bin/env bash

# https://github.com/ablab/spades/blob/spades_3.14.1/README.md

renice -n 19 $$

source ../20200609-Tools/spades/paths.txt

spades.py \
  -1 UV/b_fastq/ERR3283546.sra_1.fastq \
  -2 UV/b_fastq/ERR3283546.sra_2.fastq \
  --careful \
  -o UV/c_spades
