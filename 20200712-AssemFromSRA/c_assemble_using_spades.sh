#!/usr/bin/env bash

# https://github.com/ablab/spades/blob/spades_3.14.1/README.md

source ../20200609-Tools/spades/paths.txt

spades.py \
  -1 download/ERR3283546.sra_1.fastq \
  -2 download/ERR3283546.sra_2.fastq \
  --careful \
  -o spades_output
