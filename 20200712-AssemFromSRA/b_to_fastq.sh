#!/usr/bin/env bash

# Uses the NCBI SRA Tools
source ../20200609-Tools/sratools/paths.txt || exit

I=$(pwd)/UV/a_sra
O=$(pwd)/UV/b_fastq

cd $O || exit

cd $I || exit
fasterq-dump -O $O ERR3283546.sra
