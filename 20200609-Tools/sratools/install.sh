#!/usr/bin/env bash

# Instructions from
# https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit

wget -nc --output-document sratoolkit.tar.gz http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -vxzf sratoolkit.tar.gz && rm sratoolkit.tar.gz

echo "# DO:" > paths.txt
echo "export PATH=\"\$PATH:${PWD}/$(find sra*/bin | sort | head -n 1)\"" >> paths.txt

cat paths.txt
