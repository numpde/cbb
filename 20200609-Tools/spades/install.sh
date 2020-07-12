#!/usr/bin/env bash

# http://cab.spbu.ru/software/spades/

wget -nc http://cab.spbu.ru/files/release3.14.1/SPAdes-3.14.1-Linux.tar.gz -O spades.tar.gz
tar -vxzf spades.tar.gz && rm spades.tar.gz

echo "# DO:" > paths.txt
echo "export PATH=\"\$PATH:${PWD}/$(find */bin | sort | head -n 1)\"" >> paths.txt

cat paths.txt
