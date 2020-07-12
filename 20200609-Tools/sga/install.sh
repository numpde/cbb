#!/usr/bin/env bash

# SGA requires bamtools, libsparsehash-dev, zlib1g-dev

BASE=${PWD}
LOCAL=$BASE/local

[[ -e "sga" ]] || git clone https://github.com/jts/sga

[[ -e "sparsehash" ]] || (
  git clone https://github.com/sparsehash/sparsehash.git; \
  cd $BASE/sparsehash; \
  ./configure --prefix=$LOCAL; \
  make install;
)


cd sga/src

./autogen.sh || exit
./configure --prefix=$LOCAL --with-sparsehash=$LOCAL --with-bamtools=/usr || exit
make
make install
