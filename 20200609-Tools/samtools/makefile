# Installing samtools

ROOT=$(shell pwd)

HTSLIB=${ROOT}/htslib
SAMTOOLS=${ROOT}/samtools
LOCAL=${ROOT}/local
PATHSFILE=${ROOT}/paths.txt

default:
	@echo You can install samtools with
	@echo 
	@echo sudo apt install samtools
	@echo 
	@echo or from source:
	@echo 
	@echo To install samtools to ${LOCAL}, do: cat makefile
	@echo Check out the file ${PATHSFILE} after installation

clone: force
	git clone https://github.com/samtools/htslib ${HTSLIB}
	git clone https://github.com/samtools/samtools ${SAMTOOLS}


dependencies:
	sudo apt install zlib1g-dev
	sudo apt install libbz2-1.0 libbz2-dev
	sudo apt install liblzma-dev
	sudo apt-get install libcurl4-openssl-dev
	sudo apt-get install autoconf autogen

samtools: force samtools_install samtools_paths;

samtools_install:
	cd ${SAMTOOLS} && \
	pwd && \
	autoheader && \
	autoconf -Wno-syntax && \
	./configure --without-curses --prefix=${LOCAL} && \
	make && \
	make install;
	
samtools_paths:
	@echo "# DO:" > ${PATHSFILE}
	@echo "export PATH=\"\$$PATH:${LOCAL}/bin\"" >> ${PATHSFILE}
	@echo "export MANPATH=\"\$$(manpath -g):${LOCAL}/share/man\"" >> ${PATHSFILE}
	

force: ;

