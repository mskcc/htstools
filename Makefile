# This needs to point to wherever HTSLIB is installed
HTSDIR=/opt/common/CentOS_6-dev/htslib/v1.3.2

# This needs to be a recent version of g++ (4.9+ works)
GXX=/opt/common/CentOS_6/gcc/gcc-4.9.3/bin/g++

#########################################################

CXX=$(GXX) -std=c++11 -I$(HTSDIR) -L$(HTSDIR) -lhts

all: snp-pileup dnafrags ppflag-fixer

snp-pileup: snp-pileup.cpp

#dnafrags: dnafrags.cpp

ppflag-fixer: ppflag-fixer.cpp

clean:
	rm -f snp-pileup dnafrags ppflag-fixer

