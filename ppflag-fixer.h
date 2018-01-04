#ifndef _MAIN_H_
#define _MAIN_H_

#include <ctime>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <iostream>
#include <sstream>
#include <vector>

#include "htslib/bgzf.h"
#include "htslib/hfile.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/faidx.h" 
#include "htslib/kstring.h"
#include "htslib/khash_str2int.h"
#include "htslib/synced_bcf_reader.h"
#include "htslib/vcf.h"

// bam flags
#define READ_PAIRED 0x1
#define READ_MAPPED_PROPER_PAIR 0x2
#define READ_UNMAPPED 0x4
#define MATE_UNMAPPED 0x8
#define READ_REVERSE_STRAND 0x10
#define MATE_REVERSE_STRAND 0x20
#define READ_FIRST_IN_PAIR 0x40
#define READ_SECOND_IN_PAIR 0x80
#define READ_SECONDARY_ALIGNMENT 0x100
#define READ_FAILED_QUALITY_CHECKS 0x200
#define READ_PCR_OPTICAL_DUPLICATE 0x400
#define READ_SUPPLEMENTARY 0x800
#define READ_PRIMARY 0x900

using namespace std;

struct arguments {
	vector<char*> args;
	uint64_t max_tlen;
	bool progress;
};

#endif
