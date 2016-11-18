#ifndef _MAIN_H_
#define _MAIN_H_

#include <map>
#include <vector>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include "htslib/hfile.h"
#include "htslib/hts.h"
#include "htslib/sam.h"

using namespace std;

// flags
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

struct position_info {
	uint64_t reads;
	vector<uint64_t> lengths;
};

struct chromosome_info {
	char * name;
	map<uint64_t, position_info> positions;
	vector<uint64_t> position_keys;
};

struct arguments {
	char *args[2];
	bool include_all;
	bool verbose;
	int mapq_min;
	int length_min;
	int length_max;
	int reads_min;
	bool fill_in;
	bool rnext_ineq;
	bool progress;
	int bin_size;
	bool silent;
	char * output;
	void (*outFunc)(arguments, FILE *, map<uint64_t, position_info>, vector<uint64_t>, bam_hdr_t *, int32_t);
};

void flushOut(arguments arguments, FILE * output_file, map<uint64_t, position_info> positions, vector<uint64_t> position_keys, bam_hdr_t *hdr, int32_t current_chromosome);
double calcMedian(vector<uint64_t> lengths);
int program_main(arguments arguments);

#endif
