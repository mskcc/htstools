#include <algorithm>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <ctype.h>
#include <inttypes.h>
#include <map>
#include <stdint.h>
#include <string>
#include <vector>

#include <argp.h>

#include "dnafrags.h"

using namespace std;

static char args_doc[] = "<sequence file>";

static struct argp_option options[] = {
  {"include-all", 'i', 0, 0,  "Include all reads, even if they fail the vendor's quality check or are duplicates." },
  {"verbose", 'v', 0, 0, "Show detailed messages."},
  {"bin-size", 'b', "SIZE", 0, "Bins reads into bins of SIZE bases. Set to 0 to disable. (default: 50)" },
  {"mapq-min", 'm', "QUALITY", 0, "Discards all reads with a MAPQ less than QUALITY. (default: 50)" },
  {"length-min", 'l', "LENGTH", 0, "Discards all reads with a length less than LENGTH. (default: none)" },
  {"length-max", 'L', "LENGTH", 0, "Discards all reads with a length greater than LENGTH. (default: none)" },
  {"output", 'o', "FILE", 0, "If specified, will output all data in CSV format to the given file path. (default: none)" },
  {"reads-min", 'R', "READS", 0, "Discard all bins with less than READS reads. (default: 0)" },
  {"fill-in", 'f', 0, 0, "Fill in the spots between data points with 0 values."},
  {"rnext-ineq", 'r', 0, 0, "Include reads where RNEXT is not '='."},
  {"progress", 'p', 0, 0, "Keep track of progress through the file. This requires the file to be indexed."},
  { 0 }
};

static error_t parse_opt (int key, char *arg, struct argp_state *state) {
	struct arguments *arguments = (struct arguments*) state->input;

	switch (key) {
	case 'i':
		arguments->include_all = true;
		break;

	case 'v':
		arguments->verbose = true;
		break;

    case 'm':
		arguments->mapq_min = atoi(arg);
		break;

    case 'l':
		arguments->length_min = atoi(arg);
		break;

    case 'L':
		arguments->length_max = atoi(arg);
		break;

    case 'R':
		arguments->reads_min = atoi(arg);
		break;

    case 'o':
		arguments->output = arg;
		break;

	case 'f':
		arguments->fill_in = true;
		break;

	case 'r':
		arguments->rnext_ineq = true;
		break;

	case 'p':
		arguments->progress = true;
		break;

    case 'b':
		arguments->bin_size = atoi(arg);
		break;

	case ARGP_KEY_ARG:
		if (state->arg_num >= 1) {
			// Too many arguments.
			argp_usage (state);
		}
		
		arguments->args[state->arg_num] = arg;
		
		break;
		
	case ARGP_KEY_END:
		if (state->arg_num < 1) {
			// Not enough arguments.
			argp_usage (state);
		}
		break;
		
	default:
		return ARGP_ERR_UNKNOWN;
	}
	return 0;
}

static struct argp argp = { options, parse_opt, args_doc, 0 };

int roundPosition(int numToRound, int bin_size, int bin_size_half) {
    int remainder = numToRound % bin_size;
    if (remainder == 0)
        return numToRound;

	if (remainder < bin_size_half) {
		// then round down
		return numToRound - remainder;
	} else {
		// round up
		return numToRound + bin_size - remainder;
	}
}

double calcMedian(vector<uint64_t> lengths) {
	// get median
	double median = 0.0;
	sort(lengths.begin(), lengths.end());
	uint64_t median_index = (lengths.size() / 2);
	if (lengths.size() == 0) {
		median = 0.0;
	} else if (lengths.size() % 2 == 0) {
		median = (lengths[median_index] + lengths[median_index - 1]) / 2.0;
	} else {
		median = lengths[median_index];
	}
	return median;
}

void flushOut(arguments arguments, FILE * output_file, map<uint64_t, position_info> positions, vector<uint64_t> position_keys, bam_hdr_t *hdr, int32_t current_chromosome) {
	char * chr_name = hdr->target_name[current_chromosome];
	if (!strncmp("chr", chr_name, 3)) {
		chr_name += 3;
	}
	for (uint64_t i = 0; i < position_keys.size(); i++) {
		uint64_t position_key = position_keys[i];
		position_info info = positions[position_key];
		double median = calcMedian(info.lengths);
		if (info.reads < arguments.reads_min) {
			continue; // skip it
		}
		if (arguments.output) {
			if (arguments.fill_in) {
				// check if there are in betweens
				uint64_t lastPosition = 0;
				if (i != 0) {
					lastPosition = position_keys[i - 1];
				}
				uint64_t shouldBeLast = position_keys[i] - 50;
				if (i == 0) {
					shouldBeLast = 0;
				}
				if (shouldBeLast != lastPosition) {
					uint64_t wrotePos = lastPosition + 50;
					while (wrotePos != position_key) {
						fprintf(output_file, "%" PRIu64 ",0,0\n", wrotePos);
						wrotePos += 50;
					}
				}
			}
			// output it to the csv
			fprintf(output_file, "%" PRIu64 ",%" PRIu64 ",%.1f,%s\n", position_key, info.reads, median, chr_name);
		} else {
			printf("%" PRIu64 " in chromosome %s had %" PRIu64 " segment(s), with a median length of %.1f.\n", position_key, chr_name, info.reads, median);
		}
	}
}

int program_main(arguments arguments) {
	clock_t start = clock();
	double duration;

	// try to load sequence file
	hFILE *hfp = hopen(arguments.args[0], "r");
	htsFormat fmt;
	
	if (!hfp) {
		printf("Failed to read sequence file %s.\n", arguments.args[0]);
		return 1;
	}

	hts_detect_format(hfp, &fmt);

	if (arguments.verbose) {
		printf("Detected format: %s\n", hts_format_description(&fmt));
	}
	
	samFile *in = hts_hopen(hfp, arguments.args[0], "r");
	bam_hdr_t *hdr = NULL;
	if (!in) {
		printf("Failed to read file. Is it corrupt?\n");
		return 1;
	}
	hdr = sam_hdr_read(in);
	if (!hdr) {
		printf("Failed to read file header. Is it corrupt?\n");
		return 1;
	}

	hts_idx_t* idx = sam_index_load(in, arguments.args[0]);
	if (arguments.progress && !idx) {
		printf("Failed to read index, and --progress is set.\n");
		return 1;
	}

/*	for (int i = 0; i < hdr->n_targets; ++i) {
		printf("%s, length %d\n", hdr->target_name[i], hdr->target_len[i]);
   	}*/

	if (arguments.output) {
		// try to see if we can open the file
		FILE * test_handle = fopen(arguments.output, "r");
		if (test_handle && strcmp(arguments.output, "/dev/null")) {
			// if we can and it's not /dev/null, complain it exists
			fclose(test_handle);
			printf("Output file %s already exists!\n", arguments.output);
			return 1;
		}
		// don't need to fclose() it here because it should be NULL
	}
	
	FILE * output_file = (arguments.output ? fopen(arguments.output, "w+") : NULL);
	uint64_t nopos_skipped = 0;
	uint64_t rnext_skipped = 0;
	uint64_t mapq_skipped = 0;
	uint64_t len_skipped = 0;
	uint64_t skipped = 0;
	uint64_t total = 0;
    int32_t current_chromosome = -1;
	map<uint64_t, position_info> positions = map<uint64_t, position_info>();
	vector<uint64_t> position_keys = vector<uint64_t>();
	uint64_t position = 0;
	bool has_flushed = false;
	uint64_t index_total = 0;
	uint64_t prog_multiple = 0;
	float last_progress = 0.0;
	int bin_size = arguments.bin_size;
	int bin_size_half = bin_size / 2;
	
	if (arguments.progress) {
		for (int i = 0; i < hdr->n_targets; ++i) {
			uint64_t u, v; 
			hts_idx_get_stat(idx, i, &u, &v);
			index_total += u;
			index_total += v;
		}
		index_total += hts_idx_get_n_no_coor(idx);
		hts_idx_destroy(idx);
		prog_multiple = index_total / 100;
		if (arguments.verbose) {
			printf("Reading %" PRIu64 " read(s).\n", index_total);
		}
	}
	
	if (arguments.output && !output_file) {
		printf("Failed to open output file %s for writing. Do you have proper permissions?\n", arguments.output);
		return 1;
	}

	if (!arguments.output) {
		if (!arguments.silent) {
			printf("Results:\n");
		}
	} else {
		fprintf(output_file, "Midpoint,Reads,Median length,Chromosome\n");
	}
	
	bam1_t *b = bam_init1();
	int ret; 
	while ((ret = sam_read1(in, hdr, b)) >= 0) {
		const bam1_core_t *c = &b->core;
		total++;
		if (!arguments.include_all && ((c->flag & READ_FAILED_QUALITY_CHECKS) || (c->flag & READ_PCR_OPTICAL_DUPLICATE))) {
			skipped++;
			continue;
		}
		if (c->pos == -1 || c->isize == 0) {
			nopos_skipped++;
			continue;
		}
		if (c->qual < arguments.mapq_min) {
			mapq_skipped++;
			continue;
		}
		if (arguments.length_min != -1 && c->isize < arguments.length_min) {
			len_skipped++;
			continue;
		}
		if (arguments.length_max != -1 && c->isize > arguments.length_max) {
			len_skipped++;
			continue;
		}
		if (!arguments.rnext_ineq && c->mtid != c->tid) {
			rnext_skipped++;
			continue;
		}
		if (current_chromosome == -1) {
			// you are the current chromosome!
			current_chromosome = c->tid;
			has_flushed = false;
		}
		if (c->tid != current_chromosome) {
			// you are a different chromosome!
			// write the data out!
			(*arguments.outFunc)(arguments, output_file, positions, position_keys, hdr, current_chromosome);
			// reset structures!
			positions = map<uint64_t, position_info>();
			position_keys = vector<uint64_t>();
			has_flushed = false;
			// and set it
			current_chromosome = c->tid;
		}
		
		uint64_t start_position = c->pos;
		if (c->isize < 1) {//(!(c->flag & READ_FIRST_IN_PAIR)) {
			// you aren't the primary, ignored
			continue;
		}
		uint64_t total_length = c->isize;

		//uint64_t end_position = start_position + total_length;
		uint64_t midpoint = start_position + (total_length / 2);
		//printf("%d\n", total_length);

		//printf("%d -> %d\n", start_position, roundPosition(start_position, 50));
		uint64_t markPos = midpoint;
		if (bin_size != 0) {
			markPos = roundPosition(midpoint, bin_size, bin_size_half);
		}
		if (positions[markPos].reads == 0) {
			// create it
			position_keys.push_back(markPos);
			positions[markPos] = position_info();
			positions[markPos].lengths = vector<uint64_t>();
		}

		positions[markPos].reads++;
		positions[markPos].lengths.push_back(total_length);

		float progress = ((float)total/index_total);
		if (arguments.progress && (int)(last_progress*100) != (int)(progress*100)) {
			int barWidth = 70;
			
			std::cout << "[";
			int pos = barWidth * progress;
			for (int i = 0; i < barWidth; ++i) {
				if (i < pos) std::cout << "=";
				else if (i == pos) std::cout << ">";
				else std::cout << " ";
			}
			std::cout << "] " << int(progress * 100.0) << " %\r";
			std::cout.flush();

			last_progress = progress;
//			printf("%0.0f%% complete (%" PRIu64 " out of %" PRIu64 ")\n", ((float)total/index_total) * 100, total, index_total);
		}
	} 
	bam_destroy1(b);

	duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;

	if (!arguments.silent) {
		printf("Skipped %" PRIu64 " read(s) because of vendor quality check failure.\n", skipped);
		printf("Skipped %" PRIu64 " read(s) because below quality threshold.\n", mapq_skipped);
		if (arguments.length_max > -1) {
			printf("Skipped %" PRIu64 " read(s) because above length threshold.\n", len_skipped);
		}
		printf("Skipped %" PRIu64 " read(s) because of missing POS/TLEN field.\n", nopos_skipped);
		printf("Skipped %" PRIu64 " read(s) because of invalid RNEXT field.\n", rnext_skipped);
		printf("Total %" PRIu64 " read(s) (processed %" PRIu64 ", took %f seconds).\n", total, (total - skipped - mapq_skipped - len_skipped - nopos_skipped), duration);
	}
	if (!has_flushed) {
		(*arguments.outFunc)(arguments, output_file, positions, position_keys, hdr, current_chromosome);
	}
	if (arguments.output) {
		printf("Output results to %s.\n", arguments.output);
		fclose(output_file);
	}
}

int main(int argc, char ** argv) {	
	struct arguments arguments;

	arguments.include_all = false;
	arguments.bin_size = 50;
	arguments.length_min = -1;
	arguments.length_max = -1;
	arguments.reads_min = 0;
	arguments.mapq_min = 50;
	arguments.output = 0;
	arguments.verbose = false;
	arguments.fill_in = false;
	arguments.rnext_ineq = false;
	arguments.progress = false;
	arguments.silent = false;
	arguments.outFunc = &flushOut;
	
	argp_parse (&argp, argc, argv, 0, 0, &arguments);

	return program_main(arguments);
}
