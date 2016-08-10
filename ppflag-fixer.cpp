#include "ppflag-fixer.h"

#include <argp.h>

static char args_doc[] = "<input file> <output file>";

static struct argp_option options[] = {
	{"max-tlen", 'm', "LENGTH", 0, "Sets a maximum bound of LENGTH on all fragments; any greater and they won't be marked as proper pair."},
	{"progress", 'p', 0, 0, "Keep track of progress through the file. This requires the file to be indexed."},
	{ 0 }
};

static error_t parse_opt (int key, char *arg, struct argp_state *state) {
	struct arguments *arguments = (struct arguments*) state->input;

	switch (key) {
    case 'm':
		arguments->max_tlen = atoi(arg);
		break;
		
	case 'p':
		arguments->progress = true;
		break;

	case ARGP_KEY_ARG:
		if (state->arg_num > 2) {
			// Too many arguments.
			argp_usage (state);
		}
		arguments->args.push_back(arg);
		break;

	case ARGP_KEY_END:
		if (state->arg_num < 2) {
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


int program_main(arguments arguments) {
	clock_t start = clock();

	// try to load sequence file
	hFILE *in_hfp = hopen(arguments.args[0], "r");
	
	if (!in_hfp) {
		printf("Failed to read sequence file %s.\n", arguments.args[0]);
		return 1;
	}
	
	samFile *in = hts_hopen(in_hfp, arguments.args[0], "r");
	if (!in) {
		printf("Failed to read file. Is it corrupt?\n");
		return 1;
	}

	// try to open output file
	hFILE *out_hfp = hopen(arguments.args[1], "w");
	
	if (!out_hfp) {
		printf("Failed to open output file %s.\n", arguments.args[0]);
		return 1;
	}
	
	samFile *out = hts_hopen(out_hfp, arguments.args[1], "wbz");
	if (!out) {
		printf("Failed to open output file. This shouldn't happen?\n");
		return 1;
	}

	// write out header
	bam_hdr_t *hdr = sam_hdr_read(in);
	if (!hdr) {
		printf("Failed to read header of input file. Is it corrupt?\n");
		return 1;
	}
	if (sam_hdr_write(out, hdr) != 0) {
		printf("Failed to write header to output file.\n");
		return 1;
	}

	// read index totals
	hts_idx_t* idx = sam_index_load(in, arguments.args[0]);
	if (arguments.progress && !idx) {
		printf("Failed to read index, and --progress is set.\n");
		return 1;
	}

	// loop over everything
	bam1_t *b = bam_init1();
	int r;
	uint64_t total = 0;
	uint64_t edited = 0;
	uint64_t index_total = 0;
	float last_progress = 0.0;

	if (arguments.progress) {
		for (int i = 0; i < hdr->n_targets; ++i) {
			uint64_t u, v; 
			hts_idx_get_stat(idx, i, &u, &v);
			index_total += u;
			index_total += v;
		}
		index_total += hts_idx_get_n_no_coor(idx);
		hts_idx_destroy(idx);
		printf("Reading %" PRIu64 "read(s).\n", index_total);
	}
	
	while ((r = sam_read1(in, hdr, b)) >= 0) {
		total++;

		bam1_core_t *c = &b->core;

		// edit the flags
		if ((arguments.max_tlen == -1 || abs(c->isize) <= arguments.max_tlen) && // if there's no max tlen or we are under the max
			(c->mtid == c->tid) && // check if rnext is = 
			(c->flag & READ_PAIRED) && 
			(!(c->flag & READ_UNMAPPED)) && 
			(!(c->flag & MATE_UNMAPPED)) && 
			(!(c->flag & READ_SUPPLEMENTARY)) && 
			(!(c->flag & READ_SECONDARY_ALIGNMENT)) && 
			(!(c->flag & READ_MAPPED_PROPER_PAIR)) &&
			((c->isize < 0 && (c->flag & READ_REVERSE_STRAND)) || //tlen is < 0, this is the mate, should be reverse.
			 (c->isize > 0 && (c->flag & MATE_REVERSE_STRAND)))
			) {
			c->flag |= READ_MAPPED_PROPER_PAIR;
			edited++;
		}
		
		if (sam_write1(out, hdr, b) < 0) {
			printf("Failed to write record, aborting.\n");
			return 1;
		}

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
		}
	}
	if (r < -1) {
		printf("Truncated file!\n");
		return 1;
	}
	bam_destroy1(b);

	hts_close(out);
	/*if (!hclose(out_hfp)) {
		// I don't really care if this fails
		// This is just here so g++ doesn't get mad
		}*/

	printf("%" PRIu64 " reads were changed.\n", edited);
	
	double duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
	printf("Finished in %f seconds.\n", duration);

	return 0;
}

int main(int argc, char ** argv) {
	struct arguments arguments;

	arguments.args = vector<char *>();
	arguments.max_tlen = -1;
	arguments.progress = false;

	argp_parse (&argp, argc, argv, 0, 0, &arguments);

	return program_main(arguments);
}
