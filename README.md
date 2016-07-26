# snp-pileup
This application will, given a VCF file containing SNP locations, output for each SNP the counts of the reference nucleotide, alternative nucleotide, errors, and deletions.

# Installation
First, HTSlib must be installed on your system. To do that, [download](http://www.htslib.org/download/) it and follow the "Building and installing" instructions on that page.

Then, download this code, extract it, `cd` to where you extracted it, and run the following:
```sudo ldconfig # only needs to be run the first time
./compile.sh
sudo ./install.sh```

# Usage
`snp-pileup <vcf file> <output file> <sequence files...>`
Usage of snp-pileup requires a VCF file and one (or multiple) sequence files containing DNA. The sequence files should be in the BAM format, and both the VCF and all sequence files must be sorted.

## Parameters
Here is a list of all parameters `snp-pileup` accepts and information about what they do.

* `-A, --count-orphans`
This option will enable counting reads that are not properly paired.
* `-d, --max-depth=DEPTH`
This option will set the maximum number of reads per position&mdash;after `DEPTH` many reads, no more reads will be counted. Default is 4000.
* `-p, --progress`
This option will display a progress bar displaying the percentage of SNPs completed. This option requires the program to run through the VCF file an extra time (to count it); as such it is not recommended for use with large VCF files.
* `-P, --pseudo-snps=MULTIPLE"
This option will insert a record every `MULTIPLE` positions that will count how many bases are at that position, even if there is no SNP already there.
* `-q, --min-map-quality=QUALITY`
This option sets a minimum threshold for the mapping quality of a read. Reads with a mapping quality below `QUALITY` will be ignored. Default is 0.
* `-Q, --min-base-quality=QUALITY`
This option sets a minimum threshold for the quality of a nucleotide. Bases with a quality below `QUALITY` will be ignored. Default is 0.
* `-r, --min-read-counts=READS`
This option will the minimum read count for a SNP to be written to the output file. `READS` is a comma separated list of minimum read counts for each file. For example, setting `READS` to `25,20,10,0` would make a SNP require 25 reads in file 1, 20 reads in file 2, 10 reads in file 3, and 0 reads in file 4&mdash;SNPs that don't meet these requirements would be ignored. Default is 0 for all files.
* `-v, --verbose`
This option enables detailed messages.
* `-x, --ignore-overlaps`
By default, snp-pileup will try to detect where paired reads overlap and count those overlaps as a single read, instead of two different reads that are from the same fragment. This option will disable that detection.

## Limitations
SNPs where there are multiple nucleotides changing will be ignored, and all minimum thresholds (except for the minimum read count) apply equally to all files&mdash;there is no way to set them on a per-file basis.