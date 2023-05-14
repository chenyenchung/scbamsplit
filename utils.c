//
// Created by Yen-Chung Chen on 2/24/23.
//
#include <stdio.h>
#include <string.h>
#include "utils.h"
#include "htslib/sam.h"
#include "shared_const.h"

int show_usage(const char* type) {
    if (strcmp(type, "unknown") == 0) {
        fprintf(stderr, "Unknown options are provided. Please see below for available options.\n\n");
    }
    if (strcmp(type, "missing") == 0) {
        fprintf(stderr, "A value must be provided for the following options:\n");
        fprintf(stderr, "-f/--file, -m/--meta, -q/--mapq, -o/--output, -t/--filter-tag, and -u/--umi-tag.\n\n");
    }
    fprintf(stderr, "Program: %s (Parallel BAM file subset by read tags)\n", APP_NAME);
    fprintf(stderr, "Version: %s (Dependent on htslib v%s)\n\n", VERSION, hts_version());
    fprintf(stderr, "Usage: scbamsplit [-f path] [-m path] [-o path] [-q MAPQ] [-t read-tag] [-u read-tag] [-d] [-h] [-n] [-v] \n\n");
    fprintf(stderr, "\t-f/--file: the path for input bam file\n");
    fprintf(stderr, "\t-m/--meta: the path for input metadata (an unquoted two-column csv with column names)\n");
    fprintf(stderr, "\t-o/--output: the path to export bam files to (default: ./)\n");
    fprintf(stderr, "\t-q/--mapq: Minimal MAPQ threshold for output (default: 0)\n");
    fprintf(stderr, "\t-t/--filter-tag: The read tag you want to filter against (default: CB)\n");
    fprintf(stderr, "\t-u/--umi-tag: The read tag containing UMIs (default: UB)\n");
    fprintf(stderr, "\t-d/--dedup: Remove duplicated reads with the same cell barcode/UMI combination\n");
    fprintf(stderr, "\t-n/--dry-run: Only print out parameters\n");
    fprintf(stderr, "\t-v/--verbose: Print out parameters\n");
    fprintf(stderr, "\t-h/--help: Show this documentation\n");
    return 0;
}
