#include <string.h>  /* strcpy */
#include <stdlib.h>  /* malloc */
#include <stdio.h>   /* printf */
#include <getopt.h>  /* getopt */
#include <stdbool.h> /* Define boolean type */
#include <unistd.h>
#include "htslib/sam.h" /* Use htslib to interact with bam files, imports stdint.h as well */
#include "uthash.h"  /* hash table */
#include "hash.h"    /* Defining the hash tables actually used */
#include "utils.h" /* Show help and create output dir */
#include "sort.h"

#define rdump(...) read_dump(r2l, lout, l2fp, fout, __VA_ARGS__)
#define ddump(...) deduped_dump(r2l, lout, l2fp, fout, __VA_ARGS__)

// Dealing with global vars
char *OUT_PATH = "";
// Length of LEVEL_FALG is HARD-CODED in log_message() that must be updated if more levels
// Are added
char *LEVEL_FLAG[] = {"ERROR", "WARNING", "", "INFO", "", "DEBUG"};
log_level_t OUT_LEVEL = WARNING;
int64_t KEY_SIZE = 512;
int64_t RN_SIZE = 71;
int64_t CB_LENGTH = 21;
int64_t UB_LENGTH = 21;
int64_t chunk_size = 400000; // Approximately 1GB
int64_t MAX_THREADS = 1;

int main(int argc, char *argv[]) {
    // Use a flag to bypass commandline input during development

    int32_t opt;
    char* current_opt;
    int64_t mapq_thres = 0;
    int64_t out_level_raw = 0;
    int64_t mem_scale = 4;
    bool dedup = false, dryrun = false, verbose = false;
    char *bampath = NULL;
    char *metapath = NULL;
    char *oprefix = NULL;
    tag_meta_t *cb_meta = initialize_tag_meta();
    tag_meta_t *ub_meta = initialize_tag_meta();
    strcpy(ub_meta->tag_name, "UB");
    int64_t cb_field = 0;
    int64_t ub_field = 0;
    char bc_tag[3] = "CB";
    char umi_tag[3] = "UB";
    int32_t return_val = 0;

    // Commandline argument processing
    static struct option cl_opts[] = {
            {"file", required_argument, NULL, 'f'},
            {"meta", required_argument, NULL, 'm'},
            {"output", required_argument, NULL, 'o'},
            {"mapq", required_argument, NULL, 'q'},
            {"platform", required_argument, NULL, 'p'},
            {"dedup", no_argument, NULL, 'd'},
            {"cbc-location", required_argument, NULL, 'b'},
            {"cbc-length", required_argument, NULL, 'L'},
            {"umi-location", required_argument, NULL, 'u'},
            {"umi-length", required_argument, NULL, 'l'},
            {"rn-length", required_argument, NULL, 'r'},
            {"mem", required_argument, NULL, 'M'},
            {"threads", required_argument, NULL, '@'},
            {"dry-run", no_argument, NULL, 'n'},
            {"verbose", optional_argument, NULL, 'v'},
            {"help", no_argument, NULL, 'h'}
    };

    log_msg("Parsing commandline flags", DEBUG);
//    ":dq:f:m:o:b:l:u:i:nr:M:v::h"
    while ((opt = getopt_long(argc, argv, ":f:m:o:q:p:db:L:u:l:r:M:@:nv:h", cl_opts, NULL)) != -1) {
        current_opt = argv[optind - 1];
        switch (opt) {
            case 'f':
                bampath = optarg;
                break;
            case 'm':
                metapath = optarg;
                break;
            case 'o':
                oprefix = optarg;
                break;
            case 'q':
                // If optarg is not a number, it will default to 0...
                mapq_thres = strtol(optarg, NULL, 10);
                break;
            case 'p':
                set_CB(cb_meta, optarg);
                set_UB(ub_meta, optarg);
                CB_LENGTH = cb_meta->length;
                UB_LENGTH = ub_meta->length;
                break;
            case 'd':
                dedup = true;
                break;
            case 'b':
                cb_field = strtol(optarg, NULL, 10);
                if (cb_field == 0) {
                    strncpy(cb_meta->tag_name, optarg, 2);
                } else {
                    cb_meta->location = READ_NAME;
                    cb_meta->field = cb_field;
                }
                break;
            case 'L':
                CB_LENGTH = strtol(optarg, NULL, 10) + 1;
                if (CB_LENGTH != 0) {
                    cb_meta->length = CB_LENGTH;
                } else {
                    log_msg("Cell barcode length must be larger than 0", ERROR);
                    goto error_out_and_free;
                }
                break;
            case 'u':
                ub_field = strtol(optarg, NULL, 10);
                if (ub_field == 0) {
                    strncpy(ub_meta->tag_name, optarg, 2);
                } else {
                    ub_meta->location = READ_NAME;
                    ub_meta->field = ub_field;
                }
                break;
            case 'l':
                UB_LENGTH = strtol(optarg, NULL, 10) + 1;
                if (UB_LENGTH > 0) {
                    ub_meta->length = UB_LENGTH;
                } else {
                    log_msg("UMI length must be larger than 0", ERROR);
                    goto error_out_and_free;
                }
                break;
            case 'r':
                RN_SIZE = strtol(optarg, NULL, 10) + 1;
                if (RN_SIZE <= 1) {
                    log_msg("Read name length must be larger than 0", ERROR);
                    goto error_out_and_free;
                }
                break;
            case 'M':
                mem_scale = strtol(optarg, NULL, 10);
                if (mem_scale < 1) {
                    log_msg("Memory limit (-M/--mem) must be an integer and >= 1", ERROR);
                    goto error_out_and_free;
                }
                break;
            case '@':
                MAX_THREADS = strtol(optarg, NULL, 10);
                if (MAX_THREADS == 0) {
                    MAX_THREADS = 1;
                }
                break;
            case 'n':
                dryrun = true;
                break;
            case 'v':
                // Manual optional results in possible consumption of the next flag and has to be dealt
                // with
                if (optarg[0] == '-') {
                    optind = optind - 1;
                    break;
                }
                out_level_raw = strtol(optarg, NULL, 10);
                verbose = true;
                if (out_level_raw == 0) {
                    OUT_LEVEL = INFO;
                } else if (out_level_raw > 0 && out_level_raw < 5) {
                    OUT_LEVEL = out_level_raw;
                } else {
                    OUT_LEVEL = DEBUG;
                }
                break;
            case 'h':
                show_usage();
                return 0;
            case ':':
                if (strcmp(current_opt, "-v") == 0) {
                    verbose = true;
                    OUT_LEVEL = INFO;
                    break;
                }
                log_msg("%s must be provided with a value", ERROR, current_opt);
                log_msg("Please see \"scbamsplit --help\" for details", ERROR);
                goto error_out_and_free;
            default:
                log_msg("%s is not a valid option", ERROR, current_opt);
                log_msg("Please see \"scbamsplit --help\" for details", ERROR);
                error_out_and_free:
                    destroy_tag_meta(cb_meta);
                    destroy_tag_meta(ub_meta);
                return 1;
        }
    }

    // Error-out if missing arguments
    if (bampath == NULL || metapath == NULL) {
        log_msg("Error: Missing required arguments", ERROR);
        log_msg("Please see \"scbamsplit --help\" for details", ERROR);
        return 1;
    }

    // If there's no output prefix set, export files in the work directory
    if (oprefix == NULL) {
        oprefix = "./";
    }

    // Set chunk size by mem estimation
    chunk_size = chunk_size * mem_scale / MAX_THREADS;

    // If the output prefix does not end with /, add it.
    uint16_t oplen = strlen(oprefix) - 1;
    if (oprefix[oplen] != '/') {
        oprefix = strcat(oprefix, "/");
    }

    if (mapq_thres > 254) {
        fprintf(stderr, "Please note that the maximal value of MAPQ is 255.\n");
        fprintf(stderr, "There is no read that would have MAPQ **ABOVE** the current threshold (%lld) and be kept.\n",
                mapq_thres);
        return 1;
    }

    if (verbose || dryrun) {
        fprintf(stderr, "- Run condition:\n");
        fprintf(stderr, "\tInput bam: %s\n", bampath);
        fprintf(stderr, "\tInput metadata: %s\n", metapath);
        fprintf(stderr, "\tMAPQ threshold: %lld\n", mapq_thres);
        fprintf(stderr, "\tRead name length: %lldmer\n", RN_SIZE - 1);
        fprintf(stderr, "\tOutput prefix: %s\n", oprefix);
        fprintf(stderr, "\tMemory usage is estimated to be: %lldGB\n", mem_scale);
        fprintf(stderr, "\tLogging level is %d\n", OUT_LEVEL);
        print_tag_meta(cb_meta, "Cell barcode");
        print_tag_meta(ub_meta, "UMI");
        if (dedup) {
            fprintf(stderr, "\tRunning **with** deduplication.\n\n");
        } else {
            fprintf(stderr, "\tRunning **without** deduplication.\n\n");
        }
    }

    if (access(bampath, F_OK) != 0) {
        destroy_tag_meta(cb_meta);
        destroy_tag_meta(ub_meta);
        log_msg("%s not found", ERROR, bampath);
        return 1;
    }


    if (dryrun) {
        fprintf(stderr, "\t==========================================================\n");
        fprintf(stderr, "\t= This is a dry-run (-n/--dry-run). Nothing is executed. =\n");
        fprintf(stderr, "\t==========================================================\n");
        destroy_tag_meta(cb_meta);
        destroy_tag_meta(ub_meta);
        return 0;
    }

    // Create output folder if it does not exist
    // If it exists, ask the user for confirmation to prevent unexpected overwriting
    log_msg("Creating output directory", INFO);
    int mkdir_status;
    if (1 == (mkdir_status = create_directory(oprefix))) {
        log_msg("Exiting because the user declined overwrite", INFO);
        log_msg("Please provide a new path for the output directory", WARNING );
        return 0;
    } else if (-1 == mkdir_status) {
        return 1;
    }

    ////////// bam related //////////////////////////////////////////////
    // Open input bam file from CellRanger
    // Remember to close file handle!
    log_msg("Reading input BAM file: %s", INFO, bampath);
    samFile *fp = sam_open(bampath, "r");


    // Extract header
    log_msg("Reading SAM header", DEBUG);
    sam_hdr_t *header = sam_hdr_read(fp);
    if (dedup) sam_hdr_change_HD(header, "SO", "scbamsplit");


    // Iterating variables for reads in the bam file
    log_msg("Create temporary read", DEBUG);
    bam1_t *read = bam_init1();

    // Prepare a read-tag-to-label hash table from a metadata table
    log_msg("Loading barcode-cluster information from metadata: %s", INFO, metapath);
    rt2label *r2l = hash_readtag((char *) metapath);

    if (r2l == NULL) {
        log_msg("Failed to hash the metadata.", ERROR);
        return 1;
    }

    // Prepare a label-to-file-handle hash table from the above
    log_msg("Preparing output BAM files", INFO);
    label2fp *l2fp = NULL;
    l2fp = hash_labels(r2l, oprefix, header);

    // Iterate through the rt's and write to corresponding file handles.
    // Iterate through reads from input bam
    rt2label *lout;
    label2fp *fout;
    int32_t read_stat;


    if (!dedup) {
        char *current_CB;
        char *this_CB;
        char *this_UB;
        current_CB = (char *) calloc(CB_LENGTH, sizeof(char));
        this_CB = (char *) calloc(CB_LENGTH, sizeof(char));
        this_UB = (char *) calloc(UB_LENGTH, sizeof(char));
        while (0 <= (read_stat = sam_read1(fp, header, read))) {
            // Get read metadata
            int8_t ub_stat = get_UB(read, ub_meta, this_UB);
            int8_t cb_stat = get_CB(read, cb_meta, this_CB);
            int16_t mapq = (int16_t) read->core.qual;

            if (-1 == cb_stat || -1 == ub_stat || mapq < mapq_thres) {
                // Ignore reads without CB and UMI for consistency
                continue;
            }
            // Exporting process
            int8_t rdump_stat = rdump(this_CB, header, read);
            if (0 != rdump_stat) {
                log_msg("Fail to write sorted reads to individual BAM file (%s)", ERROR, lout->label);
            }
        }
        // No-dedup split done
        free(current_CB);
        free(this_CB);
        free(this_UB);
    } else {
        // Deduplication-specific code
        log_msg("Processing %lld reads per chunk", INFO, chunk_size);

        // Allocate heap memory for reads to sort
        log_msg("Preparing read chunks for sorting", DEBUG);

        char* tmpdir = process_bam(fp, header, chunk_size, oprefix, mapq_thres, cb_meta, ub_meta);
        if (strcmp(tmpdir, "1") == 0) {
            return_val = 1;
            goto early_exit;
        }

        // Done processing
        log_msg("Completed sorting all chunks", INFO);

        char * sorted_path = merge_bams(tmpdir);

        // Deduped-split
        log_msg("Open sorted file (%s) to split", INFO, sorted_path);

        int8_t dump_stat = ddump(tmpdir, sorted_path, read, bc_tag, umi_tag, cb_meta, ub_meta);
        if (1 == dump_stat) {
            log_msg("Please check the output folder to remove remaining temporary folder", WARNING);
        }
    }

    // Release and exit
early_exit:
    sam_close(fp);
    bam_hdr_destroy(header);
    label2fp *qs, *qtmp;
    HASH_ITER(hh, l2fp, qs, qtmp) {
        sam_close(qs->fp);
        HASH_DEL(l2fp, qs);
        free(qs);
    }

    // free the hash table contents
    rt2label *s, *tmp;
    HASH_ITER(hh, r2l, s, tmp) {
        HASH_DEL(r2l, s);
        free(s);
    }

    // Free temporary read
    bam_destroy1(read);
    destroy_tag_meta(cb_meta);
    destroy_tag_meta(ub_meta);

    return return_val;
}
