#include <string.h>  /* strcpy */
#include <stdlib.h>  /* malloc */
#include <stdio.h>   /* printf */
#include <getopt.h>  /* getopt */
#include <stdbool.h> /* Define boolean type */
#include "htslib/sam.h" /* Use htslib to interact with bam files, imports stdint.h as well */
#include "uthash.h"  /* hash table */
#include "hash.h"    /* Defining the hash tables actually used */
#include "utils.h" /* Show help and create output dir */
#include "sort.h"

// Dealing with global vars
char *OUT_PATH = "";
log_level_t OUT_LEVEL = DEBUG;
uint16_t KEY_SIZE = 512;
uint16_t RN_SIZE = 48;
uint8_t CB_LENGTH = 20;
uint8_t UB_LENGTH = 10;

bool dev = true;

int main(int argc, char *argv[]) {

    // Use a flag to bypass commandline input during development

    int opt;
    uint16_t mapq = 0;
    char *mapqstr = NULL;
    bool dedup = false, dryrun = false, verbose = false;
    char *bampath = NULL;
    char *metapath = NULL;
    char *oprefix = NULL;
    char filter_tag[3] = "CB";
    char umi_tag[3] = "UB";

    // Commandline argument processing
    static struct option cl_opts[] = {
            {"dedup", no_argument, NULL, 'd'},
            {"mapq", required_argument, NULL, 'q'},
            {"file", required_argument, NULL, 'f'},
            {"meta", required_argument, NULL, 'm'},
            {"output", required_argument, NULL, 'o'},
            {"filter-tag", required_argument, NULL, 't'},
            {"umi-tag", required_argument, NULL, 'u'},
            {"dry-run", no_argument, NULL, 'n'},
            {"verbose", no_argument, NULL, 'v'},
            {"help", no_argument, NULL, 'h'}
    };

    log_msg("Parsing commandline flags", DEBUG);


    while ((opt = getopt_long(argc, argv, ":dq:f:m:o:t:u:nvh", cl_opts, NULL)) != -1) {
        switch (opt) {
            case 'd':
                dedup = true;
                break;
            case 'q':
                if (optarg != NULL) {
                    mapq = strtol(optarg, &mapqstr, 10);
                }
                break;
            case 'f':
                bampath = optarg;
                break;
            case 'm':
                metapath = optarg;
                break;
            case 'o':
                oprefix = optarg;
                break;
            case 't':
                for (int i = 0; i < 2; i++) {
                    filter_tag[i] = optarg[i];
                }
                break;
            case 'u':
                for (int i = 0; i < 2; i++) {
                    umi_tag[i] = optarg[i];
                }
                break;
            case 'v':
                verbose = true;
                break;
            case 'n':
                dryrun = true;
                break;
            case 'h':
                show_usage("regular");
                return 0;
            case ':':
                show_usage("missing");
                return 1;
            case '?':
                show_usage("unknown");
                return 1;
            default:
                show_usage("regular");
                return 1;
        }
    }

    // Bypass input for testing during development
    if (dev) {
        log_msg("Dev mode on: Bypassing command line input", DEBUG);
        bampath = "../data/test.bam";
        metapath = "../data/test.txt";
        oprefix = "./huge_out/";
    }

    // Error-out if missing arguments
    if (bampath == NULL || metapath == NULL) {
        log_msg("Error: Missing required arguments", ERROR);
        show_usage("regular");
        return 1;
    }

    // If there's no output prefix set, export files in the work directory
    if (oprefix == NULL) {
        oprefix = "./";
    }

    // If the output prefix does not end with /, add it.
    uint16_t oplen = strlen(oprefix) - 1;
    if (oprefix[oplen] != '/')
        oprefix = strcat(oprefix, "/");


    if (verbose || dryrun) {
        fprintf(stderr, "- Run condition:\n\n");
        fprintf(stderr, "\tInput bam: %s\n", bampath);
        fprintf(stderr, "\tInput metadata: %s\n", metapath);
        fprintf(stderr, "\tMAPQ threshold: %d\n", mapq);
        fprintf(stderr, "\tRead tag to filter: %s\n", filter_tag);
        fprintf(stderr, "\tUMI tag to filter: %s\n", umi_tag);
        fprintf(stderr, "\tOutput prefix: %s\n", oprefix);
        if (dedup) {
            fprintf(stderr, "\tRunning **with** deduplication.\n\n");
        } else {
            fprintf(stderr, "\tRunning **without** deduplication.\n\n");
        }
    }

    if (dryrun) {
        fprintf(stderr, "\t==========================================================\n");
        fprintf(stderr, "\t= This is a dry-run (-n/--dry-run). Nothing is executed. =\n");
        fprintf(stderr, "\t==========================================================\n");
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
    };


    ////////// bam related //////////////////////////////////////////////
    // Open input bam file from CellRanger
    // Remember to close file handle!
    log_msg("Reading input BAM file: %s", INFO, bampath);
    samFile *fp = sam_open(bampath, "r");


    // Extract header
    log_msg("Reading SAM header", DEBUG);
    sam_hdr_t *header = sam_hdr_read(fp);


    // Iterating variables for reads in the bam file
    log_msg("Create temporary read", DEBUG);
    bam1_t *read = bam_init1();

    // Prepare a read-tag-to-label hash table from a metadata table
    log_msg("Loading barcode-cluster information from metadata: %s", INFO, metapath);
    struct rt2label *r2l = hash_readtag((char *) metapath);

    if (r2l == NULL) {
        log_msg("Failed to hash the metadata.", ERROR);
        return 1;
    }

    // Prepare a label-to-file-handle hash table from the above
    log_msg("Preparing output BAM files", INFO);
    struct label2fp *l2fp = NULL;
    l2fp = hash_labels(r2l, oprefix, header);

    log_msg("Starting to split files", INFO);
    // Iterate through the rt's and write to corresponding file handles.
    // Iterate through reads from input bam
    struct rt2label *lout;
    struct label2fp *fout;
    struct dedup *prev_reads = NULL, *find_reads;
    int ret;


    int64_t chunk_size = 1000000;
    log_msg("Processing %lld reads per chunk", INFO, chunk_size);

    // Allocate heap memory for reads to sort
    log_msg("Preparing read chunks for sorting", DEBUG);

    sam_read_t *chunk = chunk_init(chunk_size);
    if (NULL == chunk) {
        return -1;
    }

    char* tmpdir = process_bam(
                fp, header, chunk, chunk_size, oprefix
            );
    if (strcmp(tmpdir, "1") == 0) {
        return -1;
    }
    purge_tempdir(tmpdir);

    // Done processing
    log_msg("Completed sorting all chunks", INFO);
    chunk_destroy(chunk, chunk_size);
    free(chunk);


    /*
    while ((ret = sam_read1(fp, header, read)) >= 0) {
        // Only deal with reads with both corrected CBCs (and UMIs if dedup is on)
        // Extract corrected CBC from the read
        if (bam_aux_get(read, filter_tag) == NULL) {
            // If a read has no CBC, just continue.
            continue;
        }

        if (dedup && bam_aux_get(read, umi_tag) == NULL) {
            // If UMI is needed for dedup, ignore the reads without it.
            continue;
        }

        unsigned char * cbc = bam_aux_get(read, filter_tag) + 1;

        // Query read-tag-to-label table
        HASH_FIND_STR(r2l, (char *) cbc, lout);

        // If the read tag is legit
        if (lout) {
            // Query the read-tag-to-output table
            HASH_FIND_STR(l2fp, lout->label, fout);
            if (fout) {
                if (read->core.qual >= mapq) {
                    if (dedup) {
                        // Construct an identifier = CB + UB
                        unsigned char * umi = bam_aux_get(read, umi_tag) + 1;
                        char id[ID_LENGTH];
                        strcpy(id, (char *) cbc);
                        strcat(id, (char *) umi);
                        HASH_FIND_STR(prev_reads, (char *) id, find_reads);
                        if (find_reads) {
                            // If the CBC + UMI combination was observed before
                            // Just continue.
                            continue;
                        }
                        prev_reads = hash_cbumi(prev_reads, id);
                    }
                    int8_t write_to_bam;
                    write_to_bam = sam_write1(fout->fp, header, read);
                    if (write_to_bam < 0) {
                        fprintf(stderr, "[LOUT] Fail to write bam file.");
                    }
                }
            }
        }

    }
*/
    // Close file handles
    struct label2fp *qs, *qtmp;
    HASH_ITER(hh, l2fp, qs, qtmp) {
        sam_close(qs->fp);
//        printf("label is %s\n", qs->label);
        HASH_DEL(l2fp, qs);
        free(qs);
    }

    // free the hash table contents
    struct rt2label *s, *tmp;
    HASH_ITER(hh, r2l, s, tmp) {
        HASH_DEL(r2l, s);
        free(s);
    }

    // Close CB-UMI hash table
    if (dedup) {
        struct dedup *ds, *dtmp;
        HASH_ITER(hh, prev_reads, ds, dtmp) {
            HASH_DEL(prev_reads, ds);
            free(ds);
        }
    }


    // Close opened bam files
    bam_destroy1(read);
    bam_hdr_destroy(header);
    sam_close(fp);



    return 0;
}
