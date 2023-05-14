#include <string.h>  /* strcpy */
#include <stdlib.h>  /* malloc */
#include <stdio.h>   /* printf */
#include <getopt.h>  /* getopt */
#include <stdbool.h>
#include "htslib/sam.h" /* Use htslib to interact with bam files */
#include "uthash.h"  /* hash table */
#include "hash.h"    /* Defining the hash tables actually used */
#include "utils.h"
#include "shared_const.h" /* Defining shared constants */

int main(int argc, char *argv[]) {

    int opt;
    long int mapq = 0;
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

    // Error-out if missing arguments
    if (bampath == NULL || metapath == NULL) {
       fprintf(stderr, "Error: Missing required arguments.\n\n");
       show_usage("regular");
       return 1;
    }

    // If there's no output prefix set, export files in the work directory
    if (oprefix == NULL) {
        oprefix = "./";
    }

    if (verbose || dryrun) {
        fprintf(stderr, "- Run condition:\n\n");
        fprintf(stderr, "\tInput bam: %s\n", bampath);
        fprintf(stderr, "\tInput metadata: %s\n", metapath);
        fprintf(stderr, "\tMAPQ threshold: %ld\n", mapq);
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

    ////////// bam related//////////////////////////////////////////////

    // Open input bam file from CellRanger
    // Remember to close file handle!
    samFile *fp = sam_open(bampath, "r");

    // Extract header
    sam_hdr_t *header = sam_hdr_read(fp);

    // Iterating variables for reads in the bam file
    bam1_t *read = bam_init1();

    // Prepare a read-tag-to-label hash table from a metadata table
    struct rt2label *r2l = hash_readtag((char *) metapath);

    if (r2l == NULL) {
        printf("Failed to hash the metadata.\n");
        return 1;
    }

    // Prepare a label-to-file-handle hash table from the above
    struct label2fp *l2fp = NULL;
    l2fp = hash_labels(r2l, oprefix, header);

    // Iterate through the rt's and write to corresponding file handles.
    // Iterate through reads from input bam
    struct rt2label *lout;
    struct label2fp *fout;
    struct dedup *prev_reads = NULL, *find_reads;
    int ret;

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
                    sam_write1(fout->fp, header, read);
                }
            }
        }

    }

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
