#include <string.h>  /* strcpy */
#include <stdlib.h>  /* malloc */
#include <stdio.h>   /* printf */
#include <getopt.h>  /* getopt */
#include <stdbool.h> /* Define boolean type */
#include "htslib/sam.h" /* Use htslib to interact with bam files, imports stdint.h as well */
#include "uthash.h"  /* hash table */
#include "hash.h"    /* Defining the hash tables actually used */
#include "utils.h" /* Show help and create output dir */
#include "shared_const.h" /* Defining shared constants */
#include "sort.h"

int main(int argc, char *argv[]) {

    // Use a flag to bypass commandline input during development
    bool dev = true;

    int opt;
    uint16_t mapq = 0;
    char *mapqstr = NULL;
    bool dedup = false, dryrun = false, verbose = false;
    char *bampath = NULL;
    char *metapath = NULL;
    char *oprefix = NULL;
    char filter_tag[3] = "CB";
    char umi_tag[3] = "UB";
    char* logpath = "";

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

    log_message("Parsing commandline flags", DEBUG, logpath, OUT_LEVEL);


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
        log_message(
                "Dev mode on: Bypassing command line input",
                DEBUG,
                logpath,
                OUT_LEVEL
                );
        bampath = "../data/test.bam";
        metapath = "../data/test.txt";
        oprefix = "./huge_out/";
    }

    // Error-out if missing arguments
    if (bampath == NULL || metapath == NULL) {
        log_message(
                "Error: Missing required arguments",
                ERROR,
                logpath,
                OUT_LEVEL
        );
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
    } else {

    }

    if (dryrun) {
        fprintf(stderr, "\t==========================================================\n");
        fprintf(stderr, "\t= This is a dry-run (-n/--dry-run). Nothing is executed. =\n");
        fprintf(stderr, "\t==========================================================\n");
        return 0;
    }

    // Create output folder if it does not exist
    // If it exists, ask the user for confirmation to prevent unexpected overwriting
    log_message("Creating output directory", INFO, logpath, OUT_LEVEL);
    if (create_folder(oprefix) != 0) return 1;


    ////////// bam related //////////////////////////////////////////////
    // Open input bam file from CellRanger
    // Remember to close file handle!
    log_message("Reading input BAM file: %s", INFO, logpath, OUT_LEVEL, bampath);
    samFile *fp = sam_open(bampath, "r");


    // Extract header
    log_message("Reading SAM header", DEBUG, logpath, OUT_LEVEL);
    sam_hdr_t *header = sam_hdr_read(fp);

    // Iterating variables for reads in the bam file
    log_message("Create temporary read", DEBUG, logpath, OUT_LEVEL);
    bam1_t *read = bam_init1();

    // Prepare a read-tag-to-label hash table from a metadata table
    log_message(
            "Loading barcode-cluster information from metadata: %s",
            INFO, logpath, OUT_LEVEL,
            metapath
            );
    struct rt2label *r2l = hash_readtag((char *) metapath);

    if (r2l == NULL) {
        log_message("Failed to hash the metadata.", ERROR, logpath, OUT_LEVEL);
        return 1;
    }

    // Prepare a label-to-file-handle hash table from the above
    log_message("Preparing output BAM files", INFO, logpath, OUT_LEVEL);
    struct label2fp *l2fp = NULL;
    l2fp = hash_labels(r2l, oprefix, header);

    log_message("Starting to split files", INFO, logpath, OUT_LEVEL);
    // Iterate through the rt's and write to corresponding file handles.
    // Iterate through reads from input bam
    struct rt2label *lout;
    struct label2fp *fout;
    struct dedup *prev_reads = NULL, *find_reads;
    int ret;


    int64_t chunk_size = 1000000;
    log_message("Processing %lld reads per batch", INFO, logpath, OUT_LEVEL, chunk_size);

    // Allocate heap memory for reads to sort
    log_message("Preparing read chunks for sorting", DEBUG, logpath, OUT_LEVEL);
    sam_read *chunk;
    chunk = malloc(sizeof(sam_read) * chunk_size);
    if (chunk == NULL) {
        log_message("Insufficient memory. Please decrease chunk size used (%d)",
                    ERROR, logpath, OUT_LEVEL, chunk_size);
        return -1;
    }

    chunk_init(chunk, chunk_size);

    int64_t size_retrieved = chunk_size;

    int chunk_num = 0;
    // Create temporary file dir for sorted chunks
    char *tmpdir = create_tempdir(oprefix);

    while (size_retrieved == chunk_size) {
        chunk_num += 1;
        log_message("Chunk #%lld filling", DEBUG, logpath, OUT_LEVEL, chunk_num);
        size_retrieved = fill_chunk(fp, header, chunk, chunk_size);
        if (size_retrieved < -1) {
            log_message("Insufficient RN size (%d). Please increase RN size to at least %d",
                        ERROR, logpath, OUT_LEVEL, RN_SIZE, -size_retrieved + 1);
            return 1;
        }
        sort_chunk(chunk, size_retrieved);
        sam_hdr_change_HD(header, "SO", "unknown");

        char *tname;
        tname = malloc(sizeof(char) * (strlen(tmpdir) + 16)); // tmp[/chunkXXXXX.bam]
        char *uid[10];
        strcpy(tname, tmpdir);
        strcat(tname, "/chunk");
        sprintf(uid, "%05d.bam", chunk_num);
        strcat(tname, uid);

        htsFile* tfp = sam_open(tname, "wb");
        free(tname);

        int write_status = sam_hdr_write(tfp, header);
        if (write_status != 0) {
            log_message(
                    "Fail to write SAM header into temporary files",
                    ERROR, logpath, OUT_LEVEL
                    );
            return 1;
        }

        int write_to_bam = 0;
        for (int64_t i = 0; i < size_retrieved; i++) {
            write_to_bam = sam_write1(tfp, header, chunk[i].read);
            if (write_to_bam == -1) {
                log_message(
                        "Fail to write sorted reads into temporary files",
                        ERROR, logpath, OUT_LEVEL
                );
                return 1;
            }
        }
        sam_close(tfp);
    }
    free(tmpdir);
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
