#include <string.h>  /* strcpy */
#include <stdlib.h>  /* malloc */
#include <stdio.h>   /* printf */
#include "htslib/sam.h" /* Use htslib to interact with bam files */
#include "uthash.h"  /* hash table */
#include "hash.h"    /* Defining the hash tables actually used */
#include "shared_const.h" /* Defining shared constants */

int main(int argc, char *argv[]) {
    // Commandline argument processing
    printf("Running %s...\n", argv[0]);

    if( argc == 4) {
        printf(
                "Processing bam: %s\nProcessing metadata: %s\nOutput prefix: %s\n", argv[1], argv[2], argv[3]
                );
    }
    else if( argc > 3 ) {
        printf("Too many arguments supplied.\n");
        return 1;
    }
    else {
        printf("Three arguments expected.\n");
        return 1;
    }
    // Main code


    ////////// bam related//////////////////////////////////////////////

    // Open input bam file from CellRanger
    // Remember to close file handle!
//    samFile *fp = sam_open(bampath, "r");
    samFile *fp = sam_open(argv[1], "r");

    // Extract header
    sam_hdr_t *header = sam_hdr_read(fp);

    // Iterating variables for reads in the bam file
    bam1_t *read = bam_init1();

    // Prepare a read-tag-to-label hash table from a metadata table
    struct rt2label *r2l = rt_hash((char *) argv[2], MAX_LINE_LENGTH);

    if (r2l == NULL) {
        printf("Failed to hash the metadata.\n");
        return 1;
    }

    // Prepare a label-to-file-handle hash table from the above
    struct label2fp *l2fp = NULL;
    l2fp = hash_labels(r2l, argv[3], MAX_LINE_LENGTH, header);

    // Iterate through the rt's and write to corresponding file handles.
    // Iterate through reads from input bam

    struct rt2label *lout;
    struct label2fp *fout;
    int ret;
    while ((ret = sam_read1(fp, header, read)) >= 0) {
        if (bam_aux_get(read, "CB") != NULL) {
            // Extract corrected CBC from the read
            unsigned char * cbc = bam_aux_get(read, "CB") + 1;

            // Query read-tag-to-label table
            HASH_FIND_STR(r2l, (char *) cbc, lout);

            // If the read tag is legit
            if (lout) {
                // Query the read-tag-to-output table
                HASH_FIND_STR(l2fp, lout->label, fout);
                if (fout) {
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

    struct rt2label *s, *tmp;
    // free the hash table contents
    HASH_ITER(hh, r2l, s, tmp) {
        HASH_DEL(r2l, s);
        free(s);
    }

    // Close opened bam files
    bam_destroy1(read);
    bam_hdr_destroy(header);
    sam_close(fp);


    return 0;
}
