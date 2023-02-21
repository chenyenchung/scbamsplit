#include <string.h>  /* strcpy */
#include <stdlib.h>  /* malloc */
#include <stdio.h>   /* printf */
#include "uthash.h"
#include "htslib/sam.h"

const int MAX_LINE_LENGTH = 256;

struct rt2label {
    char rt[MAX_LINE_LENGTH];             /* key (string is WITHIN the structure) */
    char label[MAX_LINE_LENGTH];
    UT_hash_handle hh;         /* makes this structure hashable */
};

struct label2fp {
    char label[MAX_LINE_LENGTH];             /* key (string is WITHIN the structure) */
    samFile* fp;
    UT_hash_handle hh;         /* makes this structure hashable */
};

struct rt2label* rt_hash(char *path, const int MAX_LINE_LENGTH) {
    // Open meta data file
    FILE* meta_fp = fopen(path, "r");

    if (meta_fp == NULL) {
        // Exit and print error message if the file does not exist
        printf("Error: Cannot open file %s.\n", path);
        return NULL;
    }

    // Declare a read-tag-to-label hash table
    struct rt2label *r2l = NULL;

    // Read every line
    char meta_line[MAX_LINE_LENGTH]; // Temporary variable to store each line
    int first_line = 1;
    while (fgets(meta_line, MAX_LINE_LENGTH, meta_fp) != NULL) {
        // Assuming header and skip it
        if (first_line) {
            first_line--;
            continue;
        }
        // Strip the linebreak
        meta_line[strcspn(meta_line, "\n")] = 0;

        // Tokenize by comma
        char* tokens; // Temporary variable for iterating tokens
        int field_num = 0; // Counting numbers to examine if expected field numbers are present
        char trt[MAX_LINE_LENGTH]; // Temporary variable for read tag content
        char tlabel[MAX_LINE_LENGTH]; // Temporary variable for corresponding label content
        tokens = strtok(meta_line, ",");

        // Iterate through tokens
        while (tokens != NULL) {
            switch(field_num) {
                case 0:
                    // Expect the first field to be read tags
                    strcpy(trt, tokens);
                    break;
                case 1:
                    // Expect the second to be labels
                    strcpy(tlabel, tokens);
                    break;
                default:
                    printf("Error: There are %d fields in the metadata. Expecting 2.\n", field_num);
                    return NULL;

            }
            field_num ++;
            tokens = strtok(NULL, ",");
        }

        // Deal with metadata that has < 2 fields
        if (field_num < 2) {
            printf("Error: There is only %d field in the metadata.", field_num);
            return NULL;
        }

        // Prepare a new element for the hash table
        struct rt2label *s;
        // These hash tables must be freed through iteration!
        s = (struct rt2label *)malloc(sizeof *s);

        // Assign the read tag and label as a key-value pair
        strcpy(s->rt, trt);
        strcpy(s->label, tlabel);

        HASH_ADD_STR(r2l, rt, s);
    }
    fclose(meta_fp);
    return r2l;
}

struct label2fp* hash_labels(struct rt2label *r2l, const char *prefix, int MAX_LINE_LENGTH, sam_hdr_t *header) {
    struct rt2label *s; // Declare a temporary variable to iterate over the rt2label table
    struct label2fp *l2f = NULL; // Initialize l2f to hold the label2fp table
    struct label2fp *tmp, *new_l2f; // Declare temporary variables for HASH_FIND_STR

    // Loop over each element in the rt2label table
    for (s = r2l; s != NULL; s = s->hh.next) {
        HASH_FIND_STR(l2f, s->label, tmp);
        // If a read tag-associated label already has a file handle
        // Just move on
        if (tmp) {
            continue;
        }

        // If the label is new, malloc() for a new element
        // These hash tables must be freed through iteration!
        new_l2f = (struct label2fp *)malloc(sizeof *new_l2f);

        // Populate the label (key) for the new element
        strcpy(new_l2f->label, s->label);

        // Prepare output file path
        char outpath[MAX_LINE_LENGTH];
        strcpy(outpath, prefix);
        strcat(outpath, s->label);
        strcat(outpath, ".sam");

        // Create file handle from the path generated above
        // These handles must be closed manually!
        new_l2f->fp = sam_open(outpath, "w");

        // Populate header
        sam_hdr_write(new_l2f->fp, header);

        // Add the new element to the hash table
        HASH_ADD_STR(l2f, label, new_l2f);
    }
    return l2f;
}


int main(int argc, char *argv[]) {
    // Test only
//    const char path[] = "/Users/ycc/Dropbox/My Mac (10-17-233-196.dynapool.wireless.nyu.edu)/Desktop/scbamsplit/data/test.txt";
//    const char bampath[] = "/Users/ycc/Dropbox/My Mac (10-17-233-196.dynapool.wireless.nyu.edu)/Desktop/scbamsplit/data/test.bam";
//    const char prefix[] = "test_";

    // Commandline argument processing
    printf("Running %s...\n", argv[0]);

    if( argc == 4) {
        printf(
                "Processing bam: %s\nProcessing metadata: %s\nOutput prefix: %s\n", argv[1], argv[2], argv[3]
                );
    }
    else if( argc > 3 ) {
        printf("Too many arguments supplied.\n");
    }
    else {
        printf("Three arguments expected.\n");
    }

    // Main code


    ////////// bam related//////////////////////////////////////////////

    // Open input bam file from CellRanger
    // Remember to close file handle!
//    samFile *fp = sam_open(bampath, "r");
    samFile *fp = sam_open(argv[1], "r");

    // Extract header
    sam_hdr_t *header = sam_hdr_read(fp);

    // Iterating variables
    bam1_t *read = bam_init1();
    int ret;

    // Prepare output file
    // Remember to close file handle!



    ////////////////////////////////////////////////////////////////////////////////////////

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
    while ((ret = sam_read1(fp, header, read)) >= 0) {
        if (bam_aux_get(read, "CB") != NULL) {
            // Extract corrected CBC from the read
            unsigned char * cbc = bam_aux_get(read, "CB") + 1;

            // Query read-tag-to-label table
            HASH_FIND_STR(r2l, cbc, lout);

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