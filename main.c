#include <string.h>  /* strcpy */
#include <stdlib.h>  /* malloc */
#include <stdio.h>   /* printf */
#include "uthash.h"

const int MAX_LINE_LENGTH = 256;

struct rt2label {
    char rt[MAX_LINE_LENGTH];             /* key (string is WITHIN the structure) */
    char label[MAX_LINE_LENGTH];
    UT_hash_handle hh;         /* makes this structure hashable */
};

struct label2fp {
    char label[MAX_LINE_LENGTH];             /* key (string is WITHIN the structure) */
    FILE* fp;
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

struct label2fp* hash_labels(struct rt2label *r2l, const char *prefix, int MAX_LINE_LENGTH) {
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

        // Create file handle from the path generated above
        // These handles must be closed manually!
        new_l2f->fp = fopen(outpath, "w");

        // Add the new element to the hash table
        HASH_ADD_STR(l2f, label, new_l2f);
    }
    return l2f;
}


int main(int argc, char *argv[]) {
    // Test only
    const char path[] = "/Users/ycc/Dropbox/My Mac (10-17-233-196.dynapool.wireless.nyu.edu)/Desktop/scbamsplit/data/test.txt";

    // Main code
    const char prefix[] = "test_";

    // Prepare a read-tag-to-label hash table from a metadata table
    struct rt2label *r2l = rt_hash((char *) path, MAX_LINE_LENGTH);

    if (r2l == NULL) {
        printf("Failed to hash the metadata.\n");
        return 1;
    }

    // Prepare a label-to-file-handle hash table from the above
    struct label2fp *l2fp = NULL;
    l2fp = hash_labels(r2l, prefix, MAX_LINE_LENGTH);

    // Iterate through the rt's and write to corresponding file handles.
    struct rt2label *s, *tmp;
    struct label2fp *fout;
    HASH_ITER(hh, r2l, s, tmp) {
        HASH_FIND_STR(l2fp, s->label, fout);
        fprintf(fout->fp, "%s, %s\n", s->rt, s->label);
    }


    // Close file handles
    struct label2fp *qs, *qtmp;
    HASH_ITER(hh, l2fp, qs, qtmp) {
//        fprintf(qs->fp, "label = %s", qs->label);
        fclose(qs->fp);
//        printf("label is %s\n", qs->label);
        HASH_DEL(l2fp, qs);
        free(qs);
    }


    // free the hash table contents
    HASH_ITER(hh, r2l, s, tmp) {
        HASH_DEL(r2l, s);
        free(s);
    }


    return 0;
}