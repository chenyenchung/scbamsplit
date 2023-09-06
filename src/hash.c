//
// Created by Yen-Chung Chen on 2/23/23.
//
#include "hash.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include "utils.h"
rt2label* hash_readtag(char *path) {
    // Open meta data file
    FILE* meta_fp = fopen(path, "r");

    if (meta_fp == NULL) {
        // Exit and print error message if the file does not exist
        log_msg("Cannot open file (%s)", ERROR, path);
        return NULL;
    }

    // Declare a read-tag-to-label hash table
    rt2label *r2l = NULL;

    // Read every line
    char meta_line[MAX_LINE_LENGTH]; // Temporary variable to store each line
    bool first_line = true;
    char* tokens; // Temporary variable for iterating tokens
    uint32_t field_num = 0; // Counting numbers to examine if expected field numbers are present
    char trt[MAX_LINE_LENGTH]; // Temporary variable for read tag content
    char tlabel[MAX_LINE_LENGTH]; // Temporary variable for corresponding label content

    while (fgets(meta_line, MAX_LINE_LENGTH, meta_fp) != NULL) {
        // Assuming header and skip it
        if (first_line) {
            first_line = false;
            continue;
        }


        // Strip the linebreak
        meta_line[strcspn(meta_line, "\n")] = 0;

        // Tokenize by comma
        tokens = strtok(meta_line, ",");

        // Reset field number
        field_num = 0;

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
                    log_msg("There are %d fields in the metadata but only 2 are expected",
                            ERROR, field_num);
                    return NULL;

            }
            field_num ++;
            tokens = strtok(NULL, ",");
        }

        // Deal with metadata that has < 2 fields
        if (field_num < 2) {
            log_msg("There is only %d field in the metadata (expecting 2)", ERROR, field_num);
            return NULL;
        }

        // Prepare a new element for the hash table
        rt2label *s;
        // These hash tables must be freed through iteration!
        s = (rt2label *)calloc(1, sizeof(rt2label));

        // Assign the read tag and label as a key-value pair
        strcpy(s->rt, trt);
        strcpy(s->label, tlabel);

        HASH_ADD_STR(r2l, rt, s);
    }
    fclose(meta_fp);
    return r2l;
}

label2fp* hash_labels(rt2label *r2l, const char *prefix, sam_hdr_t *header) {
    rt2label *s; // Declare a temporary variable to iterate over the rt2label table
    label2fp *l2f = NULL; // Initialize l2f to hold the label2fp table
    label2fp *tmp, *new_l2f; // Declare temporary variables for HASH_FIND_STR

    // Loop over each element in the rt2label table
    for (s = r2l; s != NULL; s = s->hh.next) {
        HASH_FIND_STR(l2f, s->label, tmp);
        // If a read tag-associated label already has a file handle
        // Just move on
        if (tmp) {
            continue;
        }

        // If the label is new, calloc() for a new element
        // These hash tables must be freed through iteration!
        new_l2f = (label2fp *)calloc(1, sizeof(label2fp));

        // Populate the label (key) for the new element
        strcpy(new_l2f->label, s->label);

        // Prepare output file path
        char outpath[MAX_LINE_LENGTH];
        uint32_t label_size = strlen(s->label);
        char label_corrected[label_size + 1];

        strcpy(label_corrected, s->label);
        for (int i = 0; i < label_size; i++) {
            // Replace slashes with hyphens when creating output files
            if (label_corrected[i] == '/') {
                label_corrected[i] = '-';
            }
        }

        // Concatenate output path
        strcpy(outpath, prefix);
        strcat(outpath, label_corrected);
        strcat(outpath, ".bam");

        // Create file handle from the path generated above
        // These handles must be closed manually!
        new_l2f->fp = sam_open(outpath, "wb");

        // Populate header
        uint8_t hdr_write;
        hdr_write = sam_hdr_write(new_l2f->fp, header);

        if (hdr_write != 0) {
            log_msg("Fail to prepare individual output files", ERROR);
            return NULL;
        }

        // Add the new element to the hash table
        HASH_ADD_STR(l2f, label, new_l2f);
    }
    return l2f;
}
