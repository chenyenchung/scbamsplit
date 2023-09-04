//
// Created by Yen-Chung Chen on 8/23/23.
//

#ifndef SCBAMSPLIT_SORT_H
#define SCBAMSPLIT_SORT_H
#include "htslib/sam.h"

extern int64_t KEY_SIZE;
extern int64_t RN_SIZE;
extern int64_t CB_LENGTH;
extern int64_t UB_LENGTH;

typedef struct {
    char *key;
    bam1_t *read;
} sam_read_t;

enum location {
    READ_TAG,
    READ_NAME
};

typedef struct {
    enum location location;
    char *tag_name;
    char *sep;
    uint8_t field;
    uint8_t length;
} tag_meta_t;

int8_t fetch_tag(bam1_t *read, char *tag, char *tag_ptr);
int8_t get_CB (bam1_t *read, tag_meta_t* info, char* tag_ptr);
int8_t get_UB (bam1_t *read, tag_meta_t* info, char* tag_ptr);
void set_CB(tag_meta_t *tag_meta, char *platform);
void set_UB(tag_meta_t *tag_meta, char *platform);
tag_meta_t *initialize_tag_meta();
void print_tag_meta(tag_meta_t *tag_meta, const char *header);
void destroy_tag_meta(tag_meta_t *tag_meta);

int64_t
fill_chunk(samFile *fp, sam_hdr_t *header, sam_read_t **read_array, int64_t chunk_size, int16_t qthres,
           tag_meta_t *cb_meta, tag_meta_t *ub_meta);
void sort_chunk(sam_read_t **reads, int64_t chunk_size);
sam_read_t** chunk_init(uint32_t chunk_size);
void chunk_destroy(sam_read_t **read_array, uint32_t chunk_size);
char *process_bam(samFile *fp, sam_hdr_t *header, sam_read_t **chunk, int64_t chunk_size, char *oprefix, int16_t qthres,
                  tag_meta_t *cb_meta, tag_meta_t *ub_meta);


#endif //SCBAMSPLIT_SORT_H
