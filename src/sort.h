//
// Created by Yen-Chung Chen on 8/23/23.
//

#ifndef SCBAMSPLIT_SORT_H
#define SCBAMSPLIT_SORT_H
#include "htslib/sam.h"

extern uint16_t KEY_SIZE;
extern uint16_t RN_SIZE;
extern uint8_t CB_LENGTH;
extern uint8_t UB_LENGTH;

typedef struct {
    char *key;
    bam1_t *read;

} sam_read;

int64_t fill_chunk(
        samFile *fp, sam_hdr_t *header, sam_read read_array[],
        int64_t chunk_size
        );
int sort_chunk(sam_read read_array[], int64_t chunk_size);
void chunk_init(sam_read *read_array, uint32_t chunk_size);
void chunk_destroy(sam_read *read_array, uint32_t chunk_size);

#endif //SCBAMSPLIT_SORT_H
