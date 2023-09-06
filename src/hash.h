//
// Created by Yen-Chung Chen on 2/23/23.
//
#ifndef SCBAMSPLIT_HASH_H
#define SCBAMSPLIT_HASH_H
#include "htslib/sam.h"
#include "uthash.h"
#include "shared_const.h"

typedef struct {
    char rt[MAX_LINE_LENGTH];             /* key (string is WITHIN the structure) */
    char label[MAX_LINE_LENGTH];
    UT_hash_handle hh;         /* makes this structure hashable */
} rt2label;

typedef struct {
    char label[MAX_LINE_LENGTH];             /* key (string is WITHIN the structure) */
    samFile* fp;
    UT_hash_handle hh;         /* makes this structure hashable */
} label2fp ;

rt2label* hash_readtag(char *path);
label2fp* hash_labels(rt2label *r2l, const char *prefix, sam_hdr_t *header);


#endif //SCBAMSPLIT_HASH_H
