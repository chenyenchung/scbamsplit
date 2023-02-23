//
// Created by Yen-Chung Chen on 2/23/23.
//
#include "htslib/sam.h"
#include "shared_const.h"
#ifndef SCBAMSPLIT_HASH_H
#define SCBAMSPLIT_HASH_H

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

struct rt2label* rt_hash(char *path, const int max_length);
struct label2fp* hash_labels(struct rt2label *r2l, const char *prefix, int max_length, sam_hdr_t *header);


#endif //SCBAMSPLIT_HASH_H
