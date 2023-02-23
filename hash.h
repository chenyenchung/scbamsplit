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

struct dedup {
    char id[MAX_LINE_LENGTH];
    UT_hash_handle hh;
};

struct rt2label* hash_readtag(char *path);
struct label2fp* hash_labels(struct rt2label *r2l, const char *prefix, sam_hdr_t *header);
struct dedup* hash_cbumi(struct dedup *dedup_t, char *id);


#endif //SCBAMSPLIT_HASH_H
