//
// Created by Yen-Chung Chen on 2/24/23.
//
#include <stdbool.h>
#include <stdint.h>
#include "hash.h"
#include "htslib/sam.h"
#include "sort.h"

#ifndef SCBAMSPLIT_UTILS_H
#define SCBAMSPLIT_UTILS_H
///////////// Logging utilities ////////////////////
typedef enum {
    DEBUG = 5,
    ERROR = 0,
    WARNING = 1,
    INFO = 3
} log_level_t;

void log_message(char* message, log_level_t level, char* log_path, log_level_t OUT_LEVEL, ...);
#define log_msg(...) log_message(OUT_PATH, OUT_LEVEL, __VA_ARGS__)
///////////////////////////////////////////////////

///////////// String utilities ////////////////////
// Defining a string array with fixed element length and number
typedef struct {
    char **str_arr;
    uint32_t length;
    uint16_t str_length;
} str_vec_t;

str_vec_t* str_vec_init (int64_t n, uint16_t str_len);
int8_t str_vec_destroy(str_vec_t *ptr);
///////////////////////////////////////////////////

void show_usage();
int create_directory(char* pathname);
char * create_tempdir(char *dir);
str_vec_t * get_bams(char *tmpdir);
char * tname_init(char * tmpdir, char * prefix, int32_t uid_length, uint32_t oid);
char * merge_bams(char * tmpdir);

int8_t read_dump(rt2label *r2l, rt2label *lout,
                 label2fp *l2fp, label2fp *fout,
                 char * this_CB, sam_hdr_t *header, bam1_t *read);
int8_t deduped_dump(rt2label *r2l, rt2label *lout, label2fp *l2fp, label2fp *fout, char *tmpdir, char *sorted_path,
                    bam1_t *read, char *bc_tag, char *umi_tag, tag_meta_t *cb_meta, tag_meta_t *ub_meta);

extern log_level_t OUT_LEVEL;
extern char *OUT_PATH;
extern bool dev;
extern char* LEVEL_FLAG[6];
#endif //SCBAMSPLIT_UTILS_H
