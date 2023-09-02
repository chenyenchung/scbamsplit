//
// Created by Yen-Chung Chen on 2/24/23.
//
#include <stdbool.h>
#include <stdint.h>

#ifndef SCBAMSPLIT_UTILS_H
#define SCBAMSPLIT_UTILS_H
// Logging utilities
typedef enum {
    DEBUG = 5,
    ERROR = 0,
    WARNING = 1,
    INFO = 3
} log_level_t;
void log_message(char* message, log_level_t level, char* log_path, log_level_t OUT_LEVEL, ...);
#define log_msg(...) log_message(OUT_PATH, OUT_LEVEL, __VA_ARGS__)

typedef struct {
    char **str_arr;
    uint32_t length;
    uint16_t str_length;
} str_vec_t;

str_vec_t* str_vec_init (int64_t n, uint16_t str_len);
int8_t str_vec_free(str_vec_t *ptr);

int show_usage(const char* type);
int create_directory(char* pathname);
char * create_tempdir(char *dir);
int8_t purge_tempdir(char *tmpdir);
str_vec_t * get_bams(char *tmpdir);
char * tname_init(char * tmpdir, char * prefix, int32_t uid_length, uint32_t oid);
int8_t merge_bams(char * tmpdir);

#define EARLY_EXIT_MERGE free(key1); \
free(key2); \
bam_destroy1(r1); \
bam_destroy1(r2); \
sam_hdr_destroy(header1); \
sam_hdr_destroy(header2); \
sam_close(fp1); \
sam_close(fp2); \
sam_close(tfp); \
free(ff1); \
free(ff2);
extern log_level_t OUT_LEVEL;
extern char *OUT_PATH;
extern bool dev;
#endif //SCBAMSPLIT_UTILS_H
