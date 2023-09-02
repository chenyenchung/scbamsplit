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

int show_usage(const char* type);
int create_directory(char* pathname);
char * create_tempdir(char *dir);
int8_t purge_tempdir(char *tmpdir);

extern log_level_t OUT_LEVEL;
extern char *OUT_PATH;
extern bool dev;
#endif //SCBAMSPLIT_UTILS_H
