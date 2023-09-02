//
// Created by Yen-Chung Chen on 2/24/23.
//
#include <stdbool.h>

#ifndef SCBAMSPLIT_UTILS_H
#define SCBAMSPLIT_UTILS_H
int show_usage(const char* type);
int create_directory(char* pathname);
char * create_tempdir(char *dir);

// Logging utilities
typedef enum {
    DEBUG = 5,
    ERROR = 0,
    WARNING = 1,
    INFO = 3
} log_level_t;

extern log_level_t OUT_LEVEL;
extern char *OUT_PATH;
extern bool dev;


void log_message(char* message, log_level_t level, char* log_path, log_level_t out_level, ...);
#endif //SCBAMSPLIT_UTILS_H
