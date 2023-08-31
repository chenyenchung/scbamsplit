//
// Created by Yen-Chung Chen on 2/24/23.
//

#ifndef SCBAMSPLIT_UTILS_H
#define SCBAMSPLIT_UTILS_H
int show_usage(const char* type);
int create_folder(char* pathname);
char * create_tempdir(char *dir);
typedef enum {
    DEBUG = 5,
    ERROR = 0,
    WARNING = 1,
    INFO = 3
} log_level;

void log_message(char* message, log_level level, char* log_path, log_level out_level, ...);
#endif //SCBAMSPLIT_UTILS_H
