//
// Created by Yen-Chung Chen on 2/24/23.
//
#include <stdio.h>
#include <string.h>
#include <ctype.h> /* For tolower() */
#include "utils.h"
#include "htslib/sam.h"
#include "shared_const.h" /* For shared constants */
#include "sys/stat.h" /* stat() and mkdir() */
#include <time.h>
#include <stdarg.h>

int show_usage(const char* type) {
    if (strcmp(type, "unknown") == 0) {
        fprintf(stderr, "Unknown options are provided. Please see below for available options.\n\n");
    }
    if (strcmp(type, "missing") == 0) {
        fprintf(stderr, "A value must be provided for the following options:\n");
        fprintf(stderr, "-f/--file, -m/--meta, -q/--mapq, -o/--output, -t/--filter-tag, and -u/--umi-tag.\n\n");
    }
    fprintf(stderr, "Program: %s (Parallel BAM file subset by read tags)\n", APP_NAME);
    fprintf(stderr, "Version: %s (Dependent on htslib v%s)\n\n", VERSION, hts_version());
    fprintf(stderr, "Usage: scbamsplit [-f path] [-m path] [-o path] [-q MAPQ] [-t read-tag] [-u read-tag] [-d] [-h] [-n] [-v] \n\n");
    fprintf(stderr, "\t-f/--file: the path for input bam file\n");
    fprintf(stderr, "\t-m/--meta: the path for input metadata (an unquoted two-column csv with column names)\n");
    fprintf(stderr, "\t-o/--output: the path to export bam files to (default: ./)\n");
    fprintf(stderr, "\t-q/--mapq: Minimal MAPQ threshold for output (default: 0)\n");
    fprintf(stderr, "\t-t/--filter-tag: The read tag you want to filter against (default: CB)\n");
    fprintf(stderr, "\t-u/--umi-tag: The read tag containing UMIs (default: UB)\n");
    fprintf(stderr, "\t-d/--dedup: Remove duplicated reads with the same cell barcode/UMI combination\n");
    fprintf(stderr, "\t-n/--dry-run: Only print out parameters\n");
    fprintf(stderr, "\t-v/--verbose: Print out parameters\n");
    fprintf(stderr, "\t-h/--help: Show this documentation\n");
    return 0;
}

int create_folder(char* pathname) {
   struct stat st = {0};

   // stat() returns 0 on success and -1 on failure
   // and put the file information into a statbuf
   if (stat(pathname, &st)) {
       mkdir(pathname, 0700);
   } else {
       fprintf(stderr, "\n\nPlease note that output directory (%s) already exists.\n", pathname);
       char overwrite[3];
       char answer;

       int checking;
       checking = 1;
       while (checking == 1) {
           fprintf(stderr, "Are you sure you want to save results and potentially OVERWRITE files there? [y/n]: ");
           fgets(overwrite, sizeof(overwrite), stdin);

           // Clear the standard input to allow for the next round of input
           fflush(stdin);
           sscanf(overwrite, "%c\n", &answer);

           answer = tolower(answer);

           if (answer == 'n') {
               checking = 0;
               return 1;
           }

           if (answer != 'y'){
               fprintf(stderr, "Please only answer yes or no.\n");
           } else {
               checking = 0;
           }
       }

       }

   return 0;
}

char * get_time() {
    /**
     * Need to be freed!
     */
    time_t unixtime = time(NULL);
    struct tm *timeobj = localtime(&unixtime);

    // Expecting YYYY-MM-DD HH:MM:SS
    char *timestamp;
    timestamp = malloc(20 * sizeof(char));
    char fmt_str[8] = "";

    // Get year
    sprintf(fmt_str, "%d-", timeobj->tm_year + 1900);
    strcpy(timestamp, fmt_str);
    sprintf(fmt_str, "%02d-",timeobj->tm_mon + 1);
    strcat(timestamp, fmt_str);
    sprintf(fmt_str, "%02d ",timeobj->tm_mday);
    strcat(timestamp, fmt_str);
    sprintf(fmt_str, "%02d:",timeobj->tm_hour);
    strcat(timestamp, fmt_str);
    sprintf(fmt_str, "%02d:",timeobj->tm_min);
    strcat(timestamp, fmt_str);
    sprintf(fmt_str, "%02d",timeobj->tm_sec);
    strcat(timestamp, fmt_str);

    return timestamp;
}

void log_message(char* message, log_level level, char* log_path, log_level out_level, ...) {
    va_list args;
    va_start(args, out_level);

    FILE *logf;
    char *level_flag[] = {"[ERROR  ]", "[WARNING]", "", "[INFO   ]", "", "[DEBUG  ]"};
    if (strcmp("", log_path) != 0) {
        logf = fopen(log_path, "a");
    } else {
        logf = stderr;
    }

    // Only log when the level is met
    if (level > out_level) {return;}

    char * timestamp = get_time();

    // Tag the message
    fprintf(logf, "%s", level_flag[level]);

    // Timestamp the message
    fprintf(logf, " | ");
    fprintf(logf, timestamp);
    free(timestamp);

    // Print the message
    fprintf(logf, " | ");
    vfprintf(logf, message, args);
    va_end(args);
    fprintf(logf, "\n");

    if (strcmp("", log_path) != 0) {
        fclose(logf);
    }
}


