//
// Created by Yen-Chung Chen on 2/24/23.
//
#include <stdio.h>
#include <string.h>
#include <ctype.h> /* For tolower() */
#include "utils.h"
#include "sort.h"
#include "htslib/sam.h"
#include "shared_const.h" /* For shared constants */
#include "sys/stat.h" /* stat() and mkdir() */
#include <time.h>
#include <stdarg.h>
#include <errno.h>
#include <dirent.h>
#include <unistd.h>



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

void log_message(char* log_path, log_level_t out_level, char* message, log_level_t level, ...) {
    va_list args;
    va_start(args, level);

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
    fprintf(logf, "%s", timestamp);
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


int create_directory(char* pathname) {
    /**
     * @abstract Create a sustained directory. If a directory already exists, prompt the user for
     * confirmation.
     * @pathname A string containing relative path from the present working directory
     * @returns 0 on success, 1 if the user prompts NOT to overwrite, -1 if an error is encountered.
     */
   struct stat st = {0};

   // stat() returns 0 on success and -1 on failure
   // and put the file information into a statbuf
   if (stat(pathname, &st)) {
       int mkdir_stat = mkdir(pathname, 0700);
       if (0 != mkdir_stat) {
           int err_code = errno;
           log_msg("Fail to create directory (Error: %s)", ERROR, strerror(err_code));
           return -1;
       }
       return 0;
   } else {
       log_msg("Please note that output directory (%s) already exists.\n", WARNING, pathname);

       int checking;
       checking = 1;
       while (checking == 1) {
           fprintf(stderr, "Are you sure you want to save results and potentially OVERWRITE files there? [y/n]: ");
           unsigned char answer = (unsigned char) getc(stdin);



           // Clear the standard input to allow for the next round of input
           fflush(stdin);
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

char* create_tempdir(char *basedir) {
    char *tdir; // Name of temporary dir
    tdir = malloc(sizeof(char) * (strlen(basedir) + 5));
    strcpy(tdir, basedir);
    strcat(tdir, "tmp/");
    struct stat st = {0};
    if (stat(tdir, &st) == -1) {
        mkdir(tdir, 0700);
    }
    return tdir;
}

int8_t purge_tempdir(char *tmpdir) {
    //TODO: Need to rework!!! Do NOT use!!!
    DIR *dir_ptr;
    struct dirent *en;
    dir_ptr = opendir(tmpdir);

    char *rm_path;
    rm_path = (char *) malloc(sizeof(char) * (strlen(tmpdir) + 16));

    if (dir_ptr) {
        while((en = readdir(dir_ptr)) != NULL) {
            // Ignore . and ..
            if ('.' == en->d_name[0]) continue;

            // Construct path for individual file removal
            strcpy(rm_path, tmpdir);
            strcat(rm_path, en->d_name);
            int rm_stat = remove(rm_path);
            if (0 != rm_stat) {
                log_msg("Fail to remove temporary file %s", ERROR, en->d_name);
                return 1;
            }
        }
        free(rm_path);

        closedir(dir_ptr);
        int rmdir_stat = rmdir(tmpdir);
        if (0 != rmdir_stat) {
            log_msg("Fail to remove temporary directory (%s)", ERROR, tmpdir);
        }
        free(tmpdir);
    }
    return 0;
}

str_vec_t* str_vec_init (int64_t n, uint16_t str_len) {
    str_vec_t *svec;
    svec = (str_vec_t*) malloc(sizeof(str_vec_t));
    if (NULL == svec) return NULL;

    svec->length = n;
    svec->str_length = str_len;

    svec->str_arr = (char **) malloc(sizeof(char*) * n);

    for (int64_t i = 0; i < n; i++) {
        svec->str_arr[i] = malloc(sizeof(char) * str_len);
        if (NULL == svec->str_arr[i]) {
            for (int64_t j = 0; j < i; j++) {
                free(svec->str_arr[i]);
            }
            free(svec);
            return NULL;
        }
    }
    return svec;
}

int8_t str_vec_free(str_vec_t *ptr) {
    // Free elements of the str array
    for (uint32_t i = 0; i < ptr->length; i++) {
        free(ptr->str_arr[i]);
    }
    free(ptr->str_arr); // Free str array ptr
    free(ptr);
    return 0;
}

int64_t count_files(char *tmpdir) {
    DIR *dir_ptr;
    struct dirent *en;
    dir_ptr = opendir(tmpdir);

    if (dir_ptr) {
        int64_t file_count = 0;
        while(NULL != (en = readdir(dir_ptr))) {
            // Ignore . and ..
            if ('.' == en->d_name[0]) continue;
            file_count++;
        }
        closedir(dir_ptr);
        return file_count;
    }
    return -1;
}

str_vec_t * get_bams(char *tmpdir) {
    DIR *dir_ptr;
    struct dirent *en;
    dir_ptr = opendir(tmpdir);

    if (dir_ptr) {
        int64_t bam_count = count_files(tmpdir);
        if (-1 == bam_count) return NULL;

        str_vec_t *bam_list = str_vec_init(bam_count, 15);

        int8_t file_count = 0;
        while(NULL != (en = readdir(dir_ptr))) {
            // Ignore . and ..
            if ('.' == en->d_name[0]) continue;
            if (strlen(en->d_name) > (bam_list->str_length - 1)) return NULL;
            strcpy(bam_list->str_arr[file_count], en->d_name);
            file_count++;
        }

        closedir(dir_ptr);
        return bam_list;
    }
    return NULL;
}

char * tname_init(char * tmpdir, char * prefix, int32_t uid_length, uint32_t oid) {
    uint32_t pad_length = strlen(tmpdir) + strlen(prefix) + uid_length + 4; // 4: .bam
    char *name;
    name = (char *)malloc(sizeof(char) * pad_length + 1); // + null terminator
    char *uid;
    uid = (char *) malloc(sizeof(char) * uid_length + 4 + 1);
    char *fmt;
    fmt = (char *) malloc (sizeof (char) * 9);
    sprintf(fmt, "%%0%dd.bam", uid_length);
    strcpy(name, tmpdir);
    strcat(name, prefix);
    sprintf(uid, fmt, oid);
    strcat(name, uid);
    free(uid);
    free(fmt);
    return name;
}

int8_t merge_bam_pairs(char * tmpdir, char * f1, char * f2, uint32_t oid, char* prefix) {
    char *ff1;
    ff1 = (char *)malloc(sizeof(char) * (strlen(tmpdir) + strlen(f1) + 1));
    char *ff2;
    ff2 = (char *)malloc(sizeof(char) * (strlen(tmpdir) + strlen(f2) + 1));

    strcpy(ff1, tmpdir);
    strcat(ff1, f1);
    strcpy(ff2, tmpdir);
    strcat(ff2, f2);

    log_msg("Merging %s and %s", DEBUG, f1, f2);

    samFile *fp1 = sam_open(ff1, "r");
    sam_hdr_t *header1 = sam_hdr_read(fp1);
    samFile *fp2 = sam_open(ff2, "r");
    sam_hdr_t *header2 = sam_hdr_read(fp2);


    char *mname = tname_init(tmpdir, prefix, 5, oid);

    log_msg("Merging into %s", DEBUG, mname);

    htsFile* tfp = sam_open(mname, "wb");

    bam1_t *r1 = bam_init1();
    bam1_t *r2 = bam_init1();

    int r1stat = sam_read1(fp1, header1, r1);
    int r2stat = sam_read1(fp2, header2, r2);

    char * key1;
    char * key2;
    key1 = (char *) malloc(sizeof(char) * KEY_SIZE);
    key2 = (char *) malloc(sizeof(char) * KEY_SIZE);

    int wr_stat;
    bool f1done = false;
    bool f2done = false;

    //TODO: Figure out when headers need to be merged.
    int hdr_wstat = sam_hdr_write(tfp, header1);
    if (-1 == hdr_wstat) {
        log_msg("Fail to write header for merging BAM files (%s and %s)", ERROR, f1, f2);
        EARLY_EXIT_MERGE
        return 1;
    }

    while (r1stat >= 0 || r2stat >= 0) {
        int8_t tag_stat1 = get_tag(r1, "SK", key1);
        int8_t tag_stat2 = get_tag(r2, "SK", key2);

        if (!f1done && (f2done || strcmp(key1, key2) <= 0)) {
            wr_stat = sam_write1(tfp, header1, r1);
            r1stat = sam_read1(fp1, header1, r1);
            if (r1stat == -1) f1done = true;
        } else {
            wr_stat = sam_write1(tfp, header1, r2);
            r2stat = sam_read1(fp2, header1, r2);
            if (r2stat == -1) f2done = true;
        }

        if (wr_stat == -1) {
            log_msg("Fail to write merging reads into temporary files", ERROR);
            EARLY_EXIT_MERGE
            return 1;
        }
    }

    log_msg("Completed merging %s", DEBUG, mname);
    free(mname);
    free(key1);
    free(key2);
    bam_destroy1(r1);
    bam_destroy1(r2);
    sam_hdr_destroy(header1);
    sam_hdr_destroy(header2);
    sam_close(fp1);
    sam_close(fp2);
    sam_close(tfp);


    int rm1_stat = remove(ff1);
    int rm2_stat = remove(ff2);

    if (rm1_stat || rm2_stat) {
        log_msg("Fail to remove merged temporary files (%s and %s)", ERROR, ff1, ff2);
    }

    free(ff1);
    free(ff2);
    return 0;
}

int8_t merge_bams(char * tmpdir) {
    int64_t remaining = count_files(tmpdir);
    char *parray[] = {
            "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m",
            "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z"
    };

    if (-1 == remaining) {
        log_msg("Fail to examine remaining files in %s during merging", ERROR, tmpdir);
        return -1;
    }

    uint16_t mround = 0;
    while (remaining > 1) {
        str_vec_t *bam_vec = get_bams(tmpdir);

        // If there are an odds number of files, leave the last one to
        // the next round of merging
        uint32_t out_num = bam_vec->length / 2;
        for (uint32_t i = 0; i < out_num; i++) {
            merge_bam_pairs(tmpdir,
                            bam_vec->str_arr[i * 2],
                            bam_vec->str_arr[i * 2 + 1],
                            i, parray[mround % 26]);
        }
        str_vec_free(bam_vec);
        remaining = count_files(tmpdir);
        mround++;
    }
    str_vec_t *bam_vec = get_bams(tmpdir);

    char *sorted_path;
    char *sorted_out;
    // tmp/[sorted.bam]
    sorted_path = (char *) malloc(sizeof(char) * (strlen(tmpdir) + 11));
    sorted_out = (char *) malloc(sizeof(char) * (strlen(tmpdir) + bam_vec->str_length + 1));
    strcpy(sorted_out, tmpdir);
    strcat(sorted_out, bam_vec->str_arr[0]);
    strcpy(sorted_path, tmpdir);
    strcat(sorted_path, "sorted.bam");
    int mv_status = rename(sorted_out, sorted_path);
    free(sorted_path);
    free(sorted_out);
    str_vec_free(bam_vec);
    return 0;
}
