//
// Created by Yen-Chung Chen on 2/24/23.
//
#include "utils.h"
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include <ctype.h> /* For tolower() */
#include <time.h>
#include <stdarg.h>
#include <errno.h>
#include <dirent.h>
#include <unistd.h>
#include "htslib/sam.h"
#include "sys/stat.h" /* stat() and mkdir() */
#include "thread_pool.h"
#include "hash.h"


uint64_t max_strlen(char **strarr);

void show_usage() {
    fprintf(stderr, "Program: scbamsplit (Parallel BAM file subset by CBC/UMI)\n");
    fprintf(stderr, "Version: v0.3.1 (Dependent on htslib v1.17)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: scbamsplit -f path -m path\n");
    fprintf(stderr, "Options:\n\n");
    fprintf(stderr, "    Generic:\n");
    fprintf(stderr, "        [-o path] [-q MAPQ] [-d] [-r read name length] [-M memory usage (in GB)] [-n] [-v (verbosity)] [-h]\n");
    fprintf(stderr, "    CBC/UMI related:\n");
    fprintf(stderr, "        [-p platform] [-b CBC tag/field] [-L CBC length] [-u UMI tag/field] [-l UMI length]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -f/--file: the path for input bam file\n");
    fprintf(stderr, "    -m/--meta: the path for input metadata an unquoted two-column csv with column names)\n");
    fprintf(stderr, "    -o/--output: the path to export bam files to default: ./)\n");
    fprintf(stderr, "    -q/--mapq: Minimal MAPQ threshold for output default: 0)\n");
    fprintf(stderr, "    -p/--platform: Pre-fill locations and lengths for CBC and UMI (Supported platform: 10Xv2, 10Xv3, sciRNAseq3\n");
    fprintf(stderr, "         (e.g., for 10Xv3, both are stored as read tags. CBC is 16 mers tagged CB, while UMI is 12mers tagged UB).\n");
    fprintf(stderr, "    -d/--dedup: Remove duplicated reads with the same cell barcode/UMI combination\n");
    fprintf(stderr, "    -b/--cbc-location: If CBC is a read tag, provide the name (e.g., CB); if it is in the read name,\n");
    fprintf(stderr, "         provide the field number (e.g., 3) (default: CB)\n");
    fprintf(stderr, "    -L/--cbc-length: The length of the barcode you want to filter against (default: 20)\n");
    fprintf(stderr, "    -u/--umi-location: If UMI is a read tag, provide the name (e.g., UB); if it is in the read name,\n");
    fprintf(stderr, "        provide the field number (e.g., 3) (default: CB)\n");
    fprintf(stderr, "    -l/--umi-length: The length of the UMI default: 20)\n");
    fprintf(stderr, "    -r/--rn-length: The length of the read name (default: 70)\n");
    fprintf(stderr, "    -M/--mem: The estimated maximum amount of memory to use (In GB, default: 4)\n");
    fprintf(stderr, "    -@/--threads: Setting the number of threads to use (default: 1)\n");
    fprintf(stderr, "    -n/--dry-run: Only print out parameters\n");
    fprintf(stderr, "    -v/--verbose: Set verbosity level (1 - 5) (default: 2, 3 if -v provided without a value)\n");
    fprintf(stderr, "    -h/--help: Show this documentation\n");
}

char * get_time() {
    /**
     * Need to be freed!
     */
    time_t unixtime = time(NULL);
    struct tm *timeobj = localtime(&unixtime);

    // Expecting YYYY-MM-DD HH:MM:SS
    char *timestamp;
    timestamp = calloc(20, sizeof(char));
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
    if (strcmp("", log_path) != 0) {
        logf = fopen(log_path, "a");
    } else {
        logf = stderr;
    }

    // Only log when the level is met
    if (level > out_level) {return;}

    char * timestamp = get_time();

    // Tag the message
    char fmt[16] = "";
    sprintf(fmt, "[%%-%llus]", max_strlen(LEVEL_FLAG) - 1);
    fprintf(logf, fmt, LEVEL_FLAG[level]);

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

uint64_t max_strlen(char **strarr) {
    uint64_t tag_width = 0;
    if (NULL == strarr) {
        goto early_exit;
    }
    for (int8_t i = 0; i < 6; i++) {
        if (strlen(LEVEL_FLAG[i]) > tag_width) {
            tag_width = strlen(strarr[i]);
        }
    }
    tag_width = tag_width + 1; // NULL terminator
    early_exit:
    return tag_width;
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
    tdir = calloc((strlen(basedir) + 5), sizeof(char));
    strcpy(tdir, basedir);
    strcat(tdir, "tmp/");
    struct stat st = {0};
    if (stat(tdir, &st) == -1) {
        mkdir(tdir, 0700);
    }
    return tdir;
}

str_vec_t* str_vec_init (int64_t n, uint16_t str_len) {
    str_vec_t *svec;
    svec = (str_vec_t*) calloc(1, sizeof(str_vec_t));
    if (NULL == svec) return NULL;

    svec->length = n;
    svec->str_length = str_len;

    svec->str_arr = (char **) calloc(n, sizeof(char*));

    for (int64_t i = 0; i < n; i++) {
        svec->str_arr[i] = calloc(str_len + 1, sizeof(char));
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

str_vec_t *str_vec_copy(str_vec_t *ptr, int32_t from) {
    int32_t farray_size = ptr->length - from;
    char **farray_copy;
    farray_copy = calloc(farray_size, sizeof(char*));
    for (int32_t i = 0;  i < farray_size; i++) {
        farray_copy[i] = calloc(ptr->str_length , sizeof(char));
        strcpy(farray_copy[i], ptr->str_arr[from + i]);
    }
    str_vec_t* svec;
    svec = calloc(1, sizeof(str_vec_t));
    svec->length = farray_size;
    svec->str_length = ptr->str_length;
    svec->str_arr = farray_copy;
    return svec;
}

int8_t str_vec_destroy(str_vec_t *ptr) {
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

        str_vec_t *bam_list = str_vec_init(bam_count, 18);

        int64_t file_count = 0;
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
    name = (char *)calloc(pad_length + 1, sizeof(char)); // + null terminator
    char *uid;
    uid = (char *) calloc(uid_length + 4 + 1, sizeof(char)) ;
    char *fmt;
    fmt = (char *) calloc (9, sizeof (char));
    sprintf(fmt, "%%0%dd.bam", uid_length);
    strcpy(name, tmpdir);
    strcat(name, prefix);
    sprintf(uid, fmt, oid);
    strcat(name, uid);
    free(uid);
    free(fmt);
    return name;
}

int8_t merge_bam_nway(char *tmpdir, str_vec_t *bam_vec, uint32_t oid, char *prefix, int64_t n) {
    char **farray = bam_vec->str_arr;
    char **ffarray = calloc(n, sizeof(char*));
    samFile **fpa = calloc(n, sizeof(samFile*));
    sam_hdr_t **header_arr = calloc(n, sizeof(sam_hdr_t*));
    bam1_t **rarray = calloc(n, sizeof(bam1_t*));
    int32_t *rstat_arr = calloc(n, sizeof(int32_t));
    char **key_arr = calloc(n, sizeof(char*));
    char *key_temp = calloc(KEY_SIZE, sizeof(char));
    int8_t *tag_stat_arr = calloc(n, sizeof(int8_t));


    for (int64_t i = 0; i < n; i++) {
        ffarray[i] = calloc(strlen(tmpdir) + strlen(farray[i]) + 1, sizeof(char));
        strcpy(ffarray[i], tmpdir);
        strcat(ffarray[i], farray[i]);

        fpa[i] = sam_open(ffarray[i], "r");
        header_arr[i] = sam_hdr_read(fpa[i]);
        rarray[i] = bam_init1();
        rstat_arr[i] = sam_read1(fpa[i], header_arr[i], rarray[i]);
        key_arr[i] = calloc(KEY_SIZE, sizeof(char));
    }

    char *mname = tname_init(tmpdir, prefix, 5, oid);
    htsFile* tfp = sam_open(mname, "wb");

    int wr_stat;
    int8_t return_val = 0;

    //TODO: Figure out when headers need to be merged.
    int hdr_wstat = sam_hdr_write(tfp, header_arr[0]);
    if (-1 == hdr_wstat) {
        log_msg("Fail to write header for merging BAM files", ERROR);
        return_val = 1;
        goto release_key_and_fp_and_exit;
    }

    int32_t any_r = 1;
    bool first_item = true;
    int64_t key_min_id;
    while (any_r == 1) {
        first_item = true;
        key_min_id = 0;
        for (int8_t i = 0; i < n; i++) {
            if (rstat_arr[i] < 0) continue;
            tag_stat_arr[i] = fetch_tag(rarray[i], "SK", key_arr[i]);

            if (first_item) {
                strcpy(key_temp, key_arr[i]);
                first_item = false;
                key_min_id = i;
                continue;
            }

            if (strcmp(key_temp, key_arr[i]) > 0) {
                strcpy(key_temp, key_arr[i]);
                key_min_id = i;
            }
        }

        wr_stat = sam_write1(tfp, header_arr[0], rarray[key_min_id]);
        rstat_arr[key_min_id] = sam_read1(fpa[key_min_id], header_arr[key_min_id], rarray[key_min_id]);

        for (int64_t i = 0; i < n; i++) {
            if (rstat_arr[i] >= 0) {
                any_r = 1;
                break;
            }
            any_r = 0;
        }

        if (wr_stat == -1) {
            log_msg("Fail to write merging reads into temporary files", ERROR);
            return_val = 1;
            goto release_key_and_fp_and_exit;
        }
    }

    log_msg("Completed merging %s", DEBUG, mname);
    bool rm_err = false;
    release_key_and_fp_and_exit:
    free(mname);
    for (int64_t i = 0; i < n; i++) {
        sam_hdr_destroy(header_arr[i]);
        bam_destroy1(rarray[i]);
        sam_close(fpa[i]);
        free(key_arr[i]);
        if (unlink(ffarray[i]) != 0) rm_err = true;
        free(ffarray[i]);
        if (rm_err) {
            log_msg("Fail to remove merged temporary files ", ERROR);
        }
    }
    free(rstat_arr);
    free(ffarray);
    free(fpa);
    free(header_arr);
    free(rarray);
    free(key_arr);
    free(key_temp);
    free(tag_stat_arr);
    str_vec_destroy(bam_vec);
    sam_close(tfp);

    return return_val;
}

void pmerge_bam_nway(void *args) {
    mnway_args *cargs = (mnway_args*) args;
    merge_bam_nway(cargs->tmpdir, cargs->bam_vec, cargs->oid, cargs->prefix, cargs->n);
}

char* merge_bams(char * tmpdir) {
    int64_t remaining = count_files(tmpdir);
    char *parray[] = {
            "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m",
            "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z"
    };

    if (-1 == remaining) {
        log_msg("Fail to examine remaining files in %s during merging", ERROR, tmpdir);
        return NULL;
    }

    uint16_t mround = 0;
    char * long_prefix;
    // Make a prefix in the form of merged[a-z]00001.bam
    // so the original chunkXXXXX.bam will be merged first and
    // round a will be merged first in round b if left to the next round
    // in a previous round with odd # of members
    long_prefix = calloc(9, sizeof(char));

    tpool_t *merge_tp = tpool_create(MAX_THREADS, MAX_THREADS);
    while (remaining > 1) {
        str_vec_t *bam_vec = get_bams(tmpdir);
        strcpy(long_prefix, "merged");
        strcat(long_prefix, parray[mround % 26]);

        uint32_t batch_size = 8;
        uint32_t processed_size = 0;
        uint32_t oid = 0;
        uint32_t residual_size = 0;

        if (bam_vec->length <= batch_size) {
            str_vec_t *bam_vec_copy = str_vec_copy(bam_vec, 0);
            merge_bam_nway(tmpdir, bam_vec_copy, 1, long_prefix, bam_vec->length);
        } else {
            while (processed_size < bam_vec->length) {
                residual_size = bam_vec->length - processed_size;
                if (residual_size < batch_size) batch_size = residual_size;

                if (MAX_THREADS == 1) {
                    str_vec_t *bam_vec_copy = str_vec_copy(bam_vec, processed_size);
                    merge_bam_nway(tmpdir, bam_vec_copy, oid, long_prefix, batch_size);
                } else {
                    str_vec_t *bam_vec_copy = str_vec_copy(bam_vec, processed_size);
                    mnway_args args = {
                            .tmpdir = tmpdir,
                            .bam_vec = bam_vec_copy,
                            .oid = oid,
                            .prefix = long_prefix,
                            .n = batch_size,
                    };
                    bool tqueue = tpool_add_work(merge_tp, pmerge_bam_nway, &args, sizeof(mnway_args));
                }
                oid++;
                processed_size += batch_size;
            }
        }
        str_vec_destroy(bam_vec);
        tpool_wait(merge_tp);
        remaining = count_files(tmpdir);
        mround++;
    }
    tpool_wait(merge_tp);
    tpool_destroy(merge_tp);

    free(long_prefix);
    str_vec_t *bam_vec = get_bams(tmpdir);

    char *sorted_path;
    char *sorted_out;
    // tmp/[sorted.bam]
    sorted_path = (char *) calloc(strlen(tmpdir) + 11, sizeof(char));
    sorted_out = (char *) calloc(strlen(tmpdir) + bam_vec->str_length + 1, sizeof(char));
    strcpy(sorted_out, tmpdir);
    strcat(sorted_out, bam_vec->str_arr[0]);
    strcpy(sorted_path, tmpdir);
    strcat(sorted_path, "sorted.bam");
    int mv_status = rename(sorted_out, sorted_path);
    free(sorted_out);
    str_vec_destroy(bam_vec);
    return sorted_path;
}

int8_t read_dump(rt2label *r2l, rt2label *lout, label2fp *l2fp, label2fp *fout,
                 char * this_CB, sam_hdr_t *header, bam1_t *read) {
    int32_t write_to_bam = 0;
    HASH_FIND_STR(r2l, (char *) this_CB, lout);

    // If the CBC is found
    if (lout) {
        // Query the CBC-to-output table
        HASH_FIND_STR(l2fp, lout->label, fout);
        if (fout) {
            write_to_bam = sam_write1(fout->fp, header, read);
            if (write_to_bam < 0) {
                // Decide how to deal with writing failure outside
                return 1;
            }
        }
    }
    return 0;
}

int8_t deduped_dump(rt2label *r2l, rt2label *lout, label2fp *l2fp, label2fp *fout, char *tmpdir, char *sorted_path,
                    bam1_t *read, char *bc_tag, char *umi_tag, tag_meta_t *cb_meta, tag_meta_t *ub_meta) {
    int32_t read_stat;
    samFile *sfp = sam_open(sorted_path, "r");
    sam_hdr_t *sheader = sam_hdr_read(sfp);
    char *current_UB;
    char *this_UB;
    char *RN_keep;
    char *this_RN;
    char *current_CB;
    char *this_CB;

    current_UB = (char *) calloc(UB_LENGTH, sizeof(char));
    this_UB = (char *) calloc(UB_LENGTH, sizeof(char));
    RN_keep = (char *) calloc(RN_SIZE, sizeof(char));
    this_RN = (char *) calloc(RN_SIZE, sizeof(char));
    current_CB = (char *) calloc(CB_LENGTH, sizeof(char));
    this_CB = (char *) calloc(CB_LENGTH, sizeof(char));


    bool first_read = true;
    int32_t same_CB = 0;
    int32_t same_UB = 0;
    int32_t to_export = 0;
    int32_t write_to_bam = 0;
    int8_t return_val = 0;
    while (0 <= (read_stat = sam_read1(sfp, sheader, read))) {
        // Get read metadata
        int8_t cb_stat = get_CB(read, cb_meta, this_CB);
        if (-1 == cb_stat) {
            return_val = 1;
            log_msg("Cannot retrieve cell barcode from the sorted BAM", ERROR);
            goto free_res_and_exit;
        }
        int8_t ub_stat = get_UB(read, ub_meta, this_UB);
        if (-1 == ub_stat) {
            return_val = 1;
            log_msg("Cannot retrieve UMI from the sorted BAM", ERROR);
            goto free_res_and_exit;
        }
        char * rn_ptr;
        if (NULL == (rn_ptr = bam_get_qname(read))) {
            log_msg("Cannot retrieve read name from the sorted BAM", ERROR);
        } else {
            strcpy(this_RN, rn_ptr);
        }

        // We want to keep the first primary read of each CB-UMI combination
        // and all secondary mappings of the same read.
        // With the sorting mechanism, this will be the first read of each
        // CB-UMI combo.
        if (first_read) {
            first_read = !first_read;
            strcpy(current_CB, this_CB);
            strcpy(current_UB, this_UB);
            strcpy(RN_keep, this_RN);
        }


        // Check if we have entered the next CB-UMI combo.
        // Update the RN to keep if so.
        same_CB = strcmp(current_CB, this_CB);
        same_UB = strcmp(current_UB, this_UB);
        if (0 != same_CB || 0 != same_UB) {
            strcpy(RN_keep, this_RN);
        }
        if (0 != same_CB) {
            strcpy(current_CB, this_CB);
        }
        if (0 != same_UB) {
            strcpy(current_UB, this_UB);
        }

        // Export reads with the highest MAPQ per CB-UMI combo
        to_export = strcmp(RN_keep, this_RN);
        if (0 != to_export) {
            continue;
        }

        // Exporting process
        int8_t rdump_stat = read_dump(r2l, lout, l2fp, fout, this_CB, sheader, read);
        if (0 != rdump_stat) {
            return_val = 1;
            log_msg("Fail to write sorted reads to split BAM file (%s)", ERROR, lout->label);
            goto free_res_and_exit;
        }

    }

    int32_t rm_tmp = remove(sorted_path);
    if (0 != rm_tmp) {
        return_val = 1;
        log_msg("Fail to remove temporary file (%s)", ERROR, sorted_path);
        goto free_res_and_exit;
    }

    int32_t rm_tmpdir = rmdir(tmpdir);
    if (0 != rm_tmpdir) {
        return_val = 1;
        log_msg("Fail to remove temporary (%directory)", ERROR, tmpdir);
    }

    free_res_and_exit:
        free(current_UB);
        free(this_UB);
        free(RN_keep);
        free(this_RN);
        free(current_CB);
        free(this_CB);
        sam_close(sfp);
        sam_hdr_destroy(sheader);
        free(sorted_path);
        free(tmpdir);
    return return_val;
}
