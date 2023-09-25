//
// Created by Yen-Chung Chen on 8/23/23.
//
#include <string.h>
#include <stdbool.h>
#include <ctype.h>
#include "htslib/kstring.h"
#include "sort.h"
#include "utils.h"
#include "thread_pool.h"

int8_t fetch_tag(bam1_t *read, char *tag, char* tag_ptr) {
    // bam_aux_get_str() exists, but to return a kstring,
    // it involves realloc() and is rather slow.
    // We ask users to specify UB and CB size here, so we
    // should take advantage of this information to gain
    // some efficiency.
    char* tag_content = (char *) bam_aux_get(read, tag);
    if (NULL == tag_content) {
        // Tag not found
        return -1;
    }
    strcpy(tag_ptr, tag_content + 1);
    return 0;
}
int8_t fetch_tag2(bam1_t *read, char* tag_ptr, tag_meta_t *info) {
    // bam_aux_get_str() exists, but to return a kstring,
    // it involves realloc() and is rather slow.
    // We ask users to specify UB and CB size here, so we
    // should take advantage of this information to gain
    // some efficiency.
    char* tag_content = (char *) bam_aux_get(read, info->tag_name);
    if (NULL == tag_content) {
        // Tag not found
        return -1;
    }
    strncpy(tag_ptr, tag_content + 1, info->length - 1);
    return 0;
}

int8_t fetch_name(bam1_t *read, char* tag_ptr, tag_meta_t *info) {
    char* rn = bam_get_qname(read);

    // Need to copy -- strtok will modify the string
    char *rncpy = calloc(strlen(rn) + 1, sizeof(char));
    strcpy(rncpy, rn);
    char *token;


    uint8_t field_num = 0;
    int8_t return_val = 1;
    token = strtok(rncpy, info->sep);
    while (token != NULL) {
        field_num++;
        if (field_num == info->field) {
            strcpy(tag_ptr, token);
            return_val = 0;
        }
        token = strtok(NULL, info->sep);
    }
    free(rncpy);
    return return_val;
}
int8_t get_CB (bam1_t *read, tag_meta_t* info, char* tag_ptr) {
    int8_t exit_code = 0;
    switch (info->location) {
        case READ_TAG:
            fetch_tag2(read, tag_ptr, info);
            break;
        case READ_NAME:
            fetch_name(read, tag_ptr, info);
            break;
        default:
            log_msg("Unknown location type to fetch cell barcode", ERROR);
            exit_code = 1;
    }
    if (NULL == tag_ptr) {
        exit_code = -1;
    }
    return exit_code;
}
int8_t get_UB (bam1_t *read, tag_meta_t* info, char* tag_ptr) {
    int8_t exit_code = 0;
    switch (info->location) {
        case READ_TAG:
            fetch_tag2(read, tag_ptr, info);
            break;
        case READ_NAME:
            fetch_name(read, tag_ptr, info);
            break;
        default:
            log_msg("Unknown location type to fetch cell barcode", ERROR);
            exit_code = 1;
    }
    if (NULL == tag_ptr) {
        exit_code = -1;
    }
    return exit_code;
}


tag_meta_t *initialize_tag_meta() {
    // Allocate
    tag_meta_t *tag_meta;
    tag_meta = calloc(1, sizeof(tag_meta_t));
    tag_meta->tag_name = calloc(3, sizeof(char));
    tag_meta->sep = calloc(2, sizeof(char));

    // Set initial values
    tag_meta->location = READ_TAG;
    strcpy(tag_meta->tag_name, "CB");
    strcpy(tag_meta->sep, ",");
    tag_meta->field = 1;
    tag_meta->length = 21;

    return tag_meta;
}

void destroy_tag_meta(tag_meta_t *tag_meta) {
    free(tag_meta->sep);
    free(tag_meta->tag_name);
    free(tag_meta);
}
void print_tag_meta(tag_meta_t *tag_meta, const char *header) {
    char* location[2] = {"Read tag", "Read name"};
    if (NULL == header) {
        fprintf(stderr, "Tag Information\n");
    } else {
        fprintf(stderr, "\t%s:\n", header);
    }
    fprintf(stderr, "\t\tLocation: %s\n", location[tag_meta->location]);
    if (tag_meta->location == READ_NAME) {
        fprintf(stderr, "\t\tSep char: %s\n", tag_meta->sep);
        fprintf(stderr, "\t\tField number: %d\n", tag_meta->field);
    } else {
        fprintf(stderr, "\t\tTag name: %s\n", tag_meta->tag_name);
    }
    fprintf(stderr, "\t\tTag length: %d\n\n", tag_meta->length - 1);
}

void set_CB(tag_meta_t *tag_meta, char *platform) {
    for (uint32_t i = 0; i < strlen(platform); i++) {
        platform[i] = tolower(platform[i]);
    }
    if (strcmp("10xv2", platform) == 0) {
        tag_meta->length = 18 + 1;
    } else if (strcmp("scirnaseq3", platform) == 0) {
        tag_meta->location = READ_NAME;
        tag_meta->length = 20 + 1;
        tag_meta->field = 1;
    } else {
        // Default (platform == NULL) is 10Xv3
        tag_meta->length = 18 + 1;
    }
}

void set_UB(tag_meta_t *tag_meta, char *platform) {
    for (uint32_t i = 0; i < strlen(platform); i++) {
        platform[i] = tolower(platform[i]);
    }
    if (strcmp("10xv2", platform) == 0) {
        strcpy(tag_meta->tag_name, "UB");
        tag_meta->length = 10 + 1;
    } else if (strcmp("scirnaseq3", platform) == 0) {
        tag_meta->location = READ_NAME;
        tag_meta->length = 8 + 1;
        tag_meta->field = 2;
    } else {
        // Default (platform == NULL) is 10Xv3
        strcpy(tag_meta->tag_name, "UB");
        tag_meta->length = 12 + 1;
    }
}

int8_t is_primary(bam1_t *read, char* flag_ptr) {
    /**
     * @abstract Examine bit flag of a SAM read to see if it's a mapped
     * primary record
     * @read A pointer to bam1_t (HTSlib)
     * @returns 1 if the record is primary; 2 in other cases
     */
    bool not_secondary = !((read)->core.flag & 256); // Flag for secondary alignment is 256
    bool mapped = !((read)->core.flag & 4); // Flag for unmapped is 4

    // This is used for sorting, so we want primary to be coded 1 while others being 2
    sprintf(flag_ptr, "%d", 2 - not_secondary && mapped);
    return 0;
}

void get_MAPQ(bam1_t *read, char* val_ptr) {
    uint8_t mapq = (read)->core.qual;
    sprintf(val_ptr, "%03d", mapq);
}

int64_t fill_chunk(samFile *fp, sam_hdr_t *header, ichunk_t *ic, int16_t qthres,
                   tag_meta_t *cb_meta, tag_meta_t *ub_meta) {
    /**
     * @abstract Fill read buffer to designated size and return the index of next read to read or -1
     * when fails.
     *
     * @fp A file descriptor for a SAM/BAM file
     * @header A header pointer from sam_hdr_read()
     * @read A pointer to a sam_read_t array to be filled
     * @chunk_size An integer to indicate how large the cache chunk to be filled
     * @qthres An integer specifying the MAPQ threshold to pass to keep the read
     * @returns The number of reads that have been allocated into the chunk on success; -1 on error;
     * -[OBSERVED_READ_NAME_SIZE] when read names are not sufficiently padded
     */

    sam_read_t **read_array = ic->chunk;
    int64_t chunk_size = ic->chunk_size;

    // Note that read_array must be allocated OUTSIDE!
    bam1_t *temp_read = bam_init1();
    int64_t read_kept = 0;
    char *CB;
    char *UB;
    char *MAPQ;
    char *PR;
    char *RN; // Declare empty string of sufficient size
    CB = (char *) calloc(CB_LENGTH, sizeof(char));
    UB = (char *) calloc(UB_LENGTH, sizeof(char));
    MAPQ = (char *) calloc(4, sizeof(char)); // Max value for MAPQ is 255 per SAM spec v1
    PR = (char *) calloc(2, sizeof(char)); // The value for primary read is either 0 or 1
    RN = (char *) calloc(RN_SIZE, sizeof(char));

    // Fill the chunk until specified size or running out of reads
    while (read_kept < chunk_size) {
        if (sam_read1(fp, header, temp_read) < 0) {
            // sam_read1 returns -1 when encountering an error
            // or EOF
            if (read_kept > 0) {
                read_kept--;
                goto stop_fill_and_free;
            }
            read_kept = -1;
            goto stop_fill_and_free;
        }

        int8_t cb_stat = get_CB(temp_read, cb_meta, CB);
        int8_t ub_stat = get_UB(temp_read, ub_meta, UB);
        int8_t prim_stat = is_primary(temp_read, PR);
        get_MAPQ(temp_read, MAPQ); // MAPQ seems to be guaranteed by SAM spec
        int16_t mapq_val = temp_read->core.qual;

        // Skip reads that miss CB, UB, bitwise flag, or MAPQ
        if ((-1 == cb_stat) || (-1 == ub_stat) || (1 == prim_stat) || qthres > mapq_val) {
            // If a read has no CBC or UMI, just let it go.
            continue;
        } else if (1 == cb_stat || 1 == ub_stat) {
            // Error message is produced in fetch_tag()
            read_kept = -1;
            goto stop_fill_and_free;
        }

        bam1_t *read_copy_status = bam_copy1(read_array[read_kept]->read, temp_read);
        if (NULL == read_copy_status) {
            log_msg( "Fail to copy a BAM read. Could be an HTSlib issue?", ERROR);
            read_kept = -1;
            goto stop_fill_and_free;
        }

        // Prepare format string to pad the read name
        char rn_length[8];
        char fmt[16];

        // Get the read name to be padded to declared size - 1
        sprintf(rn_length, "%lld", RN_SIZE - 1);

        // Prepare format string for sprintf padding
        strcpy(fmt, "%-");
        strcat(fmt, rn_length);
        strcat(fmt, "s");

        // Pad read name
        int64_t obs_rn_size = strlen(bam_get_qname(temp_read));
        if (obs_rn_size > (RN_SIZE - 1)) {
            read_kept = -obs_rn_size;
            goto stop_fill_and_free;
        }
        sprintf(RN, fmt, bam_get_qname(temp_read));

        // Initialize the key string
        strcpy(read_array[read_kept]->key, "");

        // Concatenate CBC, UMI, read name, primary or secondary, MAPQ
        // Lump reads that are supposed to be the same molecule
        strcat(read_array[read_kept]->key, CB);
        strcat(read_array[read_kept]->key, UB);
        // Show primary mapping first and sort by MAPQ
        // so later exporting mechanism can just export the
        // first and all secondary mappings (if exists) by
        // read name
        strcat(read_array[read_kept]->key, PR);
        strcat(read_array[read_kept]->key, MAPQ);
        strcat(read_array[read_kept]->key, RN);

        // Append sorting key in the output BAM as tag "SK"
        int app_stat = bam_aux_append(
                read_array[read_kept]->read, "SK", 'Z',
                sizeof(char) * (strlen(read_array[read_kept]->key) + 1),
                (const uint8_t*) read_array[read_kept]->key
                );
        if (0 != app_stat) {
            log_msg("Fail to append sorting key to the read (%s)",
            WARNING, read_array[read_kept]->key);
        }
        read_kept++;
    }

    stop_fill_and_free:
        ic->read_kept = read_kept;
        bam_destroy1(temp_read);
        free(PR);
        free(MAPQ);
        free(UB);
        free(CB);
        free(RN);
    return read_kept;
}

int read_cmp(const void *a, const void *b) {
    const sam_read_t * reada = *(sam_read_t **) a;
    const sam_read_t * readb = *(sam_read_t **) b;

    return strcmp(reada->key, readb->key);
}

void sort_chunk(ichunk_t *ic) {
    sam_read_t **reads = ic->chunk;
    int64_t chunk_size = ic->read_kept;
    qsort(reads, chunk_size, sizeof(sam_read_t*), read_cmp);
}

sam_read_t** chunk_init(uint32_t chunk_size) {
    sam_read_t **chunk;
    chunk = calloc(chunk_size, sizeof(sam_read_t*));

    if (NULL == chunk) {
        log_msg("Fail to allocate memory for sorting", ERROR);
        goto exit_without_free;
    }

    uint32_t failed_iter = 0;
    for (uint32_t read_i = 0; read_i < chunk_size; read_i++) {
        chunk[read_i] = calloc(1, sizeof(sam_read_t));
        chunk[read_i]->read = bam_init1();
        chunk[read_i]->key = calloc(KEY_SIZE, sizeof(char));
        if (NULL == chunk[read_i] || NULL == chunk[read_i]->read || NULL == chunk[read_i]->key) {
            failed_iter = read_i;
            goto free_read_ptr;
        }
    }

    // A HTSlib-defined BAM read contains:
    // A bam1_t struct, which also points to a bam1_core_t
    return chunk;

    free_read_ptr:
        chunk_destroy(chunk, failed_iter);
    exit_without_free:
        return NULL;
}

void chunk_destroy(sam_read_t **read_array, uint32_t chunk_size) {
    for (uint32_t read_i = 0; read_i < chunk_size; read_i++) {
        bam_destroy1(read_array[read_i]->read);
        free(read_array[read_i]->key);
        free(read_array[read_i]);
    }
    free(read_array);
}

void sort_export_chunk(void *args_void) {
    chunk_arg_t *args =(chunk_arg_t *) args_void;
    ichunk_t* ic = args->ic;

    sort_chunk(ic);

    char *tname = tname_init(args->tmpdir, "chunk", 5, args->tid);

    htsFile* tfp = sam_open(tname, "wb");

    int write_status = sam_hdr_write(tfp, args->header);
    if (write_status != 0) {
        log_msg( "Fail to write SAM header into temporary files", ERROR);
        goto free_and_exit;
    }

    int write_to_bam = 0;
    for (int64_t i = 0; i < ic->read_kept; i++) {
        sam_read_t* cr = ic->chunk[i];
        write_to_bam = sam_write1(tfp, args->header, cr->read);
        if (write_to_bam == -1) {
            log_msg("Fail to write sorted reads into temporary files", ERROR);
            goto free_and_exit;
        }
    }

    chunkq_add(args->give_q, ic);

    free_and_exit:
    free(tname);
    sam_close(tfp);
}

char *process_bam(samFile *fp, sam_hdr_t *header, int64_t chunk_size, char *oprefix, int64_t qthres,
                  tag_meta_t *cb_meta, tag_meta_t *ub_meta) {
    /**
     * @abstract Process all reads in an opened SAM/BAM file in chunks and save sorted reads in a temporary
     * directory.
     * @fp A pointer to a samFile
     * @header A pointer to a SAM file header
     * @chunk An array of sam_read_t's
     * @chunk_size The size of the array (a 64-bit integer)
     * @oprefix Output dir prefix
     * @qthres An integer specifying the MAPQ threshold to pass to keep the read
     * @returns A string: the path for the temporary directory containing sorted chunks if succeeded; "-1" if failed.
     */

    int64_t size_retrieved = chunk_size;
    int32_t chunk_num = 0;
    // Create temporary file dir for sorted chunks
    char *tmpdir = create_tempdir(oprefix);

    chunkq_t *init_q = chunkq_create();
    chunkq_t *chunk_q = chunkq_create();
    tpool_t *sort_tp;
    if (MAX_THREADS == 1) {
        chunkq_ele_init(init_q, chunk_size);
    } else {
        sort_tp = tpool_create(MAX_THREADS, MAX_THREADS);
        for (int32_t i = 0; i < MAX_THREADS; i++) {
            chunkq_ele_init(init_q, chunk_size);
        }
    }

    while (size_retrieved == chunk_size) {
        chunk_num += 1;

        ichunk_t *this_chunk;
        log_msg("Chunk #%lld filling", DEBUG, chunk_num);
        if (init_q->qlength > 0) {
            this_chunk = chunkq_get(init_q);
        } else {
            this_chunk = chunkq_get(chunk_q);
        }

        size_retrieved = fill_chunk(fp, header, this_chunk, qthres, cb_meta, ub_meta);

        if (size_retrieved < -1) {
            log_msg("Insufficient RN size (%d).", ERROR, RN_SIZE - 1);
            log_msg("Please increase RN size (-r/--rn-length) to at least %d", ERROR, -size_retrieved + 1);
            goto free_tmpdir_exit;
        } else if (size_retrieved == -1) {
            // All reads are fetched
            goto wait_for_threads;
        }
        chunkq_add(chunk_q, this_chunk);

        if (MAX_THREADS == 1) {
            ichunk_t *this_chunk = chunkq_get(chunk_q);
            sort_chunk(this_chunk);

            char *tname = tname_init(tmpdir, "chunk", 5, chunk_num);
            htsFile* tfp = sam_open(tname, "wb");
            free(tname);

            int write_status = sam_hdr_write(tfp, header);
            if (write_status != 0) {
                log_msg( "Fail to write SAM header into temporary files", ERROR);
                sam_close(tfp); // This has to be done in the while loop
                goto free_tmpdir_exit;
            }

            int write_to_bam = 0;
            for (int64_t i = 0; i < size_retrieved; i++) {
                write_to_bam = sam_write1(tfp, header, (this_chunk->chunk)[i]->read);
                if (write_to_bam == -1) {
                    log_msg("Fail to write sorted reads into temporary files", ERROR);
                    sam_close(tfp); // This has to be done in the while loop
                    goto free_tmpdir_exit;
                }
            }
            sam_close(tfp);
            chunkq_add(chunk_q, this_chunk);
        } else {
            chunk_arg_t args = {
                    .ic = chunkq_get(chunk_q),
                    .give_q = chunk_q,
                    .header = header,
                    .tmpdir = tmpdir,
                    .tid = chunk_num,
            };
            tpool_add_work(sort_tp, sort_export_chunk, &args, sizeof(chunk_arg_t));
        }
    }

    wait_for_threads:
    if (MAX_THREADS > 1) {
        tpool_wait(sort_tp);
        tpool_destroy(sort_tp);
    }
    chunkq_destroy(init_q);
    chunkq_destroy(chunk_q);

    return tmpdir;

    free_tmpdir_exit:
    if (MAX_THREADS > 1) {
        tpool_wait(sort_tp);
        tpool_destroy(sort_tp);
    }
    free(tmpdir);
    chunkq_destroy(init_q);
    chunkq_destroy(chunk_q);
    return "1";
}
