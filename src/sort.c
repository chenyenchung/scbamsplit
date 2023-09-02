//
// Created by Yen-Chung Chen on 8/23/23.
//
#include "sort.h"
#include <string.h>
#include <stdbool.h>
#include "htslib/sam.h" /* include stdint.h with it */
#include "htslib/kstring.h"
#include "utils.h"

#define FREE_FILL_CHUNK bam_destroy1(temp_read); free(PR); free(MAPQ); free(UB); free(CB); free(RN)

int8_t get_tag(bam1_t *read, char* tag, char* tag_ptr) {
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

int8_t get_CB(bam1_t *read, char* tag, char* tag_ptr) {
    return get_tag(read, tag, tag_ptr);
}

int8_t get_UB(bam1_t *read, char* tag, char* tag_ptr) {
    return get_tag(read, tag, tag_ptr);
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

    if ((NULL == (void *) not_secondary) || (NULL == (void *) mapped)) {
        log_msg( "Fail to read bitwise flag", ERROR);
        return 1;
    }

    // This is used for sorting, so we want primary to be coded 1 while others being 2
    sprintf(flag_ptr, "%d", 2 - not_secondary && mapped);
    return 0;
}

int8_t get_MAPQ(bam1_t *read, char* val_ptr) {
    uint8_t mapq = (read)->core.qual;
    if (NULL == mapq) {
        return 1;
    }
    sprintf(val_ptr, "%03d", mapq);
    return 0;
}

int64_t fill_chunk(
        samFile *fp, sam_hdr_t *header, sam_read_t **read_array,
        int64_t chunk_size
        ) {
    /**
     * @abstract Fill read buffer to designated size and return the index of next read to read or -1
     * when fails.
     *
     * @fp A file descriptor for a SAM/BAM file
     * @header A header pointer from sam_hdr_read()
     * @read A pointer to a sam_read_t array to be filled
     * @chunk_size An integer to indicate how large the cache chunk to be filled
     * @returns The number of reads that have been allocated into the chunk on success; -1 on error;
     * -[OBSERVED_READ_NAME_SIZE] when read names are not sufficiently padded
     */

    // Note that read_array must be allocated OUTSIDE!
    bam1_t *temp_read = bam_init1();
    int64_t read_kept = 0;
    char *CB;
    char *UB;
    char *MAPQ;
    char *PR;
    char *RN; // Declare empty string of sufficient size
    CB = (char *) malloc(sizeof(char) * CB_LENGTH);
    UB = (char *) malloc(sizeof(char) * UB_LENGTH);
    MAPQ = (char *) malloc(sizeof(char) * 4); // Max value for MAPQ is 255 per SAM spec v1
    PR = (char *) malloc(sizeof(char) * 2); // The value for primary read is either 0 or 1
    RN = (char *) malloc(sizeof(char) * RN_SIZE);

    // Fill the chunk until specified size or running out of reads
    while (read_kept < chunk_size) {
        if (sam_read1(fp, header, temp_read) < 0) {
            // sam_read1 returns -1 when encountering an error
            // or EOF
            if (read_kept > 0) {
                FREE_FILL_CHUNK;
                return read_kept - 1;
            }
            FREE_FILL_CHUNK;
            return -1;
        }

        // TODO: need to be platform compatible.
        int8_t cb_stat = get_CB(temp_read, "CB", CB);
        int8_t ub_stat = get_UB(temp_read, "UB", UB);
        int8_t prim_stat = is_primary(temp_read, PR);
        int8_t mapq_stat = get_MAPQ(temp_read, MAPQ);

        // Skip reads that miss CB, UB, bitwise flag, or MAPQ
        if ((-1 == cb_stat) || (-1 == ub_stat) || (1 == prim_stat) || (1 == mapq_stat)) {
            // If a read has no CBC or UMI, just let it go.
            continue;
        } else if (1 == cb_stat || 1 == ub_stat) {
            // Error message is produced in get_tag()
            FREE_FILL_CHUNK;
            return -1;
        }

        bam1_t *read_copy_status = bam_copy1(read_array[read_kept]->read, temp_read);
        if (NULL == read_copy_status) {
            log_msg( "Fail to copy a BAM read. Could please an HTSlib issue?", ERROR);
            FREE_FILL_CHUNK;
            return -1;
        }

        // Prepare format string to pad the read name
        char rn_length[8];
        char fmt[16];

        // Get the read name to be padded to declared size - 1
        sprintf(rn_length, "%d", RN_SIZE - 1);

        // Prepare format string for sprintf padding
        strcpy(fmt, "%-");
        strcat(fmt, rn_length);
        strcat(fmt, "s");

        // Pad read name
        uint32_t obs_rn_size = strlen(bam_get_qname(temp_read));
        if (obs_rn_size > (RN_SIZE - 1)) {
            FREE_FILL_CHUNK;
            return -obs_rn_size;
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
                sizeof (char) * (strlen(read_array[read_kept]->key) + 1),
                (const uint8_t*) read_array[read_kept]->key
                );
        if (0 != app_stat) {
            log_msg("Fail to append sorting key to the read (%s)",
            WARNING, read_array[read_kept]->key);
        }

        read_kept++;
    }
    FREE_FILL_CHUNK;
    return read_kept;
}

int read_cmp(const void *a, const void *b) {
    const sam_read_t * reada = *(sam_read_t **) a;
    const sam_read_t * readb = *(sam_read_t **) b;

    return strcmp(reada->key, readb->key);
}

int sort_chunk(sam_read_t **reads, int64_t chunk_size) {
    qsort(reads, chunk_size, sizeof(sam_read_t*), read_cmp);
    return 0;
}

sam_read_t** chunk_init(uint32_t chunk_size) {
    sam_read_t **chunk;
    chunk = malloc(sizeof(sam_read_t*) * chunk_size);

    if (NULL == chunk) {
        log_msg("Fail to allocate memory for sorting", ERROR);
        return NULL;
    }

    for (uint32_t read_i = 0; read_i < chunk_size; read_i++) {
        chunk[read_i] = malloc(sizeof(sam_read_t));
        chunk[read_i]->read = bam_init1();
        chunk[read_i]->key = malloc(sizeof(char) * KEY_SIZE);
    }

    // A HTSlib-defined BAM read contains:
    // A bam1_t struct, which also points to a bam1_core_t

    return chunk;
}

void chunk_destroy(sam_read_t **read_array, uint32_t chunk_size) {
    for (uint32_t read_i = 0; read_i < chunk_size; read_i++) {
        bam_destroy1(read_array[read_i]->read);
        free(read_array[read_i]->key);
        free(read_array[read_i]);
    }
    free(read_array);
}

char* process_bam(samFile *fp, sam_hdr_t *header, sam_read_t **chunk, int64_t chunk_size, char *oprefix) {
    /**
     * @abstract Process all reads in an opened SAM/BAM file in chunks and save sorted reads in a temporary
     * directory.
     * @fp A pointer to a samFile
     * @header A pointer to a SAM file header
     * @chunk An array of sam_read_t's
     * @chunk_size The size of the array (a 64-bit integer)
     * @oprefix Output dir prefix
     * @returns A string: the path for the temporary directory containing sorted chunks if succeeded; "-1" if failed.
     */

    int64_t size_retrieved = chunk_size;
    int chunk_num = 0;
    // Create temporary file dir for sorted chunks
    char *tmpdir = create_tempdir(oprefix);


    while (size_retrieved == chunk_size) {
        chunk_num += 1;
        log_msg("Chunk #%lld filling", DEBUG, chunk_num);
        size_retrieved = fill_chunk(fp, header, chunk, chunk_size);
        if (size_retrieved < -1) {
            log_msg("Insufficient RN size (%d). Please increase RN size to at least %d", ERROR,
                    RN_SIZE, -size_retrieved + 1);
            return "1";
        } else if (size_retrieved == -1) {
            // Error message is generated in fill_chunk()
            return "1";
        }
        sort_chunk(chunk, size_retrieved);
        sam_hdr_change_HD(header, "SO", "unknown");

        char *tname = tname_init(tmpdir, "chunk", 5, chunk_num);

        htsFile* tfp = sam_open(tname, "wb");
        free(tname);

        int write_status = sam_hdr_write(tfp, header);
        if (write_status != 0) {
            log_msg( "Fail to write SAM header into temporary files", ERROR);
            return "1";
        }

        int write_to_bam = 0;
        for (int64_t i = 0; i < size_retrieved; i++) {
            write_to_bam = sam_write1(tfp, header, chunk[i]->read);
            if (write_to_bam == -1) {
                log_msg("Fail to write sorted reads into temporary files", ERROR);
                return "1";
            }
        }
        sam_close(tfp);

        // Exit early if in dev mode
//        if (dev && chunk_num > 2) break;
    }
    return tmpdir;
}
