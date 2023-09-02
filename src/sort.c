//
// Created by Yen-Chung Chen on 8/23/23.
//
#include "sort.h"
#include <string.h>
#include <stdbool.h>
#include "htslib/sam.h" /* include stdint.h with it */
#include "htslib/kstring.h"

const char* get_tag(bam1_t *read, char* tag) {
    kstring_t *tag_val;
    tag_val = malloc(sizeof(kstring_t));
    if (0 == bam_aux_get_str(read, tag, tag_val)) {
        return NULL;
    }

    char *tag_ptr;
    tag_ptr = malloc(sizeof(char) * (ks_len(tag_val) - 5));
    strcpy(tag_ptr, ks_str(tag_val) + 5);
    free(tag_val);
    return tag_ptr;
}

const char* get_CB(bam1_t *read, char* tag) {
    char* cb = get_tag(read, tag);
    if (NULL == cb) {return NULL;}
    return cb;
}

const char* get_UB(bam1_t *read, char* tag) {
    char* ub = get_tag(read, tag);
    if (NULL == ub) {return NULL;}
    return ub;
}



const char* is_primary(bam1_t *read) {
    bool not_secondary = !((read)->core.flag & 256); // Flag for secondary alignment is 256
    bool mapped = !((read)->core.flag & 4); // Flag for unmapped is 4
    char* primary;
    primary = malloc(sizeof(char) * 2);
    sprintf(primary, "%d", not_secondary && mapped);
    return primary;
}

const char* get_MAPQ(bam1_t *read) {
    uint8_t mapq = (read)->core.qual;
    char mapq_s[4] = "";
    char* mapq_spad;
    mapq_spad = malloc(sizeof(char) * 4);
    sprintf(mapq_s, "%d", mapq);
    sprintf(mapq_spad, "%-3s", mapq_s);
    return mapq_spad;
}

int64_t fill_chunk(
        samFile *fp, sam_hdr_t *header, sam_read read_array[],
        int64_t chunk_size
        ) {
    /**
     * @abstract Fill read buffer to designated size and return the index of next read to read or -1
     * when fails.
     *
     * @fp A file descriptor for a SAM/BAM file
     * @header A header pointer from sam_hdr_read()
     * @read A pointer to a sam_read array to be filled
     * @chunk_size An integer to indicate how large the cache chunk to be filled
     * @return The number of reads that have been allocated into the chunk on success; -1 on error or EOF;
     * -[OBSERVED_READ_NAME_SIZE] when read names are not sufficiently padded
     */

    // Note that read_array must be allocated OUTSIDE!
    bam1_t *temp_read = bam_init1();
    int64_t read_kept = 0;

    // Fill the chunk until specified size or running out of reads
    while (read_kept < chunk_size) {
        if (sam_read1(fp, header, temp_read) < 0) {
            // sam_read1 returns -1 when encountering an error
            // or EOF
            bam_destroy1(temp_read);
            if (read_kept > 0) {
                return read_kept - 1;
            }
            return -1;
        }

        // TODO: need to be platform compatible.
        const char *CB = get_CB(temp_read, "CB");
        const char *UB = get_UB(temp_read, "UB");

        if (CB == NULL | UB == NULL) {
            // If a read has no CBC or UMI, just let it go.
            continue;
        }

        bam_copy1(read_array[read_kept].read, temp_read);

        // Get padded read name
        char RN[RN_SIZE] = ""; // Declare empty string of sufficient size

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
            return -obs_rn_size;
        }
        sprintf(RN, fmt, bam_get_qname(temp_read));


        const char *PR = is_primary(temp_read);
        const char *MAPQ = get_MAPQ(temp_read);


        // Initialize the key string
        strcpy(read_array[read_kept].key, "");

        // Concatenate CBC, UMI, read name, primary or secondary, MAPQ
        // Lump reads that are supposed to be the same molecule
        strcat(read_array[read_kept].key, CB);
        strcat(read_array[read_kept].key, UB);
        // Show primary mapping first and sort by MAPQ
        // so later exporting mechanism can just export the
        // first and all secondary mappings (if exists) by
        // read name
        strcat(read_array[read_kept].key, PR);
        strcat(read_array[read_kept].key, MAPQ);
        strcat(read_array[read_kept].key, RN);


        free((char *) PR);
        free((char *) MAPQ);
        free((char *) UB);
        free((char *) CB);
        read_kept++;
    }
    bam_destroy1(temp_read);

    return read_kept;
}

int read_cmp(const void *a, const void *b) {
    const sam_read * reada = (sam_read *) a;
    const sam_read * readb = (sam_read *) b;

    return strcmp(reada->key, readb->key);
}

int sort_chunk(sam_read reads[], int64_t chunk_size) {
    qsort(reads, chunk_size, sizeof(sam_read), read_cmp);
    return 0;
}

void chunk_init(sam_read *read_array, uint32_t chunk_size) {
    for (uint32_t read_i = 0; read_i < chunk_size; read_i++) {
        read_array[read_i].read = bam_init1();
        read_array[read_i].key = malloc(sizeof(char) * KEY_SIZE);
    }
}

void chunk_destroy(sam_read *read_array, uint32_t chunk_size) {
    for (uint32_t read_i = 0; read_i < chunk_size; read_i++) {
        bam_destroy1(read_array[read_i].read);
        free(read_array[read_i].key);
    }
}
