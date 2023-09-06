//
// Created by Yen-Chung Chen on 9/6/23.
//

#ifndef SCBAMSPLIT_THREAD_POOL_H
#define SCBAMSPLIT_THREAD_POOL_H
#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>
#include <pthread.h>

// Adapted from https://nachtimwald.com/2019/04/12/thread-pool-in-c/amp/
typedef void (*thread_func_t)(void *arg);
struct tpool_work {
    thread_func_t func;
    void *arg;
    struct tpool_work *next;
};
typedef struct tpool_work tpool_work_t;
struct tpool {
    tpool_work_t *work_first;
    tpool_work_t *work_last;
    size_t queue_length;
    size_t queue_max;
    pthread_mutex_t work_mutex;
    pthread_cond_t work_avail_cond;
    pthread_cond_t working_cond;
    pthread_cond_t retrieve_cond;
    pthread_cond_t chunkq_return_cond;
    size_t working_count;
    size_t thread_count;
    bool stop;
};
typedef struct tpool tpool_t;

tpool_t *tpool_create(size_t num, size_t qmax);
void tpool_destroy(tpool_t *tm);
bool tpool_add_work(tpool_t *tm, thread_func_t func, void *arg, size_t arg_size);
void tpool_wait(tpool_t *tm);

struct chunkq;
typedef struct chunkq chunkq_t;
struct idv_chunk;
typedef struct idv_chunk ichunk_t;

struct chunkq {
    size_t qlength;
    size_t qmax;
    ichunk_t *first_chunk;
    ichunk_t *last_chunk;
    pthread_mutex_t ext_lock;
    pthread_cond_t ext_return_cond;
};

#include "utils.h"
struct idv_chunk {
    sam_read_t **chunk;
    int64_t chunk_size;
    int64_t read_kept;
    ichunk_t *next_chunk;
};

chunkq_t *chunkq_create();
void chunkq_destroy(chunkq_t *cq);
int8_t chunkq_ele_init(chunkq_t *cq, int64_t chunk_size);
int8_t chunkq_add(chunkq_t *cq, ichunk_t *ic);
ichunk_t *chunkq_get(chunkq_t *cq);

#endif //SCBAMSPLIT_THREAD_POOL_H
