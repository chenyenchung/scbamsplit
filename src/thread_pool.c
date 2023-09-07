//
// Created by Yen-Chung Chen on 9/6/23.
//

#include "thread_pool.h"
#include <string.h>
#include <stdlib.h>
#include <pthread.h>

static tpool_work_t *tpool_work_create(thread_func_t func, void *arg, size_t arg_size) {
    tpool_work_t *work;
    if (func == NULL) return NULL;

    work = calloc(1, sizeof(*work));
    work->func = func;
    work->arg = calloc(1, arg_size);
    memcpy(work->arg, arg, arg_size);
    work->next = NULL;
    return work;
}

static void tpool_work_destroy(tpool_work_t *work) {
    if (work == NULL) return;
    free(work->arg);
    free(work);
}

static tpool_work_t *tpool_work_get(tpool_t *tm) {
    tpool_work_t *work;
    if (tm == NULL) return NULL;

    work = tm->work_first;
    if (work == NULL) return NULL;

    if (work->next == NULL) {
        tm->work_first = NULL;
        tm->work_last = NULL;
    } else {
        tm->work_first = work->next;
    }
    return work;
}

static void *tpool_worker(void *arg) {
    tpool_t *tm = arg;
    tpool_work_t *work;

    while (1) {
        pthread_mutex_lock(&(tm->work_mutex));

        while (tm->work_first == NULL && !tm->stop) {
            pthread_cond_wait(&(tm->work_avail_cond), &(tm->work_mutex));
        }

        if (tm->stop) {
            break;
        }

        work = tpool_work_get(tm);
        tm->working_count++;
        tm->queue_length--;
        pthread_cond_broadcast(&(tm->retrieve_cond));
        pthread_mutex_unlock(&(tm->work_mutex));

        if (work != NULL) {
            work->func(work->arg);
            tpool_work_destroy(work);
        }

        pthread_mutex_lock(&(tm->work_mutex));
        tm->working_count--;
        if (!tm->stop && tm->working_count == 0 && tm->work_first == NULL) {
            pthread_cond_signal(&(tm->working_cond));
        }
        pthread_mutex_unlock(&(tm->work_mutex));
    }

    tm->thread_count--;
    pthread_cond_signal(&(tm->working_cond));
    pthread_mutex_unlock(&(tm->work_mutex));
    return NULL;
}

tpool_t *tpool_create(size_t num, size_t qmax) {
    tpool_t *tm;
    pthread_t thread;
    size_t i;

    if (num == 0) {
        num = 2;
    }

    tm = calloc(1, sizeof(*tm));
    tm->thread_count = num;
    tm->queue_length = 0;
    tm->queue_max = qmax;

    pthread_mutex_init(&(tm->work_mutex), NULL);
    pthread_cond_init(&(tm->work_avail_cond), NULL);
    pthread_cond_init(&(tm->working_cond), NULL);
    pthread_cond_init(&(tm->retrieve_cond), NULL);
    pthread_cond_init(&(tm->chunkq_return_cond), NULL);

    tm->work_first = NULL;
    tm->work_last = NULL;

    for (i = 0; i < num; i++) {
        pthread_create(&thread, NULL, tpool_worker, tm);
        pthread_detach(thread);
    }

    return tm;
}

void tpool_destroy(tpool_t *tm) {
    tpool_work_t *work;
    tpool_work_t *work2;

    if (tm == NULL) {
        return;
    }

    pthread_mutex_lock(&(tm->work_mutex));
    work = tm->work_first;

    while (work != NULL) {
        work2 = work->next;
        tpool_work_destroy(work);
        work = work2;
    }

    tm->stop = true;
    pthread_cond_broadcast(&(tm->work_avail_cond));
    pthread_mutex_unlock(&(tm->work_mutex));

    tpool_wait(tm);

    pthread_mutex_destroy(&(tm->work_mutex));
    pthread_cond_destroy(&(tm->work_avail_cond));
    pthread_cond_destroy(&(tm->working_cond));

    free(tm);
}

bool tpool_add_work(tpool_t *tm, thread_func_t func, void *arg, size_t arg_size) {
    tpool_work_t *work;

    if (tm == NULL) {
        return false;
    }

    work = tpool_work_create(func, arg, arg_size);
    if (work == NULL) {
        return false;
    }

    pthread_mutex_lock(&(tm->work_mutex));

    while (tm->queue_length ==  tm->queue_max) {
        pthread_cond_wait(&(tm->retrieve_cond), &(tm->work_mutex));
    }

    if (tm->work_first == NULL) {
        tm->work_first = work;
        tm->work_last = tm->work_first;
    } else {
        tm->work_last->next = work;
        tm->work_last = work;
    }

    tm->queue_length++;
    pthread_cond_broadcast(&(tm->work_avail_cond));
    pthread_mutex_unlock(&(tm->work_mutex));

    return true;
}

void tpool_wait(tpool_t *tm) {
    if (tm == NULL) return;

    pthread_mutex_lock(&(tm->work_mutex));
    while (1) {
        if ((!tm->stop && tm->queue_length != 0) || (!tm->stop && tm->working_count != 0) || (tm->stop && tm->thread_count != 0)) {
            // Must also check if the queue_length is 0 to prevent when the work is assigned
            // but tpool_wait acquire the lock and see no working threads.
            pthread_cond_wait(&(tm->working_cond), &(tm->work_mutex));
        } else {
            break;
        }
    }

    pthread_mutex_unlock(&(tm->work_mutex));
}



chunkq_t *chunkq_create () {
   chunkq_t *cq;
   cq = calloc(1, sizeof(chunkq_t));

   cq->qlength = 0;
   cq->first_chunk = NULL;
   cq->last_chunk = NULL;
   pthread_mutex_init(&(cq->ext_lock), NULL);
   pthread_cond_init(&(cq->ext_return_cond), NULL);
   return cq;
}

void chunkq_destroy(chunkq_t *cq) {
    pthread_mutex_lock(&(cq->ext_lock));
    while(NULL != cq->first_chunk) {
        ichunk_t *ic;
        ic = cq->first_chunk;
        cq->first_chunk = ic->next_chunk;
        chunk_destroy(ic->chunk, ic->chunk_size);
        free(ic);
    }
    pthread_mutex_unlock(&(cq->ext_lock));
    free(cq);
}

ichunk_t *create_chunk(chunkq_t *cq, int64_t chunk_size) {
    ichunk_t *ic;
    ic = calloc(1, sizeof(ichunk_t));

    sam_read_t **chunk = chunk_init(chunk_size);

    ic->chunk = chunk;
    ic->next_chunk = NULL;
    ic->chunk_size = chunk_size;
    ic->read_kept = 0;
    return ic;
}

int8_t chunkq_ele_init(chunkq_t *cq, int64_t chunk_size) {
    if (NULL == cq || chunk_size < 1) return 1;

    pthread_mutex_lock(&(cq->ext_lock));
    ichunk_t *new_chunk;
    new_chunk = create_chunk(cq, chunk_size);
   if (cq->first_chunk == NULL && cq->qlength == 0) {
       cq->first_chunk = new_chunk;
       cq->last_chunk = new_chunk;
   } else {
       cq->last_chunk->next_chunk = new_chunk;
       cq->last_chunk = new_chunk;
   }
    cq->qlength++;
    pthread_mutex_unlock(&(cq->ext_lock));
    return 0;
}
int8_t chunkq_add(chunkq_t *cq, ichunk_t *ic) {
    if (NULL == cq | NULL == ic) return 1;
    pthread_mutex_lock(&(cq->ext_lock));
    ic->next_chunk = NULL;
    if (cq->first_chunk == NULL && cq->qlength == 0) {
        cq->first_chunk = ic;
        cq->last_chunk = ic;
    } else {
        cq->last_chunk->next_chunk = ic;
        cq->last_chunk = ic;
    }
    cq->qlength++;
    pthread_cond_broadcast(&(cq->ext_return_cond));
    pthread_mutex_unlock(&(cq->ext_lock));
    return 0;
}

ichunk_t *chunkq_get(chunkq_t *cq) {
    if (NULL == cq) return NULL;
    pthread_mutex_lock(&(cq->ext_lock));
    while (cq->first_chunk == NULL) {
        pthread_cond_wait(&(cq->ext_return_cond), &(cq->ext_lock));
    }
    ichunk_t *ic = cq->first_chunk;
    if (ic->next_chunk == NULL) {
        cq->last_chunk = NULL;
    }
    cq->first_chunk = ic->next_chunk;
    cq->qlength--;
    pthread_mutex_unlock(&(cq->ext_lock));
    return ic;
}
