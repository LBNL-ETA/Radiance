#ifndef lint
static const char RCSid[] = "$Id$";
#endif
#include "tinycthread.h"
#include <stdlib.h>
#include "threadpool.h"

int thread_worker(void* arg) {
    thread_pool_t* pool = (thread_pool_t*)arg;
    while (1) {
        mtx_lock(&pool->lock);
        
        while (pool->head == NULL && !pool->stop) {
            cnd_wait(&pool->notify, &pool->lock);
        }
        
        
        if (pool->stop && pool->head == NULL) {
            mtx_unlock(&pool->lock);
            break;
        }

        
        task_node_t* task = pool->head;
        pool->head = task->next;
        if (pool->head == NULL) pool->tail = NULL;
        mtx_unlock(&pool->lock);

        
        task->func(task->arg);
        free(task);
    }
    return 0;
}


thread_pool_t* pool_create(int num_threads) {
    thread_pool_t* pool = (thread_pool_t*)malloc(sizeof(thread_pool_t));
    pool->thread_count = num_threads;
    pool->head = NULL;
    pool->tail = NULL;
    pool->stop = 0;
    mtx_init(&pool->lock, mtx_plain);
    cnd_init(&pool->notify);
    pool->threads = (thrd_t*)malloc(num_threads * sizeof(thrd_t));

    for (int i = 0; i < num_threads; i++) {
        thrd_create(&pool->threads[i], thread_worker, pool);
    }
    return pool;
}

void pool_enqueue(thread_pool_t* pool, task_func func, void* arg) {
    task_node_t* new_task = (task_node_t*)malloc(sizeof(task_node_t));
    new_task->func = func;
    new_task->arg = arg;
    new_task->next = NULL;

    mtx_lock(&pool->lock);
    if (pool->tail == NULL) {
        pool->head = new_task;
        pool->tail = new_task;
    } else {
        pool->tail->next = new_task;
        pool->tail = new_task;
    }
    cnd_signal(&pool->notify);
    mtx_unlock(&pool->lock);
}

void pool_destroy(thread_pool_t* pool) {
    mtx_lock(&pool->lock);
    pool->stop = 1;
    cnd_broadcast(&pool->notify);
    mtx_unlock(&pool->lock);

    for (int i = 0; i < pool->thread_count; i++) {
        thrd_join(pool->threads[i], NULL);
    }
    free(pool->threads);
    mtx_destroy(&pool->lock);
    cnd_destroy(&pool->notify);
    free(pool);
}
