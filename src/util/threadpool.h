

#ifndef MY_THREADPOOL_H
#define MY_THREADPOOL_H

#include "tinycthread.h"
typedef void (*task_func)(void* arg);



typedef struct task_node {
    task_func func;
    void* arg;
    struct task_node* next;
} task_node_t;

typedef struct {
    thrd_t* threads;
    int thread_count;
    task_node_t* head;
    task_node_t* tail;
    mtx_t lock;
    cnd_t notify;
    int stop;
} thread_pool_t;


extern int thread_worker(void* arg);
extern thread_pool_t* pool_create(int num_threads);
extern void pool_enqueue(thread_pool_t* pool, task_func func, void* arg);
extern void pool_destroy(thread_pool_t* pool);
#endif