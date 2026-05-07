/* RCSid $Id$ */
/*
 * Daylight coefficient matrix process routine declarations.
 *
 *	Yongqing
 */


#ifndef _CMMAP_H_ 
#define _CMMAP_H_

#include "cmatrix.h"
#include "threadpool.h"
#include <stdint.h>
#ifdef _WIN32
typedef HANDLE NativeHandle;
#define NATIVE_INVALID_HANDLE INVALID_HANDLE_VALUE
#else
#include <fcntl.h>
#include <unistd.h>
typedef int NativeHandle;
#define NATIVE_INVALID_HANDLE (-1)
#endif


#define MAX_PATH 1024
extern float exponent_table[256];
typedef struct {
	void* data;
	size_t size;
} MappedFile;



typedef unsigned char RGBE[4];

extern int decode_hdr_line_planar(unsigned char** src, unsigned char* channels[4], int width);

extern unsigned long long get_available_memory();

extern void init_exponent_table();

extern MappedFile map_single_file(const char* path);

extern MappedFile* mapmultiHDRs(char* fspec, int* outw, int* outh, int* num_images);

extern CMATRIX* obtain_scanline_matrix(NativeHandle* handles, int64_t** index_table, int start_y, int num_hdrs, int w, int num_lines);

extern void unmap_file(MappedFile* mf);

extern void unmapmultiHDRs(MappedFile* sequence, int num_images);

extern CMATRIX* obtain_scanline_matrix(NativeHandle* handles, int64_t** index_table, int start_y, int num_hdrs, int w, int num_lines);
int64_t** build_scanline_index(MappedFile* sequence, int num_hdrs, int w, int h, thread_pool_t* pool);
extern void _index_worker(void* arg);

extern NativeHandle* read_multihdrs(char* fspec, int num_images);
extern void close_hdrfiles(NativeHandle* handles, int count);

#endif