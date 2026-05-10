
#ifndef lint
static const char RCSid[] = "$Id$";
#endif

/*
* matrix routines (avx2).
*
* Yongqing
*/

#if defined(_WIN32) || defined(_WIN64)
#include <psapi.h>
#elif defined(__APPLE__)
#define NULL ((void*)0)
#include <mach/mach.h>
#include <mach/mach_host.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/sysctl.h>
#else
#include <sys/mman.h>
#include <sys/stat.h>
/* #define NULL ((void*)0) */
#endif
#include "platform.h"
#include "resolu.h"
#include "color.h"
#include "cmmap.h"
#include "standard.h"
#include "cmatrix.h"
#include "cmatrixalign.h"
#include "threadpool.h"
#include <stdint.h>

typedef struct {
	int file_index;
	MappedFile* sequence;
	int w;
	int h;
	int64_t** index_table;
	int* processed_count;
	mtx_t* count_mutex;
} BuildIndexTaskData;

typedef struct {
	BuildIndexTaskData* task_data;
	int* tasks_completed;
	mtx_t* completion_mutex;
	cnd_t* completion_cnd;
	int num_total;
} IndexWrapperData;



NativeHandle read_single_file(const char* path) {
#if defined(_WIN32) || defined(_WIN64)
	NativeHandle handle = CreateFileA(
		path, GENERIC_READ, FILE_SHARE_READ, NULL,
		OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL
	);
	if (handle == NATIVE_INVALID_HANDLE) {
		fprintf(stderr, "Error：Cannot open File: %s\n", path);
		return NULL;
	}

#else

	NativeHandle handle = open(path, O_RDONLY);
	if (handle == -1) return -1;
#endif
	return handle;
}


NativeHandle* read_multihdrs(char* fspec, int num_images)
{

	NativeHandle* sequence = malloc(sizeof(NativeHandle) * num_images);

	for (int i = 0; i < num_images; i++) {
		char path[256];
		sprintf(path, fspec, i);
		sequence[i] = read_single_file(path);
	}

	fprintf(stderr, "Mapping the HDR images finished\n");
	return sequence;
}

void close_hdrfiles(NativeHandle* handles, int count) {
	if (handles == NULL) return;

	for (int i = 0; i < count; i++) {
		if (handles[i] != NATIVE_INVALID_HANDLE) {
#ifdef _WIN32
			CloseHandle(handles[i]);
#else
			close(handles[i]);
#endif
			handles[i] = NATIVE_INVALID_HANDLE;
		}
	}
}


void init_exponent_table() {
	for (int i = 0; i < 256; i++) {
		if (i == 0) exponent_table[i] = 0.0f;
		else exponent_table[i] = (float)ldexp(1.0, i - (128 + 8));
	}
}




void build_index_task(void* arg) {

	BuildIndexTaskData* t_data = (BuildIndexTaskData*)arg;
	int c = t_data->file_index;
	int w = t_data->w;
	int h = t_data->h;
	unsigned char* dummy_mem;
	unsigned char* dummy_channels[4];
	unsigned char* ptr;
	unsigned char* end;
	unsigned char* base_ptr;
	int y;
	int current_count;


	dummy_mem = (unsigned char*)malloc(4 * w);
	dummy_channels[0] = dummy_mem;
	dummy_channels[1] = dummy_mem + w;
	dummy_channels[2] = dummy_mem + 2 * w;
	dummy_channels[3] = dummy_mem + 3 * w;

	base_ptr = (unsigned char*)t_data->sequence[c].data;
	ptr = base_ptr;
	end = ptr + t_data->sequence[c].size;


	int found_res = 0;
	for (unsigned char* p = base_ptr; p < end - 4; p++) {
		if (p[0] == '\n' && (p[1] == '-' || p[1] == '+') && p[2] == 'Y' && p[3] == ' ') {
			ptr = p + 1; 
			found_res = 1;
			break;
		}
	}

	if (found_res) {
		
		while (ptr < end && *ptr != '\n') {
			ptr++;
		}
		
		if (ptr < end && *ptr == '\n') {
			ptr++; 
		}
	} else {
		fprintf(stderr,"Handle error: Resolution string not found\n");
		
	}


	for (y = 0; y < h; y++) {
		t_data->index_table[c][y] = (int64_t)(ptr - base_ptr);
		decode_hdr_line_planar(&ptr, dummy_channels, w);
	}

	t_data->index_table[c][y] = t_data->sequence[c].size;
	//int64_t filesize = get_file_size_by_handle(handle);

#ifndef _WIN32
	madvise(t_data->sequence[c].data, t_data->sequence[c].size, MADV_DONTNEED);
#else
	current_count = 0;
	mtx_lock(t_data->count_mutex);
	(*t_data->processed_count)++;
	current_count = *t_data->processed_count;
	mtx_unlock(t_data->count_mutex);

	if (current_count % 50 == 0) {
		mtx_lock(t_data->count_mutex);
		EmptyWorkingSet(GetCurrentProcess());
		mtx_unlock(t_data->count_mutex);
	}
#endif

	free(dummy_mem);
}

void _index_worker(void* arg) {
	IndexWrapperData* wrap = (IndexWrapperData*)arg;
	build_index_task((void*)wrap->task_data);

	mtx_lock(wrap->completion_mutex);
	(*wrap->tasks_completed)++;
	if (*wrap->tasks_completed == wrap->num_total) {
		cnd_signal(wrap->completion_cnd);
	}
	mtx_unlock(wrap->completion_mutex);
}

int64_t** build_scanline_index(MappedFile* sequence, int num_hdrs, int w, int h, thread_pool_t* pool) {

	int64_t** index_table;
	int processed_count = 0;
	mtx_t count_mutex;
	BuildIndexTaskData* tasks;
	IndexWrapperData* wrappers;
	int tasks_completed = 0;
	mtx_t completion_mutex;
	cnd_t completion_cnd;
	int i, c;

	mtx_init(&count_mutex, mtx_plain);
	mtx_init(&completion_mutex, mtx_plain);
	cnd_init(&completion_cnd);

	index_table = (int64_t**)malloc(num_hdrs * sizeof(int64_t*));
	for (i = 0; i < num_hdrs; i++) {
		index_table[i] = (int64_t*)malloc((h + 1) * sizeof(int64_t));
	}

	fprintf(stderr, "Building global scanline index table with ThreadPool...\n");

	tasks = (BuildIndexTaskData*)malloc(num_hdrs * sizeof(BuildIndexTaskData));
	wrappers = (IndexWrapperData*)malloc(num_hdrs * sizeof(IndexWrapperData));

	for (c = 0; c < num_hdrs; c++) {
		tasks[c].file_index = c;
		tasks[c].sequence = sequence;
		tasks[c].w = w;
		tasks[c].h = h;
		tasks[c].index_table = index_table;
		tasks[c].processed_count = &processed_count;
		tasks[c].count_mutex = &count_mutex;

		wrappers[c].task_data = &tasks[c];
		wrappers[c].tasks_completed = &tasks_completed;
		wrappers[c].completion_mutex = &completion_mutex;
		wrappers[c].completion_cnd = &completion_cnd;
		wrappers[c].num_total = num_hdrs;

		pool_enqueue(pool, _index_worker, &wrappers[c]);
	}


	mtx_lock(&completion_mutex);
	while (tasks_completed < num_hdrs) {
		cnd_wait(&completion_cnd, &completion_mutex);
	}
	mtx_unlock(&completion_mutex);


	free(tasks);
	free(wrappers);
	mtx_destroy(&count_mutex);
	mtx_destroy(&completion_mutex);
	cnd_destroy(&completion_cnd);

	fprintf(stderr, "Index table built successfully.\n");
	return index_table;
}



MappedFile map_single_file(const char* path) {
	MappedFile mf = { NULL, 0 };

#if defined(_WIN32) || defined(_WIN64)
	HANDLE hFile = CreateFileA(path, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_FLAG_SEQUENTIAL_SCAN, NULL);
	if (hFile == INVALID_HANDLE_VALUE) return mf;

	LARGE_INTEGER liSize;
	GetFileSizeEx(hFile, &liSize);
	mf.size = (size_t)liSize.QuadPart;

	HANDLE hMap = CreateFileMapping(hFile, NULL, PAGE_READONLY, 0, 0, NULL);
	if (hMap == NULL) { CloseHandle(hFile); return mf; }

	mf.data = MapViewOfFile(hMap, FILE_MAP_READ, 0, 0, 0);

	CloseHandle(hMap);
	CloseHandle(hFile);
#else
	int fd = open(path, O_RDONLY);
	if (fd == -1) return mf;

	struct stat st;
	fstat(fd, &st);
	mf.size = st.st_size;

	mf.data = mmap(NULL, mf.size, PROT_READ, MAP_PRIVATE, fd, 0);
	if (mf.data == MAP_FAILED) { mf.data = NULL; close(fd); return mf; }
	posix_madvise(mf.data, mf.size, POSIX_MADV_SEQUENTIAL);
	close(fd);
#endif

	return mf;
}


void unmap_file(MappedFile* mf) {
	if (!mf->data) return;
#ifdef _WIN32
	UnmapViewOfFile(mf->data);
//#else
//	munmap(mf->data, mf->size);
#endif
}





MappedFile* mapmultiHDRs(char* fspec, int* outw, int* outh, int* num_images)
{
	char filepath[MAX_PATH];
	snprintf(filepath, sizeof(filepath), fspec, 0);
	FILE* fp = fopen(filepath, "rb");

	if (fp == NULL) {
		char err[1024];
		sprintf(err, "Cannot find a file matches %s\n", fspec);
		error(SYSTEM, err);
	}

	int width, height;
	if (getheader(fp, NULL, NULL) < 0 || !fscnresolu(&width, &height, fp)) {
		fprintf(stderr, "Fatal: Invalid Radiance HDR header in %s\n", filepath);
		fclose(fp);
		return NULL;
	}
	fclose(fp);

	int total_images = 0;

	*outw = width;
	*outh = height;

	while (1) {
		snprintf(filepath, sizeof(filepath), fspec, total_images);
		FILE* tf = fopen(filepath, "rb");
		if (!tf) break;
		fclose(tf);
		total_images++;
	}

	if (total_images == 0) {
		fprintf(stderr, "Fatal: No files found matching pattern %s\n", fspec);
		return NULL;
	}

	*num_images = total_images;

	MappedFile* sequence = malloc(sizeof(MappedFile) * total_images);

	for (int i = 0; i < total_images; i++) {
		char path[256];
		sprintf(path, fspec, i);
		sequence[i] = map_single_file(path);
	}

	fprintf(stderr, "Mapping the HDR images finished\n");
	return sequence;
}

void unmapmultiHDRs(MappedFile* sequence, int num_images)
{
	for (int i = 0; i < num_images; i++) {
		unmap_file(&sequence[i]);
	}
	free(sequence);
}


int decode_hdr_line_planar_f(NativeHandle handle, int64_t** index_table, int c, int absolute_y, unsigned char* channels[4], int width) {
	int64_t currlineoffset = index_table[c][absolute_y];
	int64_t nextlineOffset = index_table[c][absolute_y + 1];
	int64_t size = nextlineOffset - currlineoffset;
	char* buffer = (char*)malloc((size_t)size);

#if defined(_WIN32) || defined(_WIN64)
	OVERLAPPED ol = {0};
	ol.Offset = (DWORD)(currlineoffset & 0xFFFFFFFF);
	ol.OffsetHigh = (DWORD)((currlineoffset >> 32) & 0xFFFFFFFF);
	DWORD bytesRead = 0;

	if (!ReadFile(handle, buffer, (DWORD)size, &bytesRead, &ol)) {
		fprintf(stderr, "Failed to read line at offset %lld\n", currlineoffset);
	
		return -1;
	}
	size = bytesRead; 
#else
	ssize_t bytesRead = pread(handle, buffer, size, currlineoffset);
	if (bytesRead <= 0) {
		fprintf(stderr, "Failed to read line at offset %lld\n", (long long)currlineoffset);
		
		return -1;
	}
	size = bytesRead;
#endif
	char* p = buffer;
	if (p[0] == 2 && p[1] == 2 && !(p[2] & 128)) {
		p += 4;
		for (int i = 0; i < 4; i++) {
			int j = 0;
			while (j < width) {
				unsigned char code = *p++;
				if (code > 128) {
					code &= 127;
					memset(&channels[i][j], *p++, code);
					j += code;
				}
				else {
					memcpy(&channels[i][j], p, code);
					p += code;
					j += code;
				}
			}
		}
	}
	else {
		for (int j = 0; j < width; j++) {
			channels[0][j] = *p++; // R
			channels[1][j] = *p++; // G
			channels[2][j] = *p++; // B
			channels[3][j] = *p++; // E
		}
	}
	free(buffer);
	return 0;
}



int decode_hdr_line_planar(unsigned char** src, unsigned char* channels[4], int width) {
	unsigned char* p = *src;
	if (p[0] == 2 && p[1] == 2 && !(p[2] & 128)) {
		p += 4;
		for (int i = 0; i < 4; i++) {
			int j = 0;
			while (j < width) {
				unsigned char code = *p++;
				if (code > 128) {
					code &= 127;
					memset(&channels[i][j], *p++, code);
					j += code;
				}
				else {
					memcpy(&channels[i][j], p, code);
					p += code;
					j += code;
				}
			}
		}
	}
	else {
		for (int j = 0; j < width; j++) {
			channels[0][j] = *p++; // R
			channels[1][j] = *p++; // G
			channels[2][j] = *p++; // B
			channels[3][j] = *p++; // E
		}
	}

	*src = p;
	return 0;
}



CMATRIX* obtain_scanline_matrix(NativeHandle* handles, int64_t** index_table, int start_y, int num_hdrs, int w, int num_lines) {
	CMATRIX* hdr_matrix = cm_alloc_u(w * num_lines, num_hdrs);
	unsigned char* planar_mem = (unsigned char*)malloc(4 * w);
	unsigned char* channels[4] = { planar_mem, planar_mem + w, planar_mem + 2 * w, planar_mem + 3 * w };
	COLORV* line_buffer = (COLORV*)malloc(w * num_hdrs * 3 * sizeof(COLORV));
	
	for (int y = 0; y < num_lines; y++) {
		memset(line_buffer, 0, w * num_hdrs * 3 * sizeof(COLORV));
		int absolute_y = start_y + y;
		for (int c = 0; c < num_hdrs; c++) {

			if (decode_hdr_line_planar_f(handles[c], index_table, c, absolute_y, channels, w) == 0) {
				int x = 0;
				for (; x < w; x++) {
					unsigned char e = channels[3][x];
					if (e > 0) {
						float f = exponent_table[e];
						COLORV* p = &line_buffer[x * num_hdrs * 3 + c * 3];
						p[0] = (COLORV)(channels[0][x] * f);
						p[1] = (COLORV)(channels[1][x] * f);
						p[2] = (COLORV)(channels[2][x] * f);
					}
				}
			}
		}

		memcpy(&hdr_matrix->cmem[y * w * num_hdrs * 3], line_buffer, w * num_hdrs * 3 * sizeof(COLORV));
	}

	free(line_buffer);
	free(planar_mem);
	return hdr_matrix;
}


unsigned long long get_available_memory() {
#if defined(_WIN32) || defined(_WIN64)
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    if (GlobalMemoryStatusEx(&status)) {
        return status.ullAvailPhys;
    }
    return 0;

#elif defined(__APPLE__)
    vm_size_t page_size;
    mach_port_t mach_port = mach_host_self();
    host_page_size(mach_port, &page_size);

    vm_statistics64_data_t vm_stat;
    mach_msg_type_number_t count = HOST_VM_INFO64_COUNT;

    if (host_statistics64(mach_port, HOST_VM_INFO64, (host_info64_t)&vm_stat, &count) == KERN_SUCCESS) {
        
        unsigned long long available_pages = vm_stat.free_count + vm_stat.inactive_count;
        return available_pages * (unsigned long long)page_size;
    }
    return 0;

#else
    
    long pages = sysconf(_SC_AVPHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);

    if (pages > 0 && page_size > 0) {
        return (unsigned long long)pages * (unsigned long long)page_size;
    }
    return 0;
#endif
}
