#ifndef lint
static const char RCSid[] = "$Id$";
#endif
/*
 * Compute time-step result using Daylight Coefficient method.
 *
 *	G. Ward
 *	Intel optimizations by Yongqing Zhao
 */
#ifdef INTEL_DCTOPT
#ifdef _WIN32
	#define access _access
	#define getcwd _getcwd
#else
	#include <sys/mman.h>
	#include <sys/resource.h>
#endif
#include <immintrin.h> // AVX2 instruction
#include <stdint.h>
#endif /* INTEL_DCTOPT */
#include <ctype.h>
#include "platform.h"
#include "standard.h"
#include "cmatrix.h"
#include "resolu.h"

#ifdef INTEL_DCTOPT

#define MAGIC 0x424C4F42
#define MAX_PATH 1024

float exponent_table[256];

void init_exponent_table() {
	for (int i = 0; i < 256; i++) {
		if (i == 0) exponent_table[i] = 0.0f;
		else exponent_table[i] = (float)ldexp(1.0, i - (128 + 8));
	}
}

unsigned long long get_safe_buffer_size() {
#ifdef _WIN32
	MEMORYSTATUSEX statex;
	statex.dwLength = sizeof(statex);
	if (GlobalMemoryStatusEx(&statex)) return (unsigned long long)(statex.ullAvailPhys * 0.7);
#else
	long pages = sysconf(_SC_PHYS_PAGES);
	long page_size = sysconf(_SC_PAGE_SIZE);
	if (pages > 0 && page_size > 0) return (unsigned long long)((pages * page_size) * 0.7);
#endif
	return 1024ULL * 1024ULL * 1024ULL;
}

#ifndef _WIN32
void set_max_files(int num_files) {
	struct rlimit rl;

	if (getrlimit(RLIMIT_NOFILE, &rl) == -1) {

		perror("getrlimit");
		return;
	}

	rl.rlim_cur = num_files;
	if (num_files > rl.rlim_max) {
		rl.rlim_max = num_files;
	}

	if (setrlimit(RLIMIT_NOFILE, &rl) == -1) {
		perror("setrlimit");
	}
}
#endif

void*
create_and_map_swap(const char* temp_file, size_t size) {
#ifdef _WIN32
	HANDLE hFile = CreateFileA(temp_file, GENERIC_READ | GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
	if (hFile == INVALID_HANDLE_VALUE) return NULL;

	LARGE_INTEGER li; li.QuadPart = size;
	SetFilePointerEx(hFile, li, NULL, FILE_BEGIN);
	if (!SetEndOfFile(hFile)) {
		if (GetLastError() == ERROR_DISK_FULL) {
			printf("Out of storage: Failed to allocate space for this file.\n");
			return NULL;
		}
	}

	HANDLE hMap = CreateFileMapping(hFile, NULL, PAGE_READWRITE, 0, 0, NULL);
	if (hMap == NULL) {
		printf("CreateFileMapping failed: %d\n", GetLastError());
		return NULL;
	}

	void* p = MapViewOfFile(hMap, FILE_MAP_ALL_ACCESS, 0, 0, 0);
	CloseHandle(hMap); CloseHandle(hFile);
	return p;
#else
	int fd = open(temp_file, O_RDWR | O_CREAT | O_TRUNC, 0666);
	if (fd == -1) return NULL;
	int res = posix_fallocate(fd, 0, size);
	if (res != 0) {
		error(SYSTEM, "Out of storage: Failed to allocate space for this file.");
		close(fd);
		return NULL;
	}
	void* p = mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
	close(fd);
	return (p == MAP_FAILED) ? NULL : p;
#endif
}




CMATRIX*
allocate_and_decompress_cmatrix_mmap(const void* src_data, int w, int h, int nimg, const char* swap_path) {
	size_t total_pixels = (size_t)w * h;
	size_t total_floats = total_pixels * nimg * 3;


	size_t total_mem = sizeof(CMATRIX) + (total_floats - 3) * sizeof(COLORV);

	CMATRIX* cm = (CMATRIX*)create_and_map_swap(swap_path, total_mem);
	if (!cm) return NULL;


	cm->nrows = total_pixels;
	cm->ncols = nimg;


	COLORV* dest_base = cm->cmem;

	const uint8_t* p_src = (const uint8_t*)src_data;
	size_t floats_processed = 0;

	printf("Decompressing data sequentially to virtual mapping area (Swap)... \n");


	while (floats_processed < total_floats) {
		float val = *(float*)p_src;
		p_src += sizeof(float);

		uint32_t count = *(uint32_t*)p_src;
		p_src += sizeof(uint32_t);


		for (uint32_t i = 0; i < count && floats_processed < total_floats; i++) {
			dest_base[floats_processed++] = (COLORV)val;
		}
	}
	return cm;
}

CMATRIX*
map_as_cmatrix(const char* filename, int* orig_w, int* orig_h) {
	void* p = NULL;
	size_t fileSize = 0;

#ifdef _WIN32
	HANDLE hFile = CreateFileA(filename, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
	if (hFile == INVALID_HANDLE_VALUE) return NULL;
	LARGE_INTEGER fs; GetFileSizeEx(hFile, &fs); fileSize = (size_t)fs.QuadPart;
	HANDLE hMapping = CreateFileMapping(hFile, NULL, PAGE_READONLY, 0, 0, NULL);
	p = MapViewOfFile(hMapping, FILE_MAP_READ, 0, 0, 0);
	CloseHandle(hMapping); CloseHandle(hFile);
#else
	int fd = open(filename, O_RDONLY);
	if (fd == -1) return NULL;
	struct stat st; fstat(fd, &st); fileSize = st.st_size;
	p = mmap(NULL, fileSize, PROT_READ, MAP_SHARED, fd, 0);
	close(fd);
#endif

	if (!p || p == (void*)-1) return NULL;
	int* pInt = (int*)p;

	if (pInt[0] != MAGIC) {

#ifdef _WIN32
		UnmapViewOfFile(p);
#else
		munmap(p, fileSize);
#endif
		return NULL;
	}

	int w = pInt[1], h = pInt[2], is_comp = pInt[3], nimg = pInt[7];
	if (orig_w) *orig_w = w; if (orig_h) *orig_h = h;

	CMATRIX* cm = NULL;
	if (!is_comp) {

		cm = allocate_and_decompress_cmatrix_mmap(pInt + 8, w, h, nimg, "cmatrix_swap.tmp");


	}
	else {

		size_t pixels_per_channel = (size_t)w * h * nimg;
		size_t total_mem = sizeof(CMATRIX) + (3 * pixels_per_channel - 3) * sizeof(COLORV);

		cm = (CMATRIX*)malloc(total_mem);
		if (cm) {
			cm->nrows = w * h;
			cm->ncols = nimg;

			COLORV* r_base = &cm->cmem[0];
			COLORV* g_base = r_base + pixels_per_channel;
			COLORV* b_base = g_base + pixels_per_channel;

			float* src_floats = (float*)(pInt + 8);

			for (size_t i = 0; i < pixels_per_channel; i++) {
				r_base[i] = (COLORV)src_floats[i * 3 + 0];
				g_base[i] = (COLORV)src_floats[i * 3 + 1];
				b_base[i] = (COLORV)src_floats[i * 3 + 2];
			}
		}
	}


#ifdef _WIN32
	UnmapViewOfFile(p);
#else
	munmap(p, fileSize);
#endif

	return cm;
}


void free_mmap_cmatrix(CMATRIX* cm, int w, int h, int nimg) {
	//printf("free_mmap_cmatrix\n");
	if (!cm) return;
	size_t pixels_per_channel = (size_t)w * h * nimg;
	size_t total_mem = sizeof(CMATRIX) + (3 * pixels_per_channel - 3) * sizeof(COLORV);

#ifdef _WIN32
	UnmapViewOfFile(cm);
	DeleteFileA("cmatrix_swap.tmp");
#else
	munmap(cm, total_mem);
	unlink("cmatrix_swap.tmp");
#endif
	//printf("End free_mmap_cmatrix\n");
}


void
fwrite_rle_float(float* data, size_t num_elements, FILE* fp) {
	if (num_elements == 0) return;

	for (size_t i = 0; i < num_elements; i++) {
		float current_val = data[i];
		uint32_t count = 1;


		while (i + 1 < num_elements && data[i + 1] == current_val) {
			count++;
			i++;

			if (count == 0xFFFFFFFF) break;
		}
		fwrite(&current_val, sizeof(float), 1, fp);
		fwrite(&count, sizeof(uint32_t), 1, fp);
	}
}



int hdr2bin(char* inputpattern, char* outfile) {
	char* out_name = NULL, * pattern = NULL;
#ifdef _WIN32
	_setmaxstdio(8000);
#else
	set_max_files(9999);
#endif

	out_name = outfile;
	pattern = inputpattern;

	init_exponent_table();


	///printf("Scanning files...\n");
	int nimg = 0;
	while (1) {
		char tmp[MAX_PATH];
		sprintf(tmp, pattern, nimg);
		FILE* tf = fopen(tmp, "rb");
		if (!tf) break;
		fclose(tf);
		nimg++;
	}

	if (nimg == 0) {
		char errtxt[MAX_PATH];
		//sprintf(errtxt, "Fatal: No files found matching pattern(%s).\n", pattern);
		fprintf(stderr, "Fatal: No files found matching pattern(%s)\n", pattern);
		//error(SYSTEM, errtxt);
		return 1;
	}


	FILE** fps = (FILE**)malloc(sizeof(FILE*) * nimg);
	int width, height;
	for (int i = 0; i < nimg; i++) {
		char fname[MAX_PATH];
		sprintf(fname, pattern, i);
		//printf("%s\n", fname);
		fps[i] = fopen(fname, "rb");
		if (i == 0) {
			if (getheader(fps[i], NULL, NULL) < 0 || !fscnresolu(&width, &height, fps[i])) {
				fprintf(stderr, "Fatal: Invalid Radiance HDR header.\n");
				return 1;
			}
		}
		else {
			int tw, th;
			getheader(fps[i], NULL, NULL);
			fscnresolu(&tw, &th, fps[i]);
		}
	}

	int total_pixels = width * height;
	unsigned long long buf_limit = get_safe_buffer_size();


	size_t bytes_per_pixel = 3 * sizeof(float) + sizeof(COLR);
	int rows_per_chunk = (int)(buf_limit / (width * nimg * bytes_per_pixel));
	if (rows_per_chunk < 1) rows_per_chunk = 1;
	if (rows_per_chunk > height) rows_per_chunk = height;

	printf("Resolution: %dx%d, Images: %d\n", width, height, nimg);
	printf("Memory Buffer: %llu MB, Rows per chunk: %d\n", buf_limit / (1024 * 1024), rows_per_chunk);


	FILE* fout = fopen(out_name, "wb");
	int full_header[8] = { MAGIC, width, height, 0, 0, 0, total_pixels, nimg };
	fwrite(full_header, sizeof(int), 8, fout);


	float* chunk_buf = (float*)malloc((size_t)rows_per_chunk * width * nimg * 3 * sizeof(float));


	COLR* raw_chunk_colrs = (COLR*)malloc((size_t)rows_per_chunk * width * nimg * sizeof(COLR));

	for (int y_start = 0; y_start < height; y_start += rows_per_chunk) {
		int y_end = (y_start + rows_per_chunk > height) ? height : y_start + rows_per_chunk;
		int current_chunk_rows = y_end - y_start;
		int total_pixels_in_chunk = current_chunk_rows * width;

		printf("Processing lines %d to %d...\n", y_start, y_end);


		for (int i = 0; i < nimg; i++) {
			for (int y = 0; y < current_chunk_rows; y++) {

				size_t scanline_offset = ((size_t)i * current_chunk_rows + y) * width;

				if (freadcolrs(&raw_chunk_colrs[scanline_offset], width, fps[i]) < 0) {
					memset(&raw_chunk_colrs[scanline_offset], 0, sizeof(COLR) * width);
				}
			}
		}

		int p;
#pragma omp parallel for
		for (p = 0; p < total_pixels_in_chunk; p++) {
			for (int i = 0; i < nimg; i++) {

				size_t src_idx = (size_t)i * total_pixels_in_chunk + p;


				size_t dst_idx = ((size_t)p * nimg + i) * 3;

				unsigned char r = raw_chunk_colrs[src_idx][0];
				unsigned char g = raw_chunk_colrs[src_idx][1];
				unsigned char b = raw_chunk_colrs[src_idx][2];
				unsigned char e = raw_chunk_colrs[src_idx][3];

				float f = exponent_table[e];
				chunk_buf[dst_idx + 0] = (r + 0.5f) * f;
				chunk_buf[dst_idx + 1] = (g + 0.5f) * f;
				chunk_buf[dst_idx + 2] = (b + 0.5f) * f;
			}
		}
		size_t elements_to_write = (size_t)current_chunk_rows * width * nimg * 3;
		fwrite_rle_float(chunk_buf, elements_to_write, fout);
		//fwrite(chunk_buf, sizeof(float), (size_t)total_pixels_in_chunk * nimg * 3, fout);
	}

	printf("\nSuccess: Saved to %s\n", out_name);


	for (int i = 0; i < nimg; i++) fclose(fps[i]);
	free(fps);
	free(chunk_buf);
	free(raw_chunk_colrs);
	fclose(fout);

	return 0;
}


CMATRIX* init_basis_matrix(const char* fspec, int* out_w, int* out_h) {
	char full_path[1024];
	char workdir[1024];

	if (strstr(fspec, "%")) {
		if (getcwd(workdir, sizeof(workdir)) != NULL) {
#ifdef _WIN32
			snprintf(full_path, sizeof(full_path), "%s\\basis.bin", workdir);
#else
			snprintf(full_path, sizeof(full_path), "%s/basis.bin", workdir);
#endif

			if (hdr2bin(fspec, full_path) > 0) {
				fprintf(stderr, "Fatal: Cannot convert hdr images to binary file  %s\n", full_path);
				return 0;
			}

		}
	}
	else {
		strncpy(full_path, fspec, sizeof(full_path) - 1);
	}

	CMATRIX* basis = map_as_cmatrix(full_path, out_w, out_h);
	if (!basis) {
		fprintf(stderr, "Fatal: Cannot map basis matrix from %s\n", full_path);
	}
	return basis;
}


int*
stat_numn_onzerocol(CMATRIX* cv, int* size) {
	int K = cv->nrows;
	int C = cv->ncols;

	int* pixel_count_map = (int*)calloc(C, sizeof(int));
	int num_nonzero = 0;

	for (int c = 0; c < C; c++) {
		int col_pixel_count = 0;
		for (int k = 0; k < K; k++) {
			float* val = (float*)cm_lval(cv, k, c);
			if (val[0] != 0.0f || val[1] != 0.0f || val[2] != 0.0f) {
				col_pixel_count++;
			}
		}

		if (col_pixel_count > 0) {
			pixel_count_map[num_nonzero] = col_pixel_count;
			num_nonzero++;
		}
	}

	*size = num_nonzero;
	return pixel_count_map;
}

/* check to see if a string contains a %d or %o specification */
static int
hasNumberFormat(const char *s)
{
	if (s == NULL)
		return(0);

	while (*s) {
		while (*s != '%')
			if (!*s++)
				return(0);
		if (*++s == '%') {		/* ignore "%%" */
			++s;
			continue;
		}
		while (isdigit(*s))		/* field length */
			++s;
						/* field we'll use? */
		if ((*s == 'd') | (*s == 'i') | (*s == 'o') |
					(*s == 'x') | (*s == 'X'))
			return(1);
	}
	return(0);				/* didn't find one */
}

static int
sum_images_matrix(const char* fspec, const CMATRIX* cv, FILE* fout)
{
	static CMATRIX* cached_basis = NULL;
	static int cached_w = 0, cached_h = 0;
	static int is_initialized = 0;

	if (!is_initialized) {
		char full_path[1024];
		char workdir[1024];

		if (strstr(fspec, "%")) {
			if (getcwd(workdir, sizeof(workdir)) != NULL) {
#ifdef _WIN32
				snprintf(full_path, sizeof(full_path), "%s\\basis.bin", workdir);
#else
				snprintf(full_path, sizeof(full_path), "%s/basis.bin", workdir);
#endif

				if (hdr2bin(fspec, full_path) > 0) {
					fprintf(stderr, "Fatal: Cannot convert hdr images to binary file  %s\n", full_path);
					return 0;
				}

			}
		}
		else {
			strncpy(full_path, fspec, sizeof(full_path) - 1);
		}

		cached_basis = map_as_cmatrix(full_path, &cached_w, &cached_h);
		if (!cached_basis) {
			fprintf(stderr, "Fatal: Cannot map basis matrix from %s\n", full_path);
			return 0;
		}
		is_initialized = 1;
	}

	double checksum = 0;

	CMATRIX* pmat = cm_multiply_avx2(cached_basis, cv);
	if (!pmat)
	{
		printf("pmat is null\n");
		return 0;
	}


	int save_r = pmat->nrows;
	int save_c = pmat->ncols;
	pmat->nrows = cached_h;
	pmat->ncols = cached_w;

	fputformat((char*)cm_fmt_id[DTrgbe], fout);
	fputc('\n', fout);
	fflush(fout);

	int i = cm_write(pmat, DTrgbe, fout);


	pmat->nrows = save_r;
	pmat->ncols = save_c;

	cm_free_u(pmat);
	free_mmap_cmatrix(cached_basis, cached_w, cached_h, cached_basis->ncols);

	return i;
}

static int
sum_images_matrix_batch(CMATRIX* cached_basis,
	int cached_w,
	int cached_h,
	const char* fspec,
	const CMATRIX* cv_batch,
	int start_index,
	FILE* ofspec_base,
	const char* ofspec_pattern)
{
	CMATRIX* pmat_batch = cm_multiply_smart(cached_basis, cv_batch);
	if (!pmat_batch) {
		fprintf(stderr, "Matrix multiplication failed for batch starting at %d\n", start_index);
		return 0;
	}
	if (!pmat_batch) return 0;
	int batch_size = pmat_batch->ncols;
	for (int b = 0; b < batch_size; b++) {
		char fnbuf[1024];
		FILE* ofp = ofspec_base;
		int current_frame_id = start_index + b;
		int needs_fclose = 0;

		if (ofspec_pattern != NULL && hasNumberFormat(ofspec_pattern)) {
			sprintf(fnbuf, ofspec_pattern, current_frame_id);
			if ((ofp = fopen(fnbuf, "wb")) == NULL) continue;
			needs_fclose = 1;
			newheader("RADIANCE", ofp);
			fputnow(ofp);
			fputformat((char*)cm_fmt_id[DTrgbe], ofp);
			fprintf(ofp, "FRAME=%d\n", current_frame_id);
			fputc('\n', ofp);
			fflush(ofp);
		}

		CMATRIX* one_frame = cm_column_u(pmat_batch, b);
		int save_r = one_frame->nrows;
		int save_c = one_frame->ncols;
		one_frame->nrows = cached_h;
		one_frame->ncols = cached_w;

#pragma omp critical(hdr_encode_lock)
		{
			if (cm_write(one_frame, DTrgbe, ofp) < 0) {
				fprintf(stderr, "Error writing frame %d\n", current_frame_id);
			}
		}

		one_frame->nrows = save_r;
		one_frame->ncols = save_c;
		cm_free_u(one_frame);

		if (needs_fclose) {
			fclose(ofp);
		}
	}
	cm_free_u(pmat_batch);
	//free_mmap_cmatrix(cached_basis, cached_w, cached_h, cached_basis->ncols);
	return 1;
}

#endif /* INTEL_DCTOPT */

/* Sum together a set of images and write result to fout */
static int
sum_images(const char *fspec, const CMATRIX *cv, FILE *fout)
{
	static int	runcnt = 0;
	int	myDT = DTfromHeader;
	COLR	*scanline = NULL;
	CMATRIX	*pmat = NULL;
	int	myXR=0, myYR=0;
	int	i, y;

	if (cv->ncols != 1)
		error(INTERNAL, "expected vector in sum_images()");
	for (i = cv->nrows; i-- > 0; ) {
		const int	r = runcnt&1 ? i : cv->nrows-1 - i;
		const COLORV	*scv = cv_lval(cv,r);
		int		flat_file = 0;
		char		fname[1024];
		FILE		*fp;
		long		data_start;
		int		dt, xr, yr;
		COLORV		*psp;
		char		*err;
							/* check for zero */
		if ((scv[RED] == 0) & (scv[GRN] == 0) & (scv[BLU] == 0) &&
				(myDT != DTfromHeader) | (i > 0))
			continue;
							/* open next picture */
		sprintf(fname, fspec, r);
		if ((fp = fopen(fname, "rb")) == NULL) {
			sprintf(errmsg, "cannot open picture '%s'", fname);
			error(SYSTEM, errmsg);
		}
#ifdef getc_unlocked
		flockfile(fp);
#endif
		dt = DTfromHeader;
		if ((err = cm_getheader(&dt, NULL, NULL, NULL, NULL, fp)) != NULL)
			error(USER, err);
		if ((dt != DTrgbe) & (dt != DTxyze) ||
				!fscnresolu(&xr, &yr, fp)) {
			sprintf(errmsg, "file '%s' not a picture", fname);
			error(USER, errmsg);
		}
		if (myDT == DTfromHeader) {		/* on first one */
			myDT = dt;
			myXR = xr; myYR = yr;
			scanline = (COLR *)malloc(sizeof(COLR)*myXR);
			if (scanline == NULL)
				error(SYSTEM, "out of memory in sum_images()");
			pmat = cm_alloc(myYR, myXR);
			memset(pmat->cmem, 0, sizeof(COLOR)*myXR*myYR);
							/* finish header */
			fputformat(cm_fmt_id[myDT], fout);
			fputc('\n', fout);
			fflush(fout);
		} else if ((dt != myDT) | (xr != myXR) | (yr != myYR)) {
			sprintf(errmsg, "picture '%s' format/size mismatch",
					fname);
			error(USER, errmsg);
		}
							/* flat file check */
		if ((data_start = ftell(fp)) > 0 && fseek(fp, 0L, SEEK_END) == 0) {
			flat_file = (ftell(fp) >= data_start + sizeof(COLR)*xr*yr);
			if (fseek(fp, data_start, SEEK_SET) < 0) {
				sprintf(errmsg, "cannot seek on picture '%s'", fname);
				error(SYSTEM, errmsg);
			}
		}
		psp = pmat->cmem;
		for (y = 0; y < yr; y++) {		/* read it in */
			COLOR	col;
			int	x;
			if (flat_file ? getbinary(scanline, sizeof(COLR), xr, fp) != xr :
					freadcolrs(scanline, xr, fp) < 0) {
				sprintf(errmsg, "error reading picture '%s'",
						fname);
				error(SYSTEM, errmsg);
			}
							/* sum in scanline */
			for (x = 0; x < xr; x++, psp += 3) {
				if (!scanline[x][EXP])
					continue;	/* skip zeroes */
				colr_color(col, scanline[x]);
				multcolor(col, scv);
				addcolor(psp, col);
			}
		}
		fclose(fp);				/* done this picture */
	}
	free(scanline);
	i = cm_write(pmat, myDT, fout);			/* write picture */
	cm_free(pmat);					/* free data */
	++runcnt;					/* for zig-zagging */
	return(i);
}

/* adjust matrix dimensions according to user size(s) */
static int
alt_dim(CMATRIX *cm, int nr, int nc)
{
	if ((nr <= 0) & (nc <= 0))
		return(0);
	if ((nr == cm->nrows) & (nc == cm->ncols))
		return(0);
	if (nr > 0) {
		if (nc <= 0)
			nc = cm->nrows*cm->ncols/nr;
		if (nr*nc != cm->nrows*cm->ncols) {
			fprintf(stderr, "Bad dimensions: %dx%d != %dx%d\n",
					nr, nc, cm->nrows, cm->ncols);
			return(-1);
		}
	} else /* nc > 0 */ {
		nr = cm->nrows*cm->ncols/nc;
		if (nc*nr != cm->nrows*cm->ncols) {
			fprintf(stderr, "Bad dimensions: %d does not divide %dx%d evenly\n",
					nc, cm->nrows, cm->ncols);
			return(-1);
		}
	}
	cm->nrows = nr;
	cm->ncols = nc;
	return(1);
}

#ifndef INTEL_DCTOPT
/* check to see if a string contains a %d or %o specification */
static int
hasNumberFormat(const char *s)
{
	if (s == NULL)
		return(0);
#endif /* not INTEL_DCTOPT */

#ifndef INTEL_DCTOPT
	while (*s) {
		while (*s != '%')
			if (!*s++)
				return(0);
		if (*++s == '%') {		/* ignore "%%" */
			++s;
			continue;
		}
		while (isdigit(*s))		/* field length */
			++s;
						/* field we'll use? */
		if ((*s == 'd') | (*s == 'i') | (*s == 'o') |
					(*s == 'x') | (*s == 'X'))
			return(1);
	}
	return(0);				/* didn't find one */
}
#endif /* not INTEL_DCTOPT */

int
main(int argc, char *argv[])
{
	int		skyfmt = DTfromHeader;
	int		outfmt = DTascii;
	int		headout = 1;
	int		nsteps = 0;
	char		*ofspec = NULL;
	FILE		*ofp = stdout;
	int		xres=0, yres=0;
	CMATRIX		*cmtx;		/* component vector/matrix result */
	char		fnbuf[256];
	int		a, i;
					/* set global progname */
	fixargv0(argv[0]);
					/* get options */
	for (a = 1; a < argc && argv[a][0] == '-'; a++)
		switch (argv[a][1]) {
		case 'n':
			nsteps = atoi(argv[++a]);
			if (nsteps < 0)
				goto userr;
			skyfmt = nsteps ? DTascii : DTfromHeader;
			break;
		case 'h':
			headout = !headout;
			break;
		case 'i':
			switch (argv[a][2]) {
			case 'f':
				skyfmt = DTfloat;
				break;
			case 'd':
				skyfmt = DTdouble;
				break;
			case 'a':
				skyfmt = DTascii;
				break;
			default:
				goto userr;
			}
			break;
		case 'o':
			switch (argv[a][2]) {
			case '\0':	/* output specification (not format) */
				ofspec = argv[++a];
#ifdef INTEL_DCTOPT
				//ofspec = "results/2ph/%04d.hdr";
#endif /* INTEL_DCTOPT */
				break;
			case 'f':
				outfmt = DTfloat;
				break;
			case 'd':
				outfmt = DTdouble;
				break;
			case 'a':
				outfmt = DTascii;
				break;
			case 'c':
				outfmt = DTrgbe;
				break;
			default:
				goto userr;
			}
			break;
		case 'x':
			xres = atoi(argv[++a]);
			break;
		case 'y':
			yres = atoi(argv[++a]);
			break;
		default:
			goto userr;
		}
	if ((argc-a < 1) | (argc-a > 4))
		goto userr;

#ifdef INTEL_DCTOPT


#endif /* INTEL_DCTOPT */
	if (argc-a > 2) {			/* VTDs expression */
		CMATRIX		*smtx, *Dmat, *Tmat, *imtx;
		const char	*ccp;
						/* get sky vector/matrix */
#ifndef INTEL_DCTOPT
		smtx = cm_load(argv[a+3], 0, nsteps, skyfmt);
#else /* INTEL_DCTOPT */
		smtx = cm_load_u(argv[a+3], 0, nsteps, skyfmt);
#endif /* INTEL_DCTOPT */
		nsteps = smtx->ncols;
						/* load BSDF */
		if (argv[a+1][0] != '!' &&
				(ccp = strrchr(argv[a+1], '.')) > argv[a+1] &&
				!strcasecmp(ccp+1, "XML"))
			Tmat = cm_loadBTDF(argv[a+1]);
		else
#ifndef INTEL_DCTOPT
			Tmat = cm_load(argv[a+1], 0, 0, DTfromHeader);
#else /* INTEL_DCTOPT */
			Tmat = cm_load_u(argv[a+1], 0, 0, DTfromHeader);
#endif /* INTEL_DCTOPT */
						/* load Daylight matrix */
#ifndef INTEL_DCTOPT
		Dmat = cm_load(argv[a+2], Tmat->ncols,
#else /* INTEL_DCTOPT */
		Dmat = cm_load_u(argv[a+2], Tmat->ncols,
#endif /* INTEL_DCTOPT */
					smtx->nrows, DTfromHeader);
						/* multiply vector through */
#ifndef INTEL_DCTOPT
		imtx = cm_multiply(Dmat, smtx);
		cm_free(Dmat); cm_free(smtx);
		cmtx = cm_multiply(Tmat, imtx);
#else /* INTEL_DCTOPT */
		imtx = cm_multiply_smart(Dmat, smtx);
		cm_free_u(Dmat); cm_free_u(smtx);
		cmtx = cm_multiply_smart(Tmat, imtx);
#endif /* INTEL_DCTOPT */
		cm_free(Tmat); 
#ifndef INTEL_DCTOPT
		cm_free(imtx);
#else /* INTEL_DCTOPT */
		cm_free_u(imtx);
#endif /* INTEL_DCTOPT */
	} else {				/* sky vector/matrix only */
#ifndef INTEL_DCTOPT
		cmtx = cm_load(argv[a+1], 0, nsteps, skyfmt);
#else /* INTEL_DCTOPT */
		cmtx = cm_load_u(argv[a+1], 0, nsteps, skyfmt);
#endif /* INTEL_DCTOPT */
		nsteps = cmtx->ncols;
	}
						/* prepare output stream */
	if ((ofspec != NULL) & (nsteps == 1) && hasNumberFormat(ofspec)) {
		sprintf(fnbuf, ofspec, 0);
		ofspec = fnbuf;
	}
	if (ofspec != NULL && !hasNumberFormat(ofspec)) {
		if ((ofp = fopen(ofspec, "w")) == NULL) {
			fprintf(stderr, "%s: cannot open '%s' for output\n",
					progname, ofspec);
			return(1);
		}
		ofspec = NULL;			/* only need to open once */
	}
#ifdef INTEL_DCTOPT
	//argv[a] = "matrices/dcm3/%04d.hdr";
#endif /* INTEL_DCTOPT */
	if (hasNumberFormat(argv[a])) {		/* loading image vector(s) */
		if (outfmt != DTrgbe) {
			error(WARNING, "changing output type to -oc");
			outfmt = DTrgbe;
		}
		if (ofspec == NULL) {
			SET_FILE_BINARY(ofp);
			newheader("RADIANCE", ofp);
			printargs(argc, argv, ofp);
			fputnow(ofp);
		}
#ifndef INTEL_DCTOPT
		if (nsteps > 1)			/* multiple output frames? */
#else /* INTEL_DCTOPT */


		if (nsteps > 1)	{	/* multiple output frames? */

			int batch_size = 16;

			char full_path[1024];
			char workdir[1024];
			CMATRIX* cached_basis = NULL;
			int cached_w = 0, cached_h = 0;
			int size = 0;
			int* nonzeronum_percols = stat_numn_onzerocol(cmtx, &size);
			int allone = 1;
			for (int k = 0;k < size;k++)
			{
				if (*(nonzeronum_percols + k) != 1)
				{
					allone = 0;
					break;
				}
			}

			if (allone == 0)
			{
				cached_basis = init_basis_matrix(argv[a], &cached_w, &cached_h);
				if (cached_basis == NULL) {
					error(SYSTEM, "Failed to convert HDRs to Binary file");
					return 1;
				}
				#pragma omp parallel for private(i)
				for (i = 0; i < nsteps; i += batch_size) {

					FILE* ofp1 = NULL;
					int current_batch = (i + batch_size > nsteps) ? (nsteps - i) : batch_size;
					CMATRIX* cvec_batch = cm_get_batch(cmtx, i, current_batch);

					if (!sum_images_matrix_batch(cached_basis, cached_w, cached_h, argv[a], cvec_batch, i, ofp1, ofspec)) {
						error(SYSTEM, "Batch matrix calculation failed");

					}
					
					cm_free_u(cvec_batch);
				}

				free_mmap_cmatrix(cached_basis, cached_w, cached_h, cached_basis->ncols);

			} else {

#endif /* INTEL_DCTOPT */
			for (i = 0; i < nsteps; i++) {
				CMATRIX	*cvec = cm_column(cmtx, i);
				if (ofspec != NULL) {
					sprintf(fnbuf, ofspec, i);
					if ((ofp = fopen(fnbuf, "wb")) == NULL) {
						fprintf(stderr,
					"%s: cannot open '%s' for output\n",
							progname, fnbuf);
						return(1);
					}
					newheader("RADIANCE", ofp);
					printargs(argc, argv, ofp);
					fputnow(ofp);
				}
				fprintf(ofp, "FRAME=%d\n", i);
				if (!sum_images(argv[a], cvec, ofp))
					return(1);
				if (ofspec != NULL) {
					if (fclose(ofp) == EOF) {
						fprintf(stderr,
						"%s: error writing to '%s'\n",
							progname, fnbuf);
						return(1);
					}
					ofp = stdout;
				}
				cm_free(cvec);
			}
#ifndef INTEL_DCTOPT
		else if (!sum_images(argv[a], cmtx, ofp))
			return(1);
#else /* INTEL_DCTOPT */
			
			}
		}
		else {
				int size = 0;
				int* nonzeronum_percols = stat_numn_onzerocol(cmtx, &size);
				int allone = 1;
				for (int k = 0;k < size;k++)
				{
					if (*(nonzeronum_percols + k) != 1)
					{
						allone = 0;
						break;
					}
				}

				if (allone == 0)
				{
					if (!sum_images_matrix(argv[a], cmtx, ofp)) return(1);
				}
				else
				{
					if (!sum_images(argv[a], cmtx, ofp)) return (1);
				}
			}
#endif /* INTEL_DCTOPT */
	} else {				/* loading view matrix */
#ifndef INTEL_DCTOPT
		CMATRIX	*Vmat = cm_load(argv[a], 0, cmtx->nrows, DTfromHeader);
		CMATRIX	*rmtx = cm_multiply(Vmat, cmtx);
		cm_free(Vmat);
#else /* INTEL_DCTOPT */
		CMATRIX	*Vmat = cm_load_u(argv[a], 0, cmtx->nrows, DTfromHeader);
		CMATRIX	*rmtx = cm_multiply_smart(Vmat, cmtx);
		cm_free_u(Vmat);
#endif /* INTEL_DCTOPT */
		if (ofspec != NULL) {		/* multiple vector files? */
			const char	*wtype = (outfmt==DTascii) ? "w" : "wb";
			for (i = 0; i < nsteps; i++) {
#ifndef INTEL_DCTOPT
				CMATRIX	*rvec = cm_column(rmtx, i);
#else /* INTEL_DCTOPT */
				CMATRIX	*rvec = cm_column_u(rmtx, i);
#endif /* INTEL_DCTOPT */
				if (alt_dim(rvec, yres, xres) < 0)
					return(1);
				sprintf(fnbuf, ofspec, i);
				if ((ofp = fopen(fnbuf, wtype)) == NULL) {
					fprintf(stderr,
					"%s: cannot open '%s' for output\n",
							progname, fnbuf);
					return(1);
				}
#ifdef getc_unlocked
				flockfile(ofp);
#endif
				if (headout) {	/* header output */
					newheader("RADIANCE", ofp);
					printargs(argc, argv, ofp);
					fputnow(ofp);
					fprintf(ofp, "FRAME=%d\n", i);
					if ((outfmt != DTrgbe) & (outfmt != DTxyze)) {
						fprintf(ofp, "NROWS=%d\n", rvec->nrows);
						fprintf(ofp, "NCOLS=%d\n", rvec->ncols);
						fputncomp(3, ofp);
					}
					if ((outfmt == DTfloat) | (outfmt == DTdouble))
						fputendian(ofp);
					fputformat(cm_fmt_id[outfmt], ofp);
					fputc('\n', ofp);
				}
				cm_write(rvec, outfmt, ofp);
				if (fclose(ofp) == EOF) {
					fprintf(stderr,
						"%s: error writing to '%s'\n",
							progname, fnbuf);
					return(1);
				}
				ofp = stdout;
				cm_free(rvec);
			}
		} else {
#ifdef getc_unlocked
			flockfile(ofp);
#endif
			if (outfmt != DTascii)
				SET_FILE_BINARY(ofp);
			if (alt_dim(rmtx, yres, xres) < 0)
				return(1);
			if (headout) {		/* header output */
				newheader("RADIANCE", ofp);
				printargs(argc, argv, ofp);
				fputnow(ofp);
				if ((outfmt != DTrgbe) & (outfmt != DTxyze)) {
					fprintf(ofp, "NROWS=%d\n", rmtx->nrows);
					fprintf(ofp, "NCOLS=%d\n", rmtx->ncols);
					fputncomp(3, ofp);
				}
				if ((outfmt == DTfloat) | (outfmt == DTdouble))
					fputendian(ofp);
				fputformat(cm_fmt_id[outfmt], ofp);
				fputc('\n', ofp);
			}
			cm_write(rmtx, outfmt, ofp);
		}
#ifndef INTEL_DCTOPT
		cm_free(rmtx);
#else /* INTEL_DCTOPT */
		cm_free_u(rmtx);
#endif /* INTEL_DCTOPT */
	}
	if (fflush(ofp) == EOF) {		/* final clean-up */
		fprintf(stderr, "%s: write error on output\n", progname);
		return(1);
	}
#ifndef INTEL_DCTOPT
	cm_free(cmtx);
#else /* INTEL_DCTOPT */
	cm_free_u(cmtx);
#endif /* INTEL_DCTOPT */
	return(0);
userr:
	fprintf(stderr, "Usage: %s [-n nsteps][-o ospec][-x xr][-y yr][-i{f|d|h}][-o{f|d|c}] DCspec [skyf]\n",
				progname);
	fprintf(stderr, "   or: %s [-n nsteps][-o ospec][-x xr][-y yr][-i{f|d|h}][-o{f|d|c}] Vspec Tbsdf Dmat.dat [skyf]\n",
				progname);
	return(1);
}
