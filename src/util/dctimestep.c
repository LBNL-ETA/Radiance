#ifndef lint
static const char RCSid[] = "$Id$";
#endif
/*
 * Compute time-step result using Daylight Coefficient method.
 *
 *	G. Ward
 */
#if defined(_WIN32) || defined(_WIN64)
#include <windows.h>
#include <psapi.h>
#else
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/resource.h>
#endif
#include <ctype.h>
#include <stdlib.h>
#include "platform.h"
#include "standard.h"
#include "cmatrix.h"
#include "resolu.h"
#if defined(__AVX2__)
#include <immintrin.h>
#endif
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "color.h"
#include "cmmap.h"
#include "cmatrixavx2.h"
#include "threadpool.h"

int MAX_CHUNK_SIZE;
float exponent_table[256];

typedef struct {
	int task_index;
	int i_offset;
	int lines_to_read;
	int max_file_num;
	int cached_w;
	int total_images;
	int64_t** index_table;
	NativeHandle* handles;
	CMATRIX* smtx;
	CMATRIX* result_chunk;
	mtx_t sync_mtx;        
	cnd_t sync_cnd;        
	int is_done;           
} ImageTaskData;


void compute(ImageTaskData* data) {
	CMATRIX* partial_dc = obtain_scanline_matrix(data->handles, data->index_table, data->i_offset, data->total_images, data->cached_w, data->lines_to_read);
	//CMATRIX* partial_dc = obtain_scanline_matrix(data->index_table, data->i_offset, data->total_images, data->cached_w, data->lines_to_read);

	if (partial_dc == NULL) {
		data->result_chunk = NULL;
		return;
	}

	data->result_chunk = cm_multiply_smart(partial_dc, data->smtx);
	cm_free_u(partial_dc);
}

void compute_matrix_task(void* arg) {
	ImageTaskData* data = (ImageTaskData*)arg;
	compute(data);
	mtx_lock(&data->sync_mtx);
	data->is_done = 1;
	cnd_signal(&data->sync_cnd);
	mtx_unlock(&data->sync_mtx);
}


#ifndef _WIN32
void set_max_files(int num_files) {
    struct rlimit rl;

    if (getrlimit(RLIMIT_NOFILE, &rl) == -1) {
        perror("getrlimit");
        return;
    }


    if (num_files > rl.rlim_max) {
    
        
        rl.rlim_cur = rl.rlim_max;
        fprintf(stderr, "Warning: Requested %d files, but system hard limit is %llu. Capped.\n", 
                num_files, (unsigned long long)rl.rlim_max);
    } else {
      
        rl.rlim_cur = num_files;
    }

 
    if (setrlimit(RLIMIT_NOFILE, &rl) == -1) {
        perror("setrlimit");
    }
}
#endif

/* check to see if a string contains a %d or %o specification */
static int
hasNumberFormat(const char* s)
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



/* Sum together a set of images and write result to fout */
static int
sum_images(const char* fspec, const CMATRIX* cv, FILE* fout)
{
	static int	runcnt = 0;
	int	myDT = DTfromHeader;
	COLR* scanline = NULL;
	CMATRIX* pmat = NULL;
	int	myXR = 0, myYR = 0;
	int	i, y;

	if (cv->ncols != 1)
		error(INTERNAL, "expected vector in sum_images()");
	for (i = cv->nrows; i-- > 0; ) {
		const int	r = runcnt & 1 ? i : cv->nrows - 1 - i;
		const COLORV* scv = cv_lval(cv, r);
		int		flat_file = 0;
		char		fname[1024];
		FILE* fp;
		long		data_start;
		int		dt, xr, yr;
		COLORV* psp;
		char* err;
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
			scanline = (COLR*)malloc(sizeof(COLR) * myXR);
			if (scanline == NULL)
				error(SYSTEM, "out of memory in sum_images()");
			pmat = cm_alloc(myYR, myXR);
			memset(pmat->cmem, 0, sizeof(COLOR) * myXR * myYR);
			/* finish header */
			fputformat(cm_fmt_id[myDT], fout);
			fputc('\n', fout);
			fflush(fout);
		}
		else if ((dt != myDT) | (xr != myXR) | (yr != myYR)) {
			sprintf(errmsg, "picture '%s' format/size mismatch",
				fname);
			error(USER, errmsg);
		}
		/* flat file check */
		if ((data_start = ftell(fp)) > 0 && fseek(fp, 0L, SEEK_END) == 0) {
			flat_file = (ftell(fp) >= data_start + sizeof(COLR) * xr * yr);
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
alt_dim(CMATRIX* cm, int nr, int nc)
{
	if ((nr <= 0) & (nc <= 0))
		return(0);
	if ((nr == cm->nrows) & (nc == cm->ncols))
		return(0);
	if (nr > 0) {
		if (nc <= 0)
			nc = cm->nrows * cm->ncols / nr;
		if (nr * nc != cm->nrows * cm->ncols) {
			fprintf(stderr, "Bad dimensions: %dx%d != %dx%d\n",
				nr, nc, cm->nrows, cm->ncols);
			return(-1);
		}
	}
	else /* nc > 0 */ {
		nr = cm->nrows * cm->ncols / nc;
		if (nc * nr != cm->nrows * cm->ncols) {
			fprintf(stderr, "Bad dimensions: %d does not divide %dx%d evenly\n",
				nc, cm->nrows, cm->ncols);
			return(-1);
		}
	}
	cm->nrows = nr;
	cm->ncols = nc;
	return(1);
}



int
main(int argc, char* argv[])
{
	int		skyfmt = DTfromHeader;
	int		outfmt = DTascii;
	int		headout = 1;
	int		nsteps = 0;
	char* ofspec = NULL;
	FILE* ofp = stdout;
	int		xres = 0, yres = 0;
	CMATRIX* cmtx;		/* component vector/matrix result */
	char		fnbuf[256];
	int		a, i;
	int			r, c;
	const COLORV* mp;
	double mCacheGB = (double)(get_available_memory())/(1UL*1024*1024*1024);
	int nprcoess = 4;
	

#if defined(_WIN32)
	#define MAX_OPEN_FILES 8000
	_setmaxstdio(MAX_OPEN_FILES+3);
#else
	#define MAX_OPEN_FILES 20000
	set_max_files(MAX_OPEN_FILES+3);
#endif

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
		case 'N':
			nprcoess = atoi(argv[++a]);
			break;
		case 'm':
			mCacheGB = atof(argv[++a]);
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
				//ofspec = "results/2ph/%04d.hdr";
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

	init_exponent_table();

	if ((argc - a < 1) | (argc - a > 4))
		goto userr;

	if (argc - a > 2) {			/* VTDs expression */
		CMATRIX* smtx, * Dmat, * Tmat, * imtx;
		const char* ccp;
		/* get sky vector/matrix */
		smtx = cm_load_u(argv[a + 3], 0, nsteps, skyfmt);
		nsteps = smtx->ncols;
		/* load BSDF */
		if (argv[a + 1][0] != '!' &&
			(ccp = strrchr(argv[a + 1], '.')) > argv[a + 1] &&
			!strcasecmp(ccp + 1, "XML"))
			Tmat = cm_loadBTDF(argv[a + 1]);
		else
			Tmat = cm_load_u(argv[a + 1], 0, 0, DTfromHeader);
		/* load Daylight matrix */
		Dmat = cm_load_u(argv[a + 2], Tmat->ncols,
			smtx->nrows, DTfromHeader);
		/* multiply vector through */
		imtx = cm_multiply_smart(Dmat, smtx);
		cm_free_u(Dmat); cm_free_u(smtx);
		cmtx = cm_multiply_smart(Tmat, imtx);
		cm_free(Tmat);
		cm_free_u(imtx);
	}
	else {				/* sky vector/matrix only */
		cmtx = cm_load_u(argv[a + 1], 0, nsteps, skyfmt);
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
	//argv[a] = "dc/%04d.hdr";
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

		if (nsteps > 1) {	/* multiple output frames? */

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
				int cached_w = 0, cached_h = 0;
				int total_images = 0;
			
				MappedFile* map = mapmultiHDRs(argv[a], &cached_w, &cached_h, &total_images);

				if (map == NULL)
				{
					error(SYSTEM, "Failed to map cretae mapped files\n");
					return 1;
				}

				int num_cores = nprcoess;
				thread_pool_t* pool = pool_create(num_cores);

				int64_t ** index_table = build_scanline_index(map, total_images, cached_w, cached_h, pool);
				
				unmapmultiHDRs(map, total_images);
				
				NativeHandle* handles = read_multihdrs(argv[a], total_images);

				unsigned long long available_mem = mCacheGB * (1UL * 1024 * 1024 * 1024);
				unsigned long long MAX_MEMORY = available_mem / num_cores * 0.9;
	

				int batch_line = MAX_MEMORY / (MAX_OPEN_FILES * cached_w * 3 * sizeof(float));

				if (batch_line > 4)
				{
					batch_line = 4;
				}
				else if (batch_line < 1)
				{
					batch_line = 1;
					fprintf(stderr, "The specified RAM may be not enough, please increase the spcified RAM or decrease the number of the processes\n");
				}
				
				fprintf(stderr, "batch_line is %d\n", batch_line);
				int flush_counter = 0;
				for (int j = 0; j < nsteps; j += MAX_OPEN_FILES) {
					int max_file_num = MAX_OPEN_FILES;
					if (j + MAX_OPEN_FILES > nsteps) max_file_num = nsteps - j;

					FILE* out_files[MAX_OPEN_FILES];
					for (int k = 0; k < max_file_num; k++) {
						char path[256];
						sprintf(path, ofspec, j + k);
						out_files[k] = fopen(path, "wb");
						if (out_files[k] == NULL) {
							fprintf(stderr, "%s: cannot open '%s' for output\n", progname, path);
						}
						else {
							setvbuf(out_files[k], NULL, _IOFBF, 1024 * 1024);
						}
					}

					CMATRIX* smtx = cm_get_batch(cmtx, j, max_file_num);

					int total_tasks = (cached_h + batch_line - 1) / batch_line;
					ImageTaskData* tasks = (ImageTaskData*)malloc(total_tasks * sizeof(ImageTaskData));

					int task_idx = 0;
					for (int i = 0; i < cached_h; i += batch_line, task_idx++) {
						int lines_to_read = batch_line;
						if (i + batch_line > cached_h) lines_to_read = cached_h - i;

						tasks[task_idx].task_index = task_idx;
						tasks[task_idx].i_offset = i;
						tasks[task_idx].lines_to_read = lines_to_read;
						tasks[task_idx].max_file_num = max_file_num;
						tasks[task_idx].cached_w = cached_w;
						tasks[task_idx].total_images = total_images;
						tasks[task_idx].index_table = index_table;
						tasks[task_idx].smtx = smtx;
						tasks[task_idx].handles = handles;

						mtx_init(&tasks[task_idx].sync_mtx, mtx_plain);
						cnd_init(&tasks[task_idx].sync_cnd);
						tasks[task_idx].is_done = 0;
					}


					int max_in_flight = nprcoess; 
					int in_flight = 0;     
					int submit_idx = 0;   

					
					for (int process_idx = 0; process_idx < total_tasks; process_idx++) {
						while (submit_idx < total_tasks && in_flight < max_in_flight) {
							pool_enqueue(pool, compute_matrix_task, &tasks[submit_idx]);
							submit_idx++;
							in_flight++; 
						}

						
						mtx_lock(&tasks[process_idx].sync_mtx);
						while (tasks[process_idx].is_done == 0) {
							cnd_wait(&tasks[process_idx].sync_cnd, &tasks[process_idx].sync_mtx);
						}
						CMATRIX* result_chunk = tasks[process_idx].result_chunk;
						mtx_unlock(&tasks[process_idx].sync_mtx);

						if (result_chunk != NULL) {
							if(process_idx % 100 == 0)
							{
								fprintf(stderr,"process %d is computing...\n",process_idx);
							}
							int lines_to_read = tasks[process_idx].lines_to_read;
							int pixels_per_frame = lines_to_read * cached_w;
							COLORV* thread_frame_buffer = (COLORV*)malloc(pixels_per_frame * 3 * sizeof(COLORV));
							

							int num_cols = result_chunk->ncols;
							COLORV* src_mem = result_chunk->cmem;

							for (int k = 0; k < max_file_num; k++) {
								FILE* ofp1 = out_files[k];
								if (process_idx == 0 && ofspec != NULL && hasNumberFormat(ofspec)) {
									newheader("RADIANCE", ofp1);
									fprintf(ofp1, "FRAME=%d\n", j + k);
									printargs(argc, argv, ofp1);
									fputnow(ofp1);
									fputformat((char*)cm_fmt_id[DTrgbe], ofp1);
									fputc('\n', ofp1);
									fflush(ofp1);
									fprtresolu(cached_w, cached_h, ofp1);
								}

								for (int p = 0; p < pixels_per_frame; p++) {
									thread_frame_buffer[p * 3 + 0] = src_mem[(p * num_cols + k) * 3 + 0];
									thread_frame_buffer[p * 3 + 1] = src_mem[(p * num_cols + k) * 3 + 1];
									thread_frame_buffer[p * 3 + 2] = src_mem[(p * num_cols + k) * 3 + 2];
								}

								
								COLORV* mp = thread_frame_buffer;
								for (int r = 0; r < lines_to_read; r++, mp += 3 * cached_w) {
									if (fwritescan((COLOR*)mp, cached_w, ofp1) < 0) {
										fprintf(stderr, "Error writing frame %d\n", k);
									}
								}
							}

							free(thread_frame_buffer);
							cm_free_u(result_chunk);
						}
						else {
							
							fprintf(stderr, "Warning: Task %d failed to compute\n", process_idx);
						}

						
						mtx_destroy(&tasks[process_idx].sync_mtx);
						cnd_destroy(&tasks[process_idx].sync_cnd);
						in_flight--; 
					} 

					
					free(tasks);
					cm_free_u(smtx);

					for (int k = 0; k < max_file_num; k++) {
						if (out_files[k]) {
							fclose(out_files[k]);
						}
				
					}
				}
				close_hdrfiles(handles, total_images);
				pool_destroy(pool);
				fprintf(stderr,"Computation is finished\n");
			}
			else {

				for (i = 0; i < nsteps; i++) {
					CMATRIX* cvec = cm_column(cmtx, i);
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

			}
		}
		else {
			if (!sum_images(argv[a], cmtx, ofp)) return (1);
		}
	}
	else {				/* loading view matrix */
		CMATRIX* Vmat = cm_load_u(argv[a], 0, cmtx->nrows, DTfromHeader);
		CMATRIX* rmtx = cm_multiply_smart(Vmat, cmtx);
		cm_free_u(Vmat);
		if (ofspec != NULL) {		/* multiple vector files? */
			const char* wtype = (outfmt == DTascii) ? "w" : "wb";
			for (i = 0; i < nsteps; i++) {
				CMATRIX* rvec = cm_column_u(rmtx, i);
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
		}
		else {
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
		cm_free_u(rmtx);
	}
	if (fflush(ofp) == EOF) {		/* final clean-up */
		fprintf(stderr, "%s: write error on output\n", progname);
		return(1);
	}
	cm_free_u(cmtx);
	return(0);
userr:
	fprintf(stderr, "Usage: %s [-n nsteps][-o ospec][-x xr][-y yr][-i{f|d|h}][-o{f|d|c}] DCspec [skyf]\n",
		progname);
	fprintf(stderr, "   or: %s [-n nsteps][-o ospec][-x xr][-y yr][-i{f|d|h}][-o{f|d|c}] Vspec Tbsdf Dmat.dat [skyf]\n",
		progname);
	return(1);
}
