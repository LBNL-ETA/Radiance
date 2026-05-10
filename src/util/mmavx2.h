#ifndef MMAVX2_H
#define MMAVX2_H
#include "cmatrix.h"
#include <stdint.h>

extern CMATRIX* cm_multiply_avx2(const CMATRIX* cm1, const CMATRIX* cm2);
extern CMATRIX* obtain_scanline_matrix(NativeHandle* handles, int64_t** index_table, int start_y, int num_hdrs, int w, int num_lines);
#endif
