#ifndef __HELPERS__
#define __HELPERS__

#include <math.h>
#include <float.h>
#include "utils.h"

void vec_set_zero(real_t *v, int n);
void vec_copy(real_t *dst, const real_t *src, int n);
real_t dot(const real_t *a, const real_t *b, int n);
void axpy(real_t *y, real_t a, const real_t *x, int n);
void matvec_dense(const real_t *A, const real_t *x, real_t *y, int n);
real_t norm_inf_diff(const real_t *a, const real_t *b, int n);
real_t norm2(const real_t *v, int n);

#endif