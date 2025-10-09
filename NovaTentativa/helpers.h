#ifndef HELPERS_H
#define HELPERS_H

#include "utils.h"  // define real_t

/** Preenche o vetor v com zeros (v[:] = 0). */
void vector_fill_zero(real_t *v, int n);

/** Copia src -> dst (tamanho n). */
void vector_copy(real_t *dst, const real_t *src, int n);

/** Produto interno (dot) entre 'a' e 'b'. */
real_t dot_product(const real_t *a, const real_t *b, int n);

/** y := y + a * x (AXPY). */
void vector_axpy(real_t *y, real_t a, const real_t *x, int n);

/** y = A * x (matriz densa A armazenada por linhas). */
void matvec_multiply_dense(const real_t *A, const real_t *x, real_t *y, int n);

/** Norma infinito da diferença: ||a - b||_inf. */
real_t norm_inf_difference(const real_t *a, const real_t *b, int n);

/** Norma L2 (euclidiana): ||v||_2. */
real_t norm_l2(const real_t *v, int n);

#endif // HELPERS_H
