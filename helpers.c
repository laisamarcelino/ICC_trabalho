#include <math.h>
#include <float.h>
#include "utils.h"
#include "helpers.h"
// ---------- helpers básicos ----------

// Setar os valores de um vetor como 0
void vec_set_zero(real_t *v, int n) {
  for (int i = 0; i < n; ++i) 
    v[i] = 0.0;
}

// Copiar os valores de um vetor para outro vetor
void vec_copy(real_t *dst, const real_t *src, int n) {
  for (int i = 0; i < n; ++i) 
    dst[i] = src[i];
}

// Produto interno de um vetor A por um vetor B
real_t dot(const real_t *a, const real_t *b, int n) {
  real_t s = 0.0;
  for (int i = 0; i < n; ++i) 
    s += a[i]*b[i];
  return s;
}

// Soma y := y + a*x (combinação linear de vetores)
void axpy(real_t *y, real_t a, const real_t *x, int n) {
  for (int i = 0; i < n; ++i) 
    y[i] += a * x[i];
}

// Multiplica matriz densa A por vetor x, resultado em y (y = A*x)
void matvec_dense(const real_t *A, const real_t *x, real_t *y, int n) {
  for (int i = 0; i < n; ++i) {
    real_t s = 0.0;
    const real_t *Ai = &A[i*(size_t)n];
    for (int j = 0; j < n; ++j) 
      s += Ai[j] * x[j];
    y[i] = s;
  }
}

// Calcula a norma infinito da diferença entre dois vetores (máximo valor absoluto das diferenças)
real_t norm_inf_diff(const real_t *a, const real_t *b, int n) {
  real_t m = 0.0;
  for (int i = 0; i < n; ++i) {
    real_t d = fabs(a[i] - b[i]);
    if (d > m) 
      m = d;
  }
  return m;
}

// Calcula a norma 2 (euclidiana) de um vetor (raiz da soma dos quadrados dos elementos)
real_t norm2(const real_t *v, int n) {
  return sqrt(dot(v, v, n));
}