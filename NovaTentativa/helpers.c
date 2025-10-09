#include <math.h>
#include <float.h>
#include "utils.h"
#include "helpers.h"

// ---------- Utilidades de vetor/matriz ----------

// Preenche o vetor v com zeros (v[:] = 0)
void vector_fill_zero(real_t *v, int n) {
  for (int i = 0; i < n; ++i) v[i] = 0.0;
}

// Copia src -> dst
void vector_copy(real_t *dst, const real_t *src, int n) {
  for (int i = 0; i < n; ++i) dst[i] = src[i];
}

// Produto interno (dot) entre 'a' e 'b'
real_t dot_product(const real_t *a, const real_t *b, int n) {
  real_t s = 0.0;
  for (int i = 0; i < n; ++i) s += a[i]*b[i];
  return s;
}

// y := y + a * x  (AXPY)
void vector_axpy(real_t *y, real_t a, const real_t *x, int n) {
  for (int i = 0; i < n; ++i) y[i] += a * x[i];
}

// y = A * x  (matriz densa 'A' armazenada linha a linha)
void matvec_multiply_dense(const real_t *A, const real_t *x, real_t *y, int n) {
  for (int i = 0; i < n; ++i) {
    real_t s = 0.0;
    const real_t *Ai = &A[(size_t)i * (size_t)n];
    for (int j = 0; j < n; ++j) s += Ai[j] * x[j];
    y[i] = s;
  }
}

// ||a - b||_inf (máximo |a_i - b_i|)
real_t norm_inf_difference(const real_t *a, const real_t *b, int n) {
  real_t m = 0.0;
  for (int i = 0; i < n; ++i) {
    real_t d = fabs(a[i] - b[i]);
    if (d > m) m = d;
  }
  return m;
}

// ||v||_2 (euclidiana)
real_t norm_l2(const real_t *v, int n) {
  return sqrt(dot_product(v, v, n));
}
