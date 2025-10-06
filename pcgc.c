#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "pcgc.h"
#include "utils.h"

// ---------- helpers básicos ----------

static inline void vec_set_zero(real_t *v, int n) {
  for (int i = 0; i < n; ++i) v[i] = 0.0;
}

static inline void vec_copy(real_t *dst, const real_t *src, int n) {
  for (int i = 0; i < n; ++i) dst[i] = src[i];
}

static inline real_t dot(const real_t *a, const real_t *b, int n) {
  real_t s = 0.0;
  for (int i = 0; i < n; ++i) s += a[i]*b[i];
  return s;
}

static inline void axpy(real_t *y, real_t a, const real_t *x, int n) {
  // y := y + a*x
  for (int i = 0; i < n; ++i) y[i] += a * x[i];
}

static inline void matvec_dense(const real_t *A, const real_t *x, real_t *y, int n) {
  for (int i = 0; i < n; ++i) {
    real_t s = 0.0;
    const real_t *Ai = &A[i*(size_t)n];
    for (int j = 0; j < n; ++j) s += Ai[j] * x[j];
    y[i] = s;
  }
}

static inline real_t norm_inf_diff(const real_t *a, const real_t *b, int n) {
  real_t m = 0.0;
  for (int i = 0; i < n; ++i) {
    real_t d = fabs(a[i] - b[i]);
    if (d > m) m = d;
  }
  return m;
}

static inline real_t norm2(const real_t *v, int n) {
  return sqrt(dot(v, v, n));
}

// ---------- CG sem pré-condicionador (ω = -1) ----------

int cg_no_prec(const real_t *A, const real_t *b, real_t *x,
               int n, int maxit, real_t eps,
               rtime_t *t_iter, rtime_t *t_res,
               real_t *res_norm_out)
{
  if (!A || !b || !x || n <= 0 || maxit <= 0 || eps < 0.0) {
    fprintf(stderr, "cg_no_prec: parâmetros inválidos\n");
    return -1;
  }

  // aloca vetores auxiliares
  real_t *r  = (real_t*) malloc((size_t)n * sizeof(real_t));
  real_t *p  = (real_t*) malloc((size_t)n * sizeof(real_t));
  real_t *Ap = (real_t*) malloc((size_t)n * sizeof(real_t));
  real_t *x_old = (real_t*) malloc((size_t)n * sizeof(real_t));
  if (!r || !p || !Ap || !x_old) {
    fprintf(stderr, "cg_no_prec: falha de alocação\n");
    free(r); free(p); free(Ap); free(x_old);
    return -1;
  }

  // x0 = 0
  vec_set_zero(x, n);

  // r0 = b - A*x0 = b
  vec_copy(r, b, n);

  // p0 = r0  (sem pré-condicionador)
  vec_copy(p, r, n);

  real_t rho = dot(r, r, n); // r^T r
  if (rho <= DBL_EPSILON) {
    // já convergiu em x0
    if (t_iter) *t_iter = 0.0;
    // resíduo final (≈ ||b||2, pois x=0)
    if (t_res) {
      *t_res = timestamp();
      *t_res = timestamp() - *t_res;
    }
    if (res_norm_out) *res_norm_out = sqrt(rho);
    free(r); free(p); free(Ap); free(x_old);
    return 0;
  }

  int it = 0;
  rtime_t t_loop = timestamp();

  for (it = 1; it <= maxit; ++it) {
    // Ap = A * p
    matvec_dense(A, p, Ap, n);

    real_t denom = dot(p, Ap, n); // p^T A p
    if (!(denom > 0.0) || !isfinite(denom)) {
      fprintf(stderr, "cg_no_prec: breakdown numérico (p^T A p = %g)\n", denom);
      free(r); free(p); free(Ap); free(x_old);
      return -1;
    }

    real_t alpha = rho / denom;

    // x_{k+1} = x_k + alpha p_k
    vec_copy(x_old, x, n);
    axpy(x, alpha, p, n);

    // r_{k+1} = r_k - alpha Ap
    for (int i = 0; i < n; ++i) r[i] -= alpha * Ap[i];

    // critério do trabalho: ||x_new - x_old||_inf < eps
    real_t aprox_err = norm_inf_diff(x, x_old, n);
    if (aprox_err < eps) break;

    real_t rho_new = dot(r, r, n); // (r_{k+1})^T (r_{k+1})
    if (!isfinite(rho_new)) {
      fprintf(stderr, "cg_no_prec: rho_new inválido\n");
      free(r); free(p); free(Ap); free(x_old);
      return -1;
    }

    real_t beta = rho_new / rho;

    // p_{k+1} = r_{k+1} + beta p_k
    for (int i = 0; i < n; ++i) p[i] = r[i] + beta * p[i];

    rho = rho_new;
  }

  t_loop = timestamp() - t_loop;
  if (t_iter) *t_iter = (it > 0 ? t_loop / (real_t)it : 0.0);

  // resíduo final e tempo_residuo
  if (t_res || res_norm_out) {
    if (t_res) *t_res = timestamp();

    // r = b - A x
    matvec_dense(A, x, Ap, n); // reutiliza Ap como A*x
    for (int i = 0; i < n; ++i) r[i] = b[i] - Ap[i];

    if (res_norm_out) *res_norm_out = norm2(r, n);
    if (t_res) *t_res = timestamp() - *t_res;
  }

  free(r); free(p); free(Ap); free(x_old);
  return it;
}

// ---------- CG com pré-condicionador Jacobi (ω = 0.0) ----------

// Função auxiliar para resolver M z = r, onde M = D (diagonal de A)
// M é o pré-condicionador (M = D)
// z é o resíduo pré-condicionado atual (D⁻¹ × r_atual)
// r é o resíduo atual: r = b - A*x

/*
Exemplo:
M[0,0] * z[0] = r[0]  →  z[0] = r[0] / M[0,0] = r[0] / A[0,0]
M[1,1] * z[1] = r[1]  →  z[1] = r[1] / M[1,1] = r[1] / A[1,1]  
M[2,2] * z[2] = r[2]  →  z[2] = r[2] / M[2,2] = r[2] / A[2,2]
*/
static inline void jacobi_solve(const real_t *A, const real_t *r, real_t *z, int n) {
  for (int i = 0; i < n; ++i) {
    real_t aii = A[i*(size_t)n + i]; // elemento diagonal A[i][i]
    if (fabs(aii) <= DBL_EPSILON) {
      fprintf(stderr, "jacobi_solve: elemento diagonal A[%d][%d] = %g muito pequeno\n", i, i, aii);
      z[i] = 0.0; // fallback
    } else {
      z[i] = r[i] / aii; // resíduo pré-condicionado = resíduo atual / A[i][i]
    }
  }
}

// A: matriz SPD (simétrica positiva definida) n×n linearizada
// b: vetor termo independente de tamanho n
// x: vetor solução (saída) de tamanho n, inicializado em x0 = 0
// n: dimensão do sistema
// maxit: número máximo de iterações
// eps: tolerância para critério de parada ||x_new - x_old||_inf < eps
// t_iter: tempo médio por iteração (saída, em ms)
// t_res: tempo do cálculo do resíduo final (saída, em ms)
// res_norm_out: norma do resíduo final ||b - Ax||_2 (saída, opcional)
int cg_jacobi_prec(const real_t *A, const real_t *b, real_t *x,
                   int n, int maxit, real_t eps,
                   rtime_t *t_iter, rtime_t *t_res,
                   real_t *res_norm_out)
{
  if (!A || !b || !x || n <= 0 || maxit <= 0 || eps < 0.0) {
    fprintf(stderr, "cg_jacobi_prec: parâmetros inválidos\n");
    return -1;
  }

  // aloca vetores auxiliares
  real_t *r  = (real_t*) malloc((size_t)n * sizeof(real_t));
  real_t *z  = (real_t*) malloc((size_t)n * sizeof(real_t));
  real_t *p  = (real_t*) malloc((size_t)n * sizeof(real_t));
  real_t *Ap = (real_t*) malloc((size_t)n * sizeof(real_t));
  real_t *x_old = (real_t*) malloc((size_t)n * sizeof(real_t));
  if (!r || !z || !p || !Ap || !x_old) {
    fprintf(stderr, "cg_jacobi_prec: falha de alocação\n");
    free(r); free(z); free(p); free(Ap); free(x_old);
    return -1;
  }

  // x0 = 0
  vec_set_zero(x, n);

  // r0 = b - A*x0 = b
  vec_copy(r, b, n);

  // z0 = M^{-1} * r0 (resolve M z = r, onde M = D)
  jacobi_solve(A, r, z, n);

  // p0 = z0 (direção inicial)
  vec_copy(p, z, n);

  real_t rho = dot(r, z, n); // r^T z (produto interno modificado)
  if (rho <= DBL_EPSILON) {
    // já convergiu em x0
    if (t_iter) *t_iter = 0.0;
    // resíduo final (≈ ||b||2, pois x=0)
    if (t_res) {
      *t_res = timestamp();
      *t_res = timestamp() - *t_res;
    }
    if (res_norm_out) *res_norm_out = norm2(r, n);
    free(r); free(z); free(p); free(Ap); free(x_old);
    return 0;
  }

  int it = 0;
  rtime_t t_loop = timestamp();

  for (it = 1; it <= maxit; ++it) {
    // Ap = A * p
    matvec_dense(A, p, Ap, n);

    real_t denom = dot(p, Ap, n); // p^T A p
    if (!(denom > 0.0) || !isfinite(denom)) {
      fprintf(stderr, "cg_jacobi_prec: breakdown numérico (p^T A p = %g)\n", denom);
      free(r); free(z); free(p); free(Ap); free(x_old);
      return -1;
    }

    real_t alpha = rho / denom;

    // x_{k+1} = x_k + alpha p_k
    vec_copy(x_old, x, n);
    axpy(x, alpha, p, n);

    // r_{k+1} = r_k - alpha Ap
    for (int i = 0; i < n; ++i) r[i] -= alpha * Ap[i];

    // critério do trabalho: ||x_new - x_old||_inf < eps
    real_t aprox_err = norm_inf_diff(x, x_old, n);
    if (aprox_err < eps) break;

    // z_{k+1} = M^{-1} * r_{k+1}
    jacobi_solve(A, r, z, n);

    real_t rho_new = dot(r, z, n); // r_{k+1}^T z_{k+1}
    if (!isfinite(rho_new)) {
      fprintf(stderr, "cg_jacobi_prec: rho_new inválido\n");
      free(r); free(z); free(p); free(Ap); free(x_old);
      return -1;
    }

    real_t beta = rho_new / rho;

    // p_{k+1} = z_{k+1} + beta p_k
    for (int i = 0; i < n; ++i) p[i] = z[i] + beta * p[i];

    rho = rho_new;
  }

  t_loop = timestamp() - t_loop;
  if (t_iter) *t_iter = (it > 0 ? t_loop / (real_t)it : 0.0);

  // resíduo final e tempo_residuo
  if (t_res || res_norm_out) {
    if (t_res) *t_res = timestamp();

    // r = b - A x
    matvec_dense(A, x, Ap, n); // reutiliza Ap como A*x
    for (int i = 0; i < n; ++i) r[i] = b[i] - Ap[i];

    if (res_norm_out) *res_norm_out = norm2(r, n);
    if (t_res) *t_res = timestamp() - *t_res;
  }

  free(r); free(z); free(p); free(Ap); free(x_old);
  return it;
}
