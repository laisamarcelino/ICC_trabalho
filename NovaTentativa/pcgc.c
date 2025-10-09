#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "pcgc.h"
#include "utils.h"
#include "helpers.h"

// ---------- CG sem pré-condicionador (ω = -1) ----------
// Resolve Ax=b com parada em ||x_{k+1}-x_k||_inf < eps.
// Mede: tempo médio por iteração (t_iter) e tempo p/ resíduo final (t_res).
int cg_solve_no_preconditioner(const real_t *A, const real_t *b, real_t *x,
                               int n, int maxit, real_t eps,
                               rtime_t *t_iter, rtime_t *t_res,
                               real_t *res_norm_l2_out,
                               real_t *approx_err_inf_out)
{
  if (!A || !b || !x || n <= 0 || maxit <= 0 || eps < 0.0) {
    fprintf(stderr, "cg_solve_no_preconditioner: parâmetros inválidos\n");
    return -1;
  }

  real_t *r = (real_t*) malloc((size_t)n * sizeof(real_t));
  real_t *p = (real_t*) malloc((size_t)n * sizeof(real_t));
  real_t *Ap = (real_t*) malloc((size_t)n * sizeof(real_t));
  real_t *x_old = (real_t*) malloc((size_t)n * sizeof(real_t));

  if (!r || !p || !Ap || !x_old) {
    fprintf(stderr, "cg_solve_no_preconditioner: falha de alocação\n");
    free(r); free(p); free(Ap); free(x_old);
    return -1;
  }

  // [CG-1] x0 = 0
  vector_fill_zero(x, n);

  // [CG-2] r0 = b - A x0  (x0=0 => r0=b)
  vector_copy(r, b, n);

  // [CG-3] p0 = r0
  vector_copy(p, r, n);

  // rho_k = r_k^T r_k
  real_t rho = dot_product(r, r, n);

  // Valor padrão para norma aprox (caso trivial)
  if (approx_err_inf_out) *approx_err_inf_out = 0.0;

  // Caso trivial: já satisfeito Ax=b
  if (rho <= DBL_EPSILON) {
    if (t_iter) *t_iter = 0.0;

    // Tempo do resíduo final (aqui r=b)
    if (t_res || res_norm_l2_out) {
      rtime_t t0 = t_res ? timestamp() : 0.0;
      if (res_norm_l2_out) *res_norm_l2_out = sqrt(rho); // ||r||_2
      if (t_res) *t_res = timestamp() - t0;
    }

    free(r); free(p); free(Ap); free(x_old);
    return 0;
  }

  int it = 0;
  rtime_t t_loop0 = timestamp();

  real_t last_approx_err = INFINITY; // manter última ||Δx||_inf

  // [CG-4] Laço principal
  for (it = 1; it <= maxit; ++it) {
    // [CG-5] Ap_k = A * p_k
    matvec_multiply_dense(A, p, Ap, n);

    // denom = p_k^T A p_k
    real_t denom = dot_product(p, Ap, n);
    if (!(denom > 0.0) || !isfinite(denom)) {
      fprintf(stderr, "cg_solve_no_preconditioner: breakdown numérico (p^T A p = %g)\n", denom);
      free(r); free(p); free(Ap); free(x_old);
      return -1;
    }

    // [CG-6] alpha_k
    real_t alpha = rho / denom;

    // [CG-7] x_{k+1} = x_k + alpha p_k
    vector_copy(x_old, x, n);
    vector_axpy(x, alpha, p, n);

    // [CG-8] r_{k+1} = r_k - alpha Ap_k
    for (int i = 0; i < n; ++i) r[i] -= alpha * Ap[i];

    // [CG-9] critério: ||x_{k+1}-x_k||_inf < eps
    real_t aprox_err = norm_inf_difference(x, x_old, n);
    last_approx_err = aprox_err;
    if (aprox_err < eps) break;

    // rho_{new} = r_{k+1}^T r_{k+1}
    real_t rho_new = dot_product(r, r, n);
    if (!isfinite(rho_new)) {
      fprintf(stderr, "cg_solve_no_preconditioner: rho_new inválido\n");
      free(r); free(p); free(Ap); free(x_old);
      return -1;
    }

    // [CG-10] beta_k
    real_t beta = rho_new / rho;

    // [CG-11] p_{k+1} = r_{k+1} + beta p_k
    for (int i = 0; i < n; ++i) p[i] = r[i] + beta * p[i];

    rho = rho_new;
  }

  // Tempo médio por iteração
  rtime_t t_loop = timestamp() - t_loop0;
  if (t_iter) *t_iter = (it > 0) ? (t_loop / (real_t)it) : 0.0;

  // [CG-12] r = b - A x; ||r||_2 e tempo p/ computá-lo
  if (t_res || res_norm_l2_out) {
    rtime_t t0 = t_res ? timestamp() : 0.0;

    matvec_multiply_dense(A, x, Ap, n);
    for (int i = 0; i < n; ++i) r[i] = b[i] - Ap[i];

    if (res_norm_l2_out) *res_norm_l2_out = norm_l2(r, n);
    if (t_res) *t_res = timestamp() - t0;
  }

  if (approx_err_inf_out) *approx_err_inf_out = last_approx_err;

  free(r); free(p); free(Ap); free(x_old);
  return it; // nº de iterações realizadas
}
