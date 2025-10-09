#ifndef PCGC_H
#define PCGC_H

#include "utils.h"  // real_t, rtime_t

#ifdef __cplusplus
extern "C" {
#endif

// SEM pré-condicionador (ω = -1)
int cg_solve_no_preconditioner(const real_t *A, const real_t *b, real_t *x,
                               int n, int maxit, real_t eps,
                               rtime_t *t_iter, rtime_t *t_res,
                               real_t *res_norm_l2_out,
                               real_t *approx_err_inf_out);

// Jacobi (ω = 0) — opcional aqui, mantenha o protótipo se for implementar
int cg_solve_jacobi(const real_t *A, const real_t *b, real_t *x,
                    int n, int maxit, real_t eps,
                    rtime_t *t_iter, rtime_t *t_res,
                    real_t *res_norm_l2_out,
                    real_t *approx_err_inf_out);

// Gauss-Seidel (ω = 1) — opcional
int cg_solve_gauss_seidel(const real_t *A, const real_t *b, real_t *x,
                          int n, int maxit, real_t eps,
                          rtime_t *t_iter, rtime_t *t_res,
                          real_t *res_norm_l2_out,
                          real_t *approx_err_inf_out);

// SSOR (1 < ω < 2) — opcional
int cg_solve_ssor(const real_t *A, const real_t *b, real_t *x,
                  int n, int maxit, real_t eps, real_t omega,
                  rtime_t *t_iter, rtime_t *t_res,
                  real_t *res_norm_l2_out,
                  real_t *approx_err_inf_out);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // PCGC_H
