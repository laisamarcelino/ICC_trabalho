#ifndef __PCGC_H__
#define __PCGC_H__

#include "utils.h"   // real_t, rtime_t
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

// CG sem pré-condicionador (ω = -1 → M = I)
// Entradas:
//   A  : matriz SPD (linearizada, n×n) — use a ASP = AᵗA
//   b  : vetor termo independente — use bsp = Aᵗb
//   x  : solução (saída); será escrito a partir de x0 = 0
//   n  : dimensão
//   maxit : máx. iterações
//   eps: tolerância no critério ||x^k - x^{k-1}||_inf < eps
// Saídas medição:
//   t_iter: tempo MÉDIO por iteração (ms)
//   t_res : tempo do cálculo do resíduo final (ms)
//   res_norm_out: ||b - A x_final||_2 (opcional; pode ser NULL)
// Retorno: número de iterações realizadas (>0), ou
//          0 se já convergiu em x0, ou
//         -1 em quebra numérica / denom <= 0 / NaN etc.
int cg_no_prec(const real_t *A, const real_t *b, real_t *x,
               int n, int maxit, real_t eps,
               rtime_t *t_iter, rtime_t *t_res,
               real_t *res_norm_out);

// CG com pré-condicionador de Jacobi (ω = 0.0 → M = D)
// Onde D é a matriz diagonal de A
int cg_jacobi_prec(const real_t *A, const real_t *b, real_t *x,
                   int n, int maxit, real_t eps,
                   rtime_t *t_iter, rtime_t *t_res,
                   real_t *res_norm_out);

#ifdef __cplusplus
}
#endif
#endif // __PCGC_H__
