// Versão compacta funcional até agora

#ifndef PCG_H
#define PCG_H

#include <stddef.h>
#include "sislin.h"
#include "helpers.h"
#include "utils.h"

typedef enum {
  PCG_PRECOND_NONE = 0,
  PCG_PRECOND_JACOBI,
  PCG_PRECOND_SGS,
  PCG_PRECOND_SSOR
} pcg_precond_t;

/* Solver unificado: CG “puro” (NONE) ou PCG Jacobi (M=JACOBI).
   Critério de parada: ||Δx||_inf < eps_inf. */
int cg_solve_none(const real_t *A, const real_t *b, real_t *x,
                  int n, int k, int maxit, real_t eps_inf,
                  pcg_precond_t M);

#endif /* PCG_H */
