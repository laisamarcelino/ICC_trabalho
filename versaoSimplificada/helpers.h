#ifndef __HELPERS_H__
#define __HELPERS_H__

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include <stddef.h>
#include "sislin.h"

#define IDX(i, j, n) ((i) * (n) + (j))

/* ---------- Operações vetoriais básicas ---------- */

real_t vet_produto(int n, const real_t *x, const real_t *y);
real_t vet_norma2(int n, const real_t *x);
void   vet_copy(int n, const real_t *fonte, real_t *dst);
void   vet_preenche(int n, real_t val, real_t *x);
void   vet_axpy(int n, real_t alpha, const real_t *x, real_t *y);
void   vet_escala(int n, real_t alpha, real_t *x);

/* ---------- Mat-vec (densa) ---------- */

// Multiplica matriz densa A (n×n) por vetor x: y = A*x
void matvet_densa(const real_t *A, const real_t *x, real_t *y, int n);

/* ---------- Pré-condicionador de Jacobi ---------- */

// Extrai diagonal e seu inverso (checando zeros)
int  extrai_diag_e_invD(const real_t *A, int n, int k, real_t *D, real_t *invD, real_t eps);

// Aplica Jacobi: y = D⁻¹ r
void aplica_jacobi(int n, const real_t *invD, const real_t *r, real_t *y);

/* ---------- Varreduras para SGS/SSOR ---------- */
/*   Mapeamento dos casos:
   - Gauss-Seidel (ω=1): diagScale = 1.0  -> (D + L) e (D + U)
   - SSOR (0<ω<2):      diagScale = 1/ω  -> (D/ω + L) e (D/ω + U)
*/
void forward_sweep_DL(const real_t *A, int n, int k,
                      real_t diagScale,
                      const real_t *rhs, real_t *t);

void backward_sweep_DU(const real_t *A, int n, int k,
                       real_t diagScale,
                       const real_t *rhs, real_t *y);

#endif // __HELPERS_H__
