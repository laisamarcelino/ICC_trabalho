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
void vet_copy(int n, const real_t *fonte, real_t *dst);
void vet_preenche(int n, real_t val, real_t *x);
void vet_axpy(int n, real_t alpha, const real_t *x, real_t *y);
void vet_escala(int n, real_t alpha, real_t *x);

/* ---------- Mat-vec (diagonal) ---------- */

void matvet_diagonais(const matdiag_t *A, const real_t *x, real_t *y);

/* ---------- Pré-condicionador de Jacobi ---------- */

// Extrai diagonal e seu inverso da diagonal
int extrai_diag_e_invD_diag(const matdiag_t *A, real_t *D, real_t *invD, real_t eps);

// Aplica Jacobi: y = D⁻¹ r
void aplica_jacobi(int n, const real_t *invD, const real_t *r, real_t *y);

/* ---------- Varreduras para SGS/SSOR ---------- *
   Mapeamento dos casos:
   - Gauss-Seidel (ω=1): escala_diag = 1.0  -> (D + L) e (D + U)
   - SSOR (0<ω<2):      escala_diag = 1/ω  -> (D/ω + L) e (D/ω + U)
*/
void varredura_progressiva_DL(const matdiag_t *A, int n, int k,
                              real_t escala_diag, const real_t *vet_dir, real_t *y);

void varredura_regressiva_DU(const matdiag_t *A, int n, int k,
                             real_t escala_diag, const real_t *vet_dir, real_t *y);

/* ---------- Cálculo do resíduo em norma L2 ---------- */
real_t residuo_l2_v2(const matdiag_t *A, const real_t *b, const real_t *x);

void liberaMatDiag(matdiag_t *A);

#endif // __HELPERS_H__
