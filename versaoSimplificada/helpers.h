#ifndef HELPERS_H
#define HELPERS_H

#include "sislin.h"
#include <stddef.h>

#define IDX(i, j, n) ((i) * (n) + (j))

/* ---------- Operações vetoriais básicas ---------- */

real_t vet_produto(int n, const real_t *x, const real_t *y);
real_t vet_norma2(int n, const real_t *x);
void vet_copy(int n, const real_t *fonte, real_t *dst);
void vet_preenche(int n, real_t val, real_t *x);
void vet_axpy(int n, real_t alpha, const real_t *x, real_t *y);
void vet_escala(int n, real_t alpha, real_t *x);


/* ---------- Mat-vec otimizado para k-diagonal ---------- */

/* z = A x, explorando a banda k (k = 2*p+1). */
void band_matvet(const real_t *A, int n, int k, const real_t *x, real_t *z);


/* ---------- Preparação e aplicação do pré-condicionador de Jacobi ---------- */

/* Extrai diagonal e seu inverso (checando zeros).
   Retorna 0 se OK; !=0 se encontrou diagonal ~0 (para você tratar no Jacobi/SSOR). */
int extrai_diag_e_invD(const real_t *A, int n, int k, real_t *D, 
            real_t *invD, real_t eps); /* TODO */

/* Aplica o Pré-condicionador de Jacobi y = M^{-1} r = D^{-1} r */
void aplica_jacobi(int n, const real_t *invD, const real_t *r, real_t *y);

/* (Opcional) SETUP “empacotado” do Jacobi: aloca D e invD e chama extract_
   Retorna 0 se OK; !=0 se falha (ex.: diagonal ~0 ou alocação). */
int setup_jacobi(const real_t *A, int n, int k,
                 real_t **D_out, real_t **invD_out, real_t eps);

                
/* ---------- Varreduras para SGS/SSOR (sem formar L/U) ---------- */
/* Resolve (D + L) t = rhs (forward sweep), usando apenas elementos da banda.
   Use diagScale=1.0 para SGS; para SSOR use diagScale = 1.0/omega ao montar (D/ω + L).
*/
void forward_sweep_DL(const real_t *A, int n, int k,
                      real_t diagScale,
                      const real_t *rhs, real_t *t); /* TODO */

/* Resolve (D + U) y = rhs (backward sweep).
   Use diagScale=1.0 para SGS; para SSOR use diagScale = 1.0/omega para (D/ω + U).
*/
void backward_sweep_DU(const real_t *A, int n, int k,
                       real_t diagScale,
                       const real_t *rhs, real_t *y); /* TODO */

#endif /* HELPERS_H */