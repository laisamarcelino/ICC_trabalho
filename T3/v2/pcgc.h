#ifndef __PCG_H__
#define __PCG_H__

#include <stddef.h>
#include "sislin.h"
#include "helpers.h"
#include "utils.h"

/* PCG -> Gradientes Conjugados Pré-Condicionados
 * Este header define:
 *  - o enum com os tipos de pré-condicionador (NONE/JACOBI/SGS/SSOR)
 *  - o “contexto” do PC (buffers e parâmetros)
 *  - as funções de setup/apply/free do PC
 *  - o solver unificado (CG puro + PCG) com ω (para GS/SSOR)
 */

typedef enum {
  PCG_PRECOND_NONE = 0,
  PCG_PRECOND_JACOBI, // Jacobi (ω=0)
  PCG_PRECOND_SGS,   // Gauss-Seidel (ω=1) 
  PCG_PRECOND_SSOR   // SSOR (0<ω<2) 
} pcg_precond_t;      // Tipo do pré-condicionador 

typedef struct {
  pcg_precond_t type;  // qual PC foi escolhido
  real_t omega;        // GS=1.0; SSOR=ω (0<ω<2)

  // Buffers de trabalho compartilhados
  real_t *y;           // y = M⁻¹ r (sempre que M != NONE)

  // Jacobi
  real_t *D, *invD;    // diagonal e seu inverso

  // GS/SSOR
  real_t *t, *u;       // workspaces para as duas varreduras (progressiva e regresiva)
} pcg_contexto_t;

// Setup/Apply/Free do pré-condicionador
int pcg_setup(const matdiag_t *A, int n, int k, pcg_precond_t M, real_t omega, pcg_contexto_t *contexto);
void pcg_apply(const matdiag_t *A, int n, int k, const pcg_contexto_t *contexto, const real_t *r, real_t *y);
void pcg_free(pcg_contexto_t *contexto);

// Solver unificado (CG puro + PCG)
int cg_solve(const matdiag_t *A, const real_t *b, real_t *x,
             int n, int k, int maxit, real_t eps_inf,
             pcg_precond_t M, real_t omega,
             real_t *norma_delta_x_inf_out);

#endif // __PCG_H__
