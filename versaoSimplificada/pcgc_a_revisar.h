#ifndef PCG_H
#define PCG_H

#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#include "helpers.h"
#include "sislin.h"   
#include "utils.h"    

/* ---------------- Tipos e parâmetros ---------------- */

/* Pré-condicionadores suportados.
   ATENÇÃO (CG): usar apenas PCs SPD. 
   - PCG_PRECOND_SGS equivale a SSOR com omega = 1.0 (simétrico).
*/
typedef enum {
  PCG_PRECOND_NONE = 0,
  PCG_PRECOND_JACOBI,
  PCG_PRECOND_SGS,     /* Gauss-Seidel simétrico (forward+backward) */
  PCG_PRECOND_SSOR     /* SSOR(omega) */
} pcg_precond_t;

/* Parâmetros do solver */
typedef struct {
  int maxit;          /* máximo de iterações */
  real_t tol;         /* tolerância (use ||r||2 ou r^T M^-1 r coerente com as métricas) */
  pcg_precond_t M;    /* tipo de pré-condicionador */
  real_t omega;       /* parâmetro do SSOR (0<ω<2). Ignorado nos demais. */
} pcg_params_t;

/* Métricas para relatório de eficiência */
typedef struct {
  int iters;          /* iterações realizadas */
  rtime_t t_setup;    /* tempo para preparar o pré-condicionador */
  rtime_t t_iter;     /* soma do tempo dentro do loop (pode dividir por iters para média) */
  rtime_t t_total;    /* tempo total (setup + loop) */
  real_t norm_r2;     /* ||r||_2 ao final */
  real_t rtMr;        /* r^T M^-1 r ao final (se aplicável) */
  real_t err_inf;     /* ||Δx||_∞, se você calcular fora/ao final */
} pcg_metrics_t;

/* Dados internos do pré-condicionador (alocados no setup e usados no loop) */
typedef struct {
  pcg_precond_t type;
  real_t omega;         /* para SSOR */
  /* Jacobi */
  real_t *D;            /* diagonal de A (opcional; útil para debug) */
  real_t *invD;         /* 1/diag(A) — uso direto no apply */
  /* Para SGS/SSOR: não guardamos L/U explicitamente; 
     as varreduras usam A e k (largura da banda). 
     Se preferir, você pode preencher ponteiros para L/U mais tarde. */
} pcg_prec_data_t;

/* ------------------- Setup / Apply ------------------- */

/* Prepara dados do pré-condicionador.
   - A: matriz cheia n×n, mas com banda k-diagonal já gerada pelo professor
   - n, k: ordem e nº de diagonais (k = 2*p + 1)
   Saída:
   - pd: estrutura preenchida (aloca invD quando preciso)
   - retorna tempo de setup (em ms) também por conveniência.
*/
rtime_t pcg_setup_precond(const real_t *A, int n, int k,
                          pcg_precond_t type, real_t omega,
                          pcg_prec_data_t *pd);  /* TODO: implementar em pcg.c */

/* y = M^{-1} r  (aplicação do pré-condicionador).
   Implementações esperadas:
   - NONE: y = r
   - JACOBI: y[i] = invD[i] * r[i]
   - SGS: forward ((D+L) t = r), scale (u = D^{-1} t), backward ((D+U) y = u)
   - SSOR: backward ((D/ω + U) t = r), scale (u = D^{-1} t),
           forward ((D/ω + L) y = u), y *= ((2-ω)/ω)
*/
void pcg_apply_precond(const real_t *A, int n, int k,
                       const pcg_prec_data_t *pd,
                       const real_t *r, real_t *y);   /* TODO: implementar em pcg.c */

/* -------------------- Solver PCG --------------------- */

/* Resolve A x = b pelo Gradiente Conjugado (único loop) com pré-condicionador M.
   - A: n×n (k-diagonal)
   - b, x: vetores (x pode chegar zerado)
   - params: escolher M e ω
   - métricas preenchidas para seu relatório
   Retorno:
   - 0 em caso de sucesso; !=0 em erro (ex.: diagonal nula para Jacobi)
*/
int pcg_solve(const real_t *A, const real_t *b, real_t *x,
              int n, int k,
              const pcg_params_t *params,
              pcg_metrics_t *out);                 /* TODO: implementar em pcg.c */

/* Libera buffers do pré-condicionador (se alocados no setup) */
void pcg_free_precond(pcg_prec_data_t *pd);        /* TODO: implementar em pcg.c */


#endif /* PCG_H */
