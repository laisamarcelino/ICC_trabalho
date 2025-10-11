/* main.c — testes com CG (NONE) e CG com Jacobi
 *
 * Compilar (exemplo):
 *   gcc -O2 -std=c11 -o cgSolver main.c pcgc.c helpers.c sislin.c utils.c -lm
 */

 #include <stdio.h>
 #include <stdlib.h>
 #include <string.h>
 #include <math.h>
 
 #include "pcgc.h"     /* pcg_precond_t e assinatura nova de cg_solve_none */
 #include "helpers.h"
 #include "sislin.h"
 
 /* assinatura da SUA função (com param de pré-condicionador) */
 int cg_solve_none(const real_t *A, const real_t *b, real_t *x,
                   int n, int k, int maxit, real_t eps_inf,
                   pcg_precond_t M);
 
 static void imprime_vetor(const char *nome, const real_t *v, int n) {
     printf("%s = [", nome);
     for (int i = 0; i < n; ++i) {
         printf("%.16g", v[i]);
         if (i < n-1) printf(", ");
     }
     printf("]\n");
 }
 
 static real_t residual_l2(const real_t *A, int n, int k,
                           const real_t *x, const real_t *b)
 {
     real_t *Ax = (real_t*) calloc((size_t)n, sizeof(real_t));
     real_t *r  = (real_t*) calloc((size_t)n, sizeof(real_t));
     band_matvet(A, n, k, x, Ax);           /* Ax = A*x */
     for (int i = 0; i < n; ++i) r[i] = b[i] - Ax[i];
     real_t rn2 = vet_norma2(n, r);
     free(Ax); free(r);
     return rn2;
 }
 
 int main(void)
 {
     /* tamanho do sistema e largura de banda */
     const int n = 8;
     const int k = 3;   /* tri-diagonal → p=1, k=2*p+1=3 */
 
     /* alocações */
     real_t *A = (real_t*) calloc((size_t)n*(size_t)n, sizeof(real_t));
     real_t *b = (real_t*) calloc((size_t)n, sizeof(real_t));
     real_t *x_none   = (real_t*) calloc((size_t)n, sizeof(real_t));
     real_t *x_jacobi = (real_t*) calloc((size_t)n, sizeof(real_t));
     if (!A || !b || !x_none || !x_jacobi) {
         fprintf(stderr, "Falha de alocacao\n");
         return 1;
     }
 
     /* === monta uma matriz SPD simples tri-diagonal ===
        A = 4*I + 1*sub/super-diagonais
     */
     for (int i = 0; i < n; ++i) {
         A[IDX(i,i,n)] = 4.0;
         if (i-1 >= 0) A[IDX(i,i-1,n)] = 1.0;
         if (i+1 <  n) A[IDX(i,i+1,n)] = 1.0;
     }
 
     /* b = [1,1,...,1] */
     for (int i = 0; i < n; ++i) b[i] = 1.0;
 
     /* parâmetros do CG: critério do TRABALHO usa eps_inf em ||Δx||_∞ */
     const int   maxit   = 1000;
     const real_t epsinf = 1e-12;
 
     /* ========================= CG (NONE) ========================= */
     int it_none = cg_solve_none(A, b, x_none, n, k, maxit, epsinf, PCG_PRECOND_NONE);
     real_t r2_none = residual_l2(A, n, k, x_none, b);
 
     printf("=== CG (NONE) ===\n");
     printf("Iteracoes: %d\n", it_none);
     printf("||r||_2   = %.8g\n", (double)r2_none);
     imprime_vetor("x_none", x_none, n);
     puts("");
 
     /* ====================== CG + Jacobi ========================== */
     int it_jac = cg_solve_none(A, b, x_jacobi, n, k, maxit, epsinf, PCG_PRECOND_JACOBI);
     real_t r2_jac = residual_l2(A, n, k, x_jacobi, b);
 
     printf("=== CG (Jacobi) ===\n");
     printf("Iteracoes: %d\n", it_jac);
     printf("||r||_2   = %.8g\n", (double)r2_jac);
     imprime_vetor("x_jacobi", x_jacobi, n);
 
     /* limpeza */
     free(A); free(b); free(x_none); free(x_jacobi);
     return 0;
 }
 