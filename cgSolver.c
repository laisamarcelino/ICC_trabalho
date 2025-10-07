#include <stdio.h>
#include <stdlib.h>
#include "sislin.h"
#include "pcgc.h"
#include "utils.h"

int main(void) {
  int n, k, maxit;
  real_t w, eps;

  // Lê os parâmetros da entrada padrão
  if (scanf("%d", &n) != 1 || n <= 10) {
    fprintf(stderr, "Erro: n deve ser > 10\n");
    return 1;
  }
  
  if (scanf("%d", &k) != 1 || k <= 1 || k % 2 == 0) {
    fprintf(stderr, "Erro: k deve ser > 1 e ímpar\n");
    return 1;
  }
  
  if (scanf("%lf", &w) != 1) {
    fprintf(stderr, "Erro ao ler ω\n");
    return 1;
  }
  
  if (scanf("%d", &maxit) != 1 || maxit <= 0) {
    fprintf(stderr, "Erro: maxit deve ser > 0\n");
    return 1;
  }
  
  if (scanf("%lf", &eps) != 1 || eps <= 0) {
    fprintf(stderr, "Erro: ε deve ser > 0\n");
    return 1;
  }

  srandom(20252);

  real_t *A = NULL, *b = NULL;
  // Gerando matriz k-diagonal
  criaKDiagonal(n, k, &A, &b);

  real_t *ASP = NULL;
  real_t *bsp = (real_t*) calloc(n, sizeof(real_t));
  rtime_t tpc = 0.0, titer = 0.0, tres = 0.0;
  // Gerando A = Aᵗ positiva
  genSimetricaPositiva(A, b, n, k, &ASP, bsp, &tpc);

  real_t *x = (real_t*) calloc(n, sizeof(real_t));
  real_t res_norm = 0.0;
  int it = -1;

  // Seleção do método baseado no valor de ω

  // ω = -1: sem pré-condicionador
  if (w == -1.0) {
    it = cg_no_prec(ASP, bsp, x, n, maxit, eps, &titer, &tres, &res_norm);
  }
  
  // ω = 0.0: pré-condicionador de Jacobi
  else if (w == 0.0) {
    it = cg_jacobi_prec(ASP, bsp, x, n, maxit, eps, &titer, &tres, &res_norm);
  }

  // ω = 1.0: pré-condicionador de Gauss-Seidel
  else if (w == 1.0) {
    // TODO: Implementar CG com pré-condicionador de Gauss-Seidel
    // it = cg_gauss_seidel_prec(ASP, bsp, x, n, maxit, eps, &titer, &tres, &res_norm);
  }

  // 1.0 < ω < 2.0: pré-condicionador SSOR
  else if (w > 1.0 && w < 2.0) {
    // TODO: Implementar CG com pré-condicionador SSOR
    // it = cg_ssor_prec(ASP, bsp, x, n, maxit, eps, w, &titer, &tres, &res_norm);
  }

  else {
    // Valor de ω inválido
    fprintf(stderr, "Erro: Valor de ω inválido (%.2f). Valores aceitos:\n", w);
    fprintf(stderr, "  ω = -1.0: sem pré-condicionador\n");
    fprintf(stderr, "  ω = 0.0:  pré-condicionador de Jacobi\n");
    fprintf(stderr, "  ω = 1.0:  pré-condicionador de Gauss-Seidel\n");
    fprintf(stderr, "  1.0 < ω < 2.0: pré-condicionador SSOR\n");
    free(A); free(b); free(ASP); free(bsp); free(x);
    return 1;
  }

  if (it < 0) {
    fprintf(stderr, "Falha no CG sem pré-condicionador.\n");
    return 1;
  }

  // Impressões no formato do enunciado
  printf("%d\n", n);
  for (int i = 0; i < n; ++i) {
    printf("%.16g%c", x[i], (i+1<n ? ' ' : '\n'));
  }
  // norma aprox ||x^k - x^{k-1}||_inf foi o critério; aqui só mostramos o resíduo final:
  printf("%.8g\n", res_norm);      // residuo (||r||2)
  printf("%.8g\n", res_norm);      // (se quiser separar, calcule/guarde norma aprox)
  printf("%.8g\n", tpc);           // tempo_pc
  printf("%.8g\n", titer);         // tempo_iter (médio)
  printf("%.8g\n", tres);          // tempo_residuo

  free(A); free(b); free(ASP); free(bsp); free(x);
  return 0;
}