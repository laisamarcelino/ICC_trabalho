#include <stdio.h>
#include <stdlib.h>
#include "sislin.h"
#include "pcgc.h"
#include "utils.h"

int main(void) {
  int n, k, maxit;
  real_t w, eps;

  // Leitura dos parâmetros
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

  // Semente para geração dos dados
  srandom(20252);

  // Geração da matriz k-diagonal e do vetor b (não-SPD)
  real_t *A = NULL, *b = NULL;
  criaKDiagonal(n, k, &A, &b);

  // Transformação p/ SPD (AᵗA ou equivalente) e cópia de b adequada
  real_t *ASP = NULL;
  real_t *bsp = (real_t*) calloc((size_t)n, sizeof(real_t));
  if (!bsp) {
    fprintf(stderr, "Erro de alocação\n");
    free(A); free(b);
    return 1;
  }

  rtime_t tpc = 0.0, titer = 0.0, tres = 0.0;
  genSimetricaPositiva(A, b, n, k, &ASP, bsp, &tpc);

  // Vetor solução
  real_t *x = (real_t*) calloc((size_t)n, sizeof(real_t));
  if (!x) {
    fprintf(stderr, "Erro de alocação\n");
    free(A); free(b); free(ASP); free(bsp);
    return 1;
  }

  // Saídas do solver
  real_t res_norm_l2 = 0.0;   // ||r||_2 ao final
  real_t approx_err_inf = 0.0; // ||x_k - x_{k-1}||_inf da última iteração
  int it = -1;

  // Seleção do método conforme ω
  if (w == -1.0) {
    it = cg_solve_no_preconditioner(ASP, bsp, x, n, maxit, eps,
                                    &titer, &tres,
                                    &res_norm_l2, &approx_err_inf);
  } else if (w == 0.0) {
    it = cg_solve_jacobi(ASP, bsp, x, n, maxit, eps,
                         &titer, &tres,
                         &res_norm_l2, &approx_err_inf);
  } else if (w == 1.0) {
    it = cg_solve_gauss_seidel(ASP, bsp, x, n, maxit, eps,
                               &titer, &tres,
                               &res_norm_l2, &approx_err_inf);
  } else if (w > 1.0 && w < 2.0) {
    it = cg_solve_ssor(ASP, bsp, x, n, maxit, eps, w,
                       &titer, &tres,
                       &res_norm_l2, &approx_err_inf);
  } else {
    fprintf(stderr, "Erro: Valor de ω inválido (%.2f). Valores aceitos:\n", w);
    fprintf(stderr, "  ω = -1.0: sem pré-condicionador\n");
    fprintf(stderr, "  ω = 0.0:  Jacobi\n");
    fprintf(stderr, "  ω = 1.0:  Gauss-Seidel\n");
    fprintf(stderr, "  1.0 < ω < 2.0: SSOR\n");
    free(A); free(b); free(ASP); free(bsp); free(x);
    return 1;
  }

  if (it < 0) {
    fprintf(stderr, "Falha no método de Gradientes Conjugados.\n");
    free(A); free(b); free(ASP); free(bsp); free(x);
    return 1;
  }

  // ---------- Saída no formato do enunciado ----------
  // 1) n
  printf("%d\n", n);

  // 2) vetor solução x (n valores, %.16g, separados por espaço)
  for (int i = 0; i < n; ++i) {
    printf("%.16g%c", x[i], (i + 1 < n ? ' ' : '\n'));
  }

  // 3) norma: ||x_k - x_{k-1}||_inf (critério de parada)
  printf("%.8g\n", approx_err_inf);

  // 4) resíduo: ||r||_2 ao final (r = b - A x)
  printf("%.8g\n", res_norm_l2);

  // 5) tempos (ms): pc, iteração média, resíduo final
  printf("%.8g\n", tpc);
  printf("%.8g\n", titer);
  printf("%.8g\n", tres);

  // limpeza
  free(A); free(b); free(ASP); free(bsp); free(x);
  return 0;
}
