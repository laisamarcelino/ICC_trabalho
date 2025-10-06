#include <stdio.h>
#include <stdlib.h>
#include "sislin.h"
#include "pcgc.h"
#include "utils.h"

int main(void) {
  int n = 5, k = 5, maxit = 1000;
  real_t w = -1.0;      // ω = -1 → sem PC
  real_t eps = 1e-8;

  srandom(20252);

  real_t *A = NULL, *b = NULL;
  criaKDiagonal(n, k, &A, &b);

  real_t *ASP = NULL;
  real_t *bsp = (real_t*) calloc(n, sizeof(real_t));
  rtime_t tpc = 0.0, titer = 0.0, tres = 0.0;
  genSimetricaPositiva(A, b, n, k, &ASP, bsp, &tpc);

  real_t *x = (real_t*) calloc(n, sizeof(real_t));
  real_t res_norm = 0.0;

  int it = cg_no_prec(ASP, bsp, x, n, maxit, eps, &titer, &tres, &res_norm);
  if (it < 0) {
    fprintf(stderr, "Falha no CG sem pré-condicionador.\n");
    return 1;
  }

  // Impressões no formato do enunciado (exemplo mínimo):
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
