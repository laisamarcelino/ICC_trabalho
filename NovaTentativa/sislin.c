#include <stdio.h>
#include <stdlib.h>   // exit, random/srandom
#include <string.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"

// Macro de indexação (linha-major): A[i][j] == A[ i*n + j ]
#define IDX(i,j,n) ((i)*(n) + (j))

// ---------------------------------------------------------------------
// Geração de entradas aleatórias (controladas por srandom no main)
// ---------------------------------------------------------------------

/**
 * Coeficiente A_{ij} para matriz k-diagonal:
 * - Diagonal principal recebe "peso" maior (≈ 2k) para ajudar diagonal-dominância.
 * - Fora da diagonal (na banda), coeficiente ~ U[0,1].
 */
static inline real_t generateRandomA(unsigned int i, unsigned int j, unsigned int k) {
  static const real_t invRandMax = 1.0 / (real_t)RAND_MAX;
  const real_t base = (i == j) ? (real_t)(k << 1) : 1.0;   // 2k ou 1
  return base * (real_t)random() * invRandMax;
}

/**
 * Termo independente b_i: escala ~ 4k * U[0,1].
 */
static inline real_t generateRandomB(unsigned int k) {
  static const real_t invRandMax = 1.0 / (real_t)RAND_MAX;
  return (real_t)(k << 2) * (real_t)random() * invRandMax; // 4k * U[0,1]
}

// ---------------------------------------------------------------------
// Criação de sistema k-diagonal denso (armazenado n x n por linhas)
// ---------------------------------------------------------------------

void criaKDiagonal(int n, int k, real_t **A, real_t **B) {
  *A = (real_t*) calloc((size_t)n * (size_t)n, sizeof(real_t));
  *B = (real_t*) calloc((size_t)n, sizeof(real_t));
  if (!*A || !*B) {
    fprintf(stderr, "Erro ao alocar memória em criaKDiagonal\n");
    free(*A); free(*B);
    exit(2);
  }

  const int p = k / 2; // nº de diagonais abaixo/acima da principal

  for (int i = 0; i < n; ++i) {
    // faixa da banda na linha i
    const int j0 = (i - p < 0) ? 0 : (i - p);
    const int j1 = (i + p >= n) ? (n - 1) : (i + p);

    for (int j = j0; j <= j1; ++j) {
      (*A)[IDX(i, j, n)] = generateRandomA((unsigned)i, (unsigned)j, (unsigned)k);
    }

    // termo independente
    (*B)[i] = generateRandomB((unsigned)k);
  }
}

// ---------------------------------------------------------------------
// Transformação para sistema SPD equivalente: A_spd = A^T A, b_spd = A^T b
// Contabiliza tempo de "pré-condicionamento" pedido na especificação.
// Custo: O(n^3) (denso). Para k-diagonal esparso idealmente faria A^T A esparso.
// ---------------------------------------------------------------------

void genSimetricaPositiva(real_t *A, real_t *b, int n, int k,
                          real_t **ASP, real_t *bsp, rtime_t *tempo)
{
  (void)k; // não é necessário aqui, mantido para assinatura pedida
  if (tempo) *tempo = timestamp();

  // Aloca A_spd = A^T A
  *ASP = (real_t*) calloc((size_t)n * (size_t)n, sizeof(real_t));
  if (!*ASP) {
    fprintf(stderr, "Erro ao alocar memória para ASP\n");
    exit(1);
  }

  // zera b_spd
  for (int i = 0; i < n; ++i) bsp[i] = 0.0;

  // A_spd = A^T A
  // ASP[i,j] = sum_l A[l,i] * A[l,j]
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      real_t sum = 0.0;
      for (int l = 0; l < n; ++l) {
        const real_t a_li = A[l * n + i]; // A[l][i]
        const real_t a_lj = A[l * n + j]; // A[l][j]
        sum += a_li * a_lj;
      }
      (*ASP)[(size_t)i * (size_t)n + (size_t)j] = sum;
    }
  }

  // b_spd = A^T b
  // bsp[i] = sum_l A[l,i] * b[l]
  for (int i = 0; i < n; ++i) {
    real_t sum = 0.0;
    for (int l = 0; l < n; ++l) {
      sum += A[l * n + i] * b[l];
    }
    bsp[i] = sum;
  }

  if (tempo) *tempo = timestamp() - *tempo;
}

// ---------------------------------------------------------------------
// Stubs para D, L, U e pré-condicionador (opcionais/bonus no trabalho).
// Mantidos para extensão futura (IC(0), SSOR etc.). Não usados no fluxo base.
// ---------------------------------------------------------------------

void geraDLU(real_t *A, int n, int k,
             real_t **D, real_t **L, real_t **U, rtime_t *tempo)
{
  if (tempo) *tempo = timestamp();
  (void)A; (void)n; (void)k;

  // Alocar matrizes n x n preenchidas com zero
  *D = (real_t*) calloc((size_t)n * (size_t)n, sizeof(real_t));
  *L = (real_t*) calloc((size_t)n * (size_t)n, sizeof(real_t));
  *U = (real_t*) calloc((size_t)n * (size_t)n, sizeof(real_t));
  if (!*D || !*L || !*U) {
    fprintf(stderr, "geraDLU: erro de alocação\n");
    free(*D); free(*L); free(*U);
    exit(1);
  }

  // Observação: implementação real deveria decompor A em D (diag), L (estritamente inf), U (estritamente sup).
  // Aqui deixamos zerado para indicar "não implementado" sem quebrar quem chama.
  if (tempo) *tempo = timestamp() - *tempo;
}

/**
 * Gera (ou prepara) M^{-1} para o pré-condicionador.
 * Para Jacobi: M = D  => M^{-1} é só o inverso da diagonal (armazenável em vetor).
 * Aqui mantemos stub genérico. Não utilizado no fluxo base sem bônus.
 */
void geraPreCond(real_t *D, real_t *L, real_t *U, real_t w, int n, int k,
                 real_t **M, rtime_t *tempo)
{
  if (tempo) *tempo = timestamp();
  (void)L; (void)U; (void)w; (void)n; (void)k;

  // Exemplo de Jacobi: M seria apenas Diagonal (copiar D).
  // Aqui, apenas duplicamos D como matriz cheia (stub).
  *M = (real_t*) calloc((size_t)n * (size_t)n, sizeof(real_t));
  if (!*M) {
    fprintf(stderr, "geraPreCond: erro de alocação\n");
    exit(1);
  }
  // Copia D (assumindo D é matriz cheia com apenas diagonal populada)
  memcpy(*M, D, (size_t)n * (size_t)n * sizeof(real_t));

  if (tempo) *tempo = timestamp() - *tempo;
}

// ---------------------------------------------------------------------
// Norma L2 do resíduo ||r||_2, r = b - A x. Mede tempo do cálculo.
// Útil para checagens e/ou impressão final.
// ---------------------------------------------------------------------

real_t calcResiduoSL(real_t *A, real_t *b, real_t *X,
                     int n, int k, rtime_t *tempo)
{
  (void)k; // não usado aqui, matriz é densa
  if (tempo) *tempo = timestamp();

  real_t *y = (real_t*) calloc((size_t)n, sizeof(real_t)); // y = A x
  real_t *r = (real_t*) calloc((size_t)n, sizeof(real_t)); // r = b - y
  if (!y || !r) {
    fprintf(stderr, "calcResiduoSL: erro de alocação\n");
    free(y); free(r);
    return -1.0;
  }

  // y = A x
  for (int i = 0; i < n; ++i) {
    real_t sum = 0.0;
    const real_t *Ai = &A[(size_t)i * (size_t)n];
    for (int j = 0; j < n; ++j) sum += Ai[j] * X[j];
    y[i] = sum;
  }

  // r = b - y  e ||r||_2
  real_t norm2 = 0.0;
  for (int i = 0; i < n; ++i) {
    r[i] = b[i] - y[i];
    norm2 += r[i] * r[i];
  }
  norm2 = sqrt(norm2);

  free(y);
  free(r);

  if (tempo) *tempo = timestamp() - *tempo;
  return norm2;
}
