#include "sislin.h"

static inline real_t generateRandomA(unsigned int i, unsigned int j, unsigned int k);
static inline real_t generateRandomB(unsigned int k);

/**
 * Função que gera os coeficientes de um sistema linear k-diagonal
 * @param i,j coordenadas do elemento a ser calculado (0<=i,j<n)
 * @param k numero de diagonais da matriz A
 */
static inline real_t generateRandomA(unsigned int i, unsigned int j, unsigned int k)
{
  static real_t invRandMax = 1.0 / (real_t)RAND_MAX;
  return ((i == j) ? (real_t)(k << 1) : 1.0) * (real_t)random() * invRandMax;
}

/**
 * Função que gera os termos independentes de um sistema linear k-diagonal
 * @param k numero de diagonais da matriz A
 */
static inline real_t generateRandomB(unsigned int k)
{
  static real_t invRandMax = 1.0 / (real_t)RAND_MAX;
  return (real_t)(k << 2) * (real_t)random() * invRandMax;
}

/* Cria matriz 'A' k-diagonal e Termos independentes B */
void criaKDiagonal(int n, int k, real_t **A, real_t **B) {
  // Aloca matriz A (n x n) linearizada e vetor B (n)
  *A = (real_t *) calloc(n * n, sizeof(real_t));
  *B = (real_t *) calloc(n, sizeof(real_t));

  if (*A == NULL || *B == NULL) {
      fprintf(stderr, "Erro ao alocar memória em criaKDiagonal\n");
      exit(1);
  }

  int offset = k / 2;  // número de diagonais abaixo/acima da principal

  for (int i = 0; i < n; i++) {
      // Preenche a linha i da matriz A
      for (int j = i - offset; j <= i + offset; j++) {
          if (j >= 0 && j < n) {
              (*A)[i * n + j] = generateRandomA(i, j, k);
          }
      }

      // Preenche o vetor B com valor aleatório
      (*B)[i] = generateRandomB(k);
  }
}

/* Gera matriz simetrica positiva */
void genSimetricaPositiva(real_t *A, real_t *b, int n, int k,
                          real_t **ASP, real_t *bsp, rtime_t *tempo)
{
  *tempo = timestamp();

  // Aloca ASP e bsp
  *ASP = (real_t *)calloc(n * n, sizeof(real_t));
  if (*ASP == NULL)
  {
    fprintf(stderr, "Erro ao alocar memória para ASP\n");
    exit(1);
  }

  for (int i = 0; i < n; i++)
    bsp[i] = 0.0;

  // Calcula ASP = Aᵗ·A  (produto matricial)
  for (int i = 0; i < n; i++)
  { // linha de Aᵗ (coluna de A)
    for (int j = 0; j < n; j++)
    { // coluna de A
      real_t sum = 0.0;
      for (int l = 0; l < n; l++)
      {                             // percorre linha l
        real_t a_li = A[l * n + i]; // Aᵗ[i][l] == A[l][i]
        real_t a_lj = A[l * n + j]; // A[l][j]
        sum += a_li * a_lj;
      }
      (*ASP)[i * n + j] = sum;
    }
  }

  // Calcula bsp = Aᵗ·b
  for (int i = 0; i < n; i++)
  { // Aᵗ[i][*]
    real_t sum = 0.0;
    for (int l = 0; l < n; l++)
    {
      sum += A[l * n + i] * b[l]; // Aᵗ[i][l] * b[l] == A[l][i] * b[l]
    }
    bsp[i] = sum;
  }

  *tempo = timestamp() - *tempo;
}

void geraDLU(real_t *A, int n, int k, real_t **D, real_t **L, real_t **U, real_t *tempo)
{
  *tempo = timestamp();

  *tempo = timestamp() - *tempo;
}

/**
 * Devolve matriz M⁻¹
 *
 */
void geraPreCond(real_t *D, real_t *L, real_t *U, real_t w, int n, int k, real_t **M, real_t *tempo)
{
  *tempo = timestamp();

  *tempo = timestamp() - *tempo;
}

real_t calcResiduoSL(real_t *A, real_t *b, real_t *X, int n, int k, real_t *tempo)
{
  *tempo = timestamp();

  real_t *r = calloc(n, sizeof(real_t));

  *tempo = timestamp() - *tempo;
}
