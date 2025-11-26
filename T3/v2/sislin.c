#include "sislin.h"
#include "helpers.h"

static void alloc_matdiag(matdiag_t *A, int n, int k);
static inline void matdiag_add_entry(matdiag_t *M, int row, int col, real_t val);

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

// Cria matriz A k-diagonal e Termos independentes B
void criaKDiagonal(int n, int k, real_t **A, real_t **B)
{
  *A = (real_t *)calloc((size_t)n * (size_t)n, sizeof(real_t));
  *B = (real_t *)calloc((size_t)n, sizeof(real_t));
  if (*A == NULL || *B == NULL)
  {
    fprintf(stderr, "Erro ao alocar memória em criaKDiagonal\n");
    exit(2);
  }

  // nº de diagonais abaixo/acima da principal
  int p = k / 2;

  for (int i = 0; i < n; i++)
  {
    // Preenche a banda da linha i
    int j0 = (i - p < 0) ? 0 : (i - p);
    int j1 = (i + p >= n) ? (n - 1) : (i + p);
    for (int j = j0; j <= j1; j++)
    {
      (*A)[IDX(i, j, n)] = generateRandomA((unsigned)i, (unsigned)j, (unsigned)k);
    }
    // Termo independente
    (*B)[i] = generateRandomB((unsigned)k);
  }
}

void criaKDiagonal_v2(int n, int k, matdiag_t *A, real_t **B)
{
  int p = k / 2; // número de diagonais acima/abaixo da principal

  alloc_matdiag(A, n, k);
  for (int d = 0; d < k; ++d)
    A->offsets[d] = d - p;

  // vetor b
  *B = (real_t *)calloc(n, sizeof(real_t));
  if (!*B)
  {
    fprintf(stderr, "Erro ao alocar vetor B\n");
    exit(2);
  }

  // preenche diagonais e vetor b
  for (int i = 0; i < n; ++i)
  {
    for (int d = 0; d < k; ++d)
    {
      int j = i + A->offsets[d];
      if (j >= 0 && j < n)
      {
        A->diag[d][i] = generateRandomA((unsigned)i, (unsigned)j, (unsigned)k);
      }
    }
    (*B)[i] = generateRandomB((unsigned)k);
  }
}

// Gera matriz simetrica positiva
void genSimetricaPositiva(real_t *A, real_t *b, int n,
                          real_t **ASP, real_t *bsp, rtime_t *tempo)
{
  *tempo = timestamp();

  // Alocação da matriz
  *ASP = (real_t *)calloc(n * n, sizeof(real_t));
  if (*ASP == NULL)
  {
    fprintf(stderr, "Erro ao alocar memória para ASP\n");
    exit(1);
  }

  // Atribui valor 0 para B simétrica positiva (Aᵗ·b) ou bsp
  for (int i = 0; i < n; i++)
    bsp[i] = 0.0;

  // Calcula ASP (A simétrica positiva (Aᵗ·A))
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

  // Calcula bsp (B simétrica positiva (Aᵗ·b))
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

void genSimetricaPositiva_diag(const matdiag_t *A, const real_t *b,
                               matdiag_t *ASP, real_t **bsp, rtime_t *tempo)
{
  if (!A || !b || !ASP || !bsp)
  {
    fprintf(stderr, "Parametros invalidos em genSimetricaPositiva_diag\n");
    exit(2);
  }

  const int n = A->n;
  const int p = A->k / 2;
  const int bandwidth = 2 * p;      // meia-largura de AᵗA
  const int k_spd = 2 * bandwidth + 1;

  alloc_matdiag(ASP, n, k_spd);
  for (int d = 0; d < k_spd; ++d)
    ASP->offsets[d] = d - bandwidth;

  *bsp = (real_t *)calloc(n, sizeof(real_t));
  if (!*bsp)
  {
    fprintf(stderr, "Erro ao alocar vetor bsp\n");
    exit(2);
  }

  if (tempo)
    *tempo = timestamp();

  // Calcula bsp = Aᵗ b reaproveitando a estrutura k-diagonal
  for (int row = 0; row < n; ++row)
  {
    for (int d = 0; d < A->k; ++d)
    {
      int col = row + A->offsets[d];
      if (col < 0 || col >= n)
        continue;

      real_t aval = A->diag[d][row];
      if (aval == 0.0)
        continue;

      (*bsp)[col] += aval * b[row];
    }
  }

  // Calcula ASP = Aᵗ A mantendo apenas as diagonais não nulas
  for (int row = 0; row < n; ++row)
  {
    for (int d1 = 0; d1 < A->k; ++d1)
    {
      int c1 = row + A->offsets[d1];
      if (c1 < 0 || c1 >= n)
        continue;

      real_t val1 = A->diag[d1][row];
      if (val1 == 0.0)
        continue;

      for (int d2 = 0; d2 < A->k; ++d2)
      {
        int c2 = row + A->offsets[d2];
        if (c2 < 0 || c2 >= n)
          continue;

        real_t val2 = A->diag[d2][row];
        if (val2 == 0.0)
          continue;

        matdiag_add_entry(ASP, c1, c2, val1 * val2);
      }
    }
  }

  if (tempo)
    *tempo = timestamp() - *tempo;
}

static void alloc_matdiag(matdiag_t *A, int n, int k)
{
  A->n = n;
  A->k = k;
  A->offsets = (int *)malloc(k * sizeof(int));
  A->diag = (real_t **)malloc(k * sizeof(real_t *));
  if (!A->offsets || !A->diag)
  {
    fprintf(stderr, "Erro ao alocar estrutura matdiag_t\n");
    exit(2);
  }

  for (int d = 0; d < k; ++d)
  {
    A->diag[d] = (real_t *)calloc(n, sizeof(real_t));
    if (!A->diag[d])
    {
      fprintf(stderr, "Erro ao alocar diagonal %d\n", d);
      exit(2);
    }
  }
}

static inline void matdiag_add_entry(matdiag_t *M, int row, int col, real_t val)
{
  if (!M || row < 0 || col < 0 || row >= M->n || col >= M->n)
    return;

  int offset = col - row;
  for (int d = 0; d < M->k; ++d)
  {
    if (M->offsets[d] == offset)
    {
      M->diag[d][row] += val;
      return;
    }
  }
}
