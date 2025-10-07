// sislin.c — Implementações auxiliares para o cgSolver
// ----------------------------------------------------
// Este módulo provê rotinas de geração/manipulação de sistemas lineares
// bandados (k-diagonais) e utilidades para pré-condicionamento.
// Todas as matrizes são armazenadas em formato denso, linearizadas
// (vetor de tamanho n*n), com acesso via IDX(i,j,n).

#include "sislin.h"
#include <float.h>   // para DBL_EPSILON

// Macro de indexação (linha-major): A[i][j] == A[ i*n + j ]
#define IDX(i,j,n) ((i)*(n) + (j))

// Limite para tratar diagonal "quase nula"
#ifndef DLIMIT
#define DLIMIT (1e-14)
#endif

// ---------- Protótipos internos (estáticos) ---------- //
static inline real_t generateRandomA(unsigned int i, unsigned int j, unsigned int k);
static inline real_t generateRandomB(unsigned int k);
static void matvec_band(const real_t *A, const real_t *x, real_t *y, int n, int k);
static void matvec_dense(const real_t *A, const real_t *x, real_t *y, int n);
static void forward_solve(const real_t *L, const real_t *b, real_t *x, int n);
static void backward_solve(const real_t *U, const real_t *b, real_t *x, int n);

// =====================================================
// Geração pseudo-aleatória dos coeficientes
// -----------------------------------------------------
// IMPORTANTE: o enunciado exige que o programa chame srandom(20252)
// uma única vez antes de usar estas funções.
// =====================================================

/**
 * Gera elemento A[i,j] de uma matriz k-diagonal.
 * - Diagonal principal recebe valores proporcionalmente maiores (para
 *   facilitar predominância diagonal média; ainda assim, A NÃO é SPD).
 * - Fora da diagonal, valores aleatórios menores.
 * Obs: 'k' é o nº total de diagonais não nulas (ímpar), logo o semi-largura
 * de banda é p = k/2.
 */
static inline real_t generateRandomA(unsigned int i, unsigned int j, unsigned int k) {
  static real_t invRandMax = 1.0 / (real_t)RAND_MAX;
  if (i == j) {
    return (real_t)(k * 2) * (real_t)random() * invRandMax;
  } else {
    return 1.0 * (real_t)random() * invRandMax;
  }
}

/**
 * Gera termo independente b[i] coerente com a escala típica de A.
 */
static inline real_t generateRandomB(unsigned int k) {
  static real_t invRandMax = 1.0 / (real_t)RAND_MAX;
  return (real_t)(k * 4) * (real_t)random() * invRandMax;
}

// =====================================================
// criaKDiagonal
// -----------------------------------------------------
// Constrói matriz A (n x n) k-diagonal (armazenada densa) e vetor b.
// O semi-largura de banda é p = k/2, e apenas colunas j em [i-p, i+p]
// são preenchidas; demais ficam zero.
// =====================================================
void criaKDiagonal(int n, int k, real_t **A, real_t **B) {
  *A = (real_t *) calloc((size_t)n * (size_t)n, sizeof(real_t));
  *B = (real_t *) calloc((size_t)n, sizeof(real_t));
  if (*A == NULL || *B == NULL) {
    fprintf(stderr, "Erro ao alocar memória em criaKDiagonal\n");
    exit(2);
  }

  int p = k / 2; // nº de diagonais abaixo/acima da principal

  for (int i = 0; i < n; i++) {
    // Preenche a banda da linha i
    int j0 = (i - p < 0) ? 0 : (i - p);
    int j1 = (i + p >= n) ? (n - 1) : (i + p);
    for (int j = j0; j <= j1; j++) {
      (*A)[IDX(i,j,n)] = generateRandomA((unsigned)i, (unsigned)j, (unsigned)k);
    }
    // Termo independente
    (*B)[i] = generateRandomB((unsigned)k);
  }
}

// =====================================================
// genSimetricaPositiva
// -----------------------------------------------------
// O método de Gradientes Conjugados (CG) requer A simétrica definida
// positiva (SPD). Se o A original não for SPD, transformamos o sistema
// Ax=b em (A^T A) x = (A^T b). A^T A é SPD (exceto casos degenerados).
//
// Observação: mesmo que A seja k-diagonal, A^T A costuma ficar mais densa
// (a largura de banda dobra aprox.). Aqui formamos A_spd densa por simplicidade.
// =====================================================
// =====================================================
// genSimetricaPositiva (corrigido: bsp alocado)
// =====================================================
void genSimetricaPositiva(real_t *A, real_t *b, int n, int k,
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


// =====================================================
// geraDLU
// -----------------------------------------------------
// Separa A em três componentes: A = L + D + U,
// onde D é diagonal, L é estritamente inferior, U é estritamente superior.
// Mantemos armazenamento denso por simplicidade, mas preenchendo apenas
// as posições relevantes (aproveitando a banda k-diagonal).
// =====================================================
void geraDLU(real_t *A, int n, int k,
             real_t **D, real_t **L, real_t **U, real_t *tempo) {
  *tempo = timestamp();

  *D = (real_t *) calloc((size_t)n * (size_t)n, sizeof(real_t));
  *L = (real_t *) calloc((size_t)n * (size_t)n, sizeof(real_t));
  *U = (real_t *) calloc((size_t)n * (size_t)n, sizeof(real_t));
  if (*D == NULL || *L == NULL || *U == NULL) {
    fprintf(stderr, "Erro ao alocar memória em geraDLU\n");
    exit(2);
  }

  int p = k / 2; // semi-largura de banda

  for (int i = 0; i < n; i++) {
    int j0 = (i - p < 0) ? 0 : (i - p);
    int j1 = (i + p >= n) ? (n - 1) : (i + p);
    for (int j = j0; j <= j1; j++) {
      real_t aij = A[IDX(i,j,n)];
      if (i == j) {
        (*D)[IDX(i,j,n)] = aij;
      } else if (j < i) {
        (*L)[IDX(i,j,n)] = aij; // abaixo da diagonal
      } else { // j > i
        (*U)[IDX(i,j,n)] = aij; // acima da diagonal
      }
    }
  }

  *tempo = timestamp() - *tempo;
}

// =====================================================
// geraPreCond
// -----------------------------------------------------
// Constrói **M^{-1}** conforme ω:
//   ω = -1   → sem pré-condicionador:      M^{-1} = I
//   ω =  0   → Jacobi:                     M^{-1} = D^{-1}
//   ω =  1   → Gauss-Seidel (opcional):    M^{-1} = (D+U)^{-1} D (D+L)^{-1}
//   1<ω< 2   → SSOR (opcional):            M^{-1} = (D+ωU)^{-1} D (D+ωL)^{-1}
//
// Observação importante: para GS/SSOR a forma mais eficiente em CG é
// aplicar y = M^{-1}r via duas substituições (sem formar M^{-1}).
// Aqui, para simplificar integração com o restante do trabalho, **formamos
// explicitamente** M^{-1} (custo O(n^3)). Em problemas muito grandes,
// recomenda-se trocar por operador-aplicação (sem formar M^{-1}).
// =====================================================
void geraPreCond(real_t *D, real_t *L, real_t *U, real_t w,
                 int n, int k, real_t **M, real_t *tempo) {
  (void)k; // não usamos k diretamente aqui
  *tempo = timestamp();

  // Aloca M^{-1}
  *M = (real_t *) calloc((size_t)n * (size_t)n, sizeof(real_t));
  if (*M == NULL) {
    fprintf(stderr, "Erro ao alocar memória para M^{-1}\n");
    exit(2);
  }

  // Caso 1: sem pré-condicionador → M^{-1} = I
  if (w < 0.0) {
    for (int i = 0; i < n; i++) {
      (*M)[IDX(i,i,n)] = 1.0;
    }
    *tempo = timestamp() - *tempo;
    return;
  }

  // Caso 2: Jacobi → M^{-1} = D^{-1} (apenas diagonal)
  if (fabs(w) <= DBL_EPSILON) {
    for (int i = 0; i < n; i++) {
      real_t dii = D[IDX(i,i,n)];
      if (fabs(dii) < DLIMIT) {
        fprintf(stderr, "Pré-condicionador Jacobi: elemento D[%d,%d] ~ 0\n", i, i);
        exit(2);
      }
      (*M)[IDX(i,i,n)] = 1.0 / dii;
    }
    *tempo = timestamp() - *tempo;
    return;
  }

  // Casos 3/4: GS (w=1) ou SSOR (1<w<2) — Opcional
  // Construímos DL = D + wL (triangular inferior) e DU = D + wU (triangular superior),
  // ambos denso/linearizado para simplicidade.
  real_t *DL = (real_t *) calloc((size_t)n * (size_t)n, sizeof(real_t));
  real_t *DU = (real_t *) calloc((size_t)n * (size_t)n, sizeof(real_t));
  if (DL == NULL || DU == NULL) {
    fprintf(stderr, "Erro ao alocar DL/DU para SSOR\n");
    exit(2);
  }

  // Monta DL e DU
  for (int i = 0; i < n; i++) {
    for (int j = 0; j <= i; j++) {
      real_t lij = L[IDX(i,j,n)];
      if (i == j) lij = 0.0; // L tem diagonal nula por definição
      DL[IDX(i,j,n)] = (i == j ? D[IDX(i,i,n)] : 0.0) + w * lij;
    }
    for (int j = i; j < n; j++) {
      real_t uij = U[IDX(i,j,n)];
      if (i == j) uij = 0.0; // U tem diagonal nula por definição
      DU[IDX(i,j,n)] = (i == j ? D[IDX(i,i,n)] : 0.0) + w * uij;
    }
  }

  // Verifica diagonais não nulas em DL e DU
  for (int i = 0; i < n; i++) {
    if (fabs(DL[IDX(i,i,n)]) < DLIMIT || fabs(DU[IDX(i,i,n)]) < DLIMIT) {
      fprintf(stderr, "SSOR/GS: diagonal ~0 em DL/DU na linha %d\n", i);
      exit(2);
    }
  }

  // Para formar M^{-1} explicitamente: para cada coluna j de M^{-1}
  // resolvemos (D+ωL) z = e_j   (substituição direta)
  //            y = D z          (multiplicação por diagonal)
  //            (D+ωU) x = y     (substituição inversa)
  // e definimos col_j(M^{-1}) = x.
  real_t *ej = (real_t *) calloc((size_t)n, sizeof(real_t));
  real_t *z  = (real_t *) calloc((size_t)n, sizeof(real_t));
  real_t *y  = (real_t *) calloc((size_t)n, sizeof(real_t));
  real_t *x  = (real_t *) calloc((size_t)n, sizeof(real_t));
  if (!ej || !z || !y || !x) {
    fprintf(stderr, "Erro ao alocar vetores temporários para SSOR/GS\n");
    exit(2);
  }

  for (int j = 0; j < n; j++) {
    // e_j
    for (int t = 0; t < n; t++) ej[t] = 0.0;
    ej[j] = 1.0;

    // (D+ωL) z = e_j
    // forward_solve espera matriz triangular inferior (com diagonal não nula)
    forward_solve(DL, ej, z, n);

    // y = D z (D é diagonal: y_i = D_ii * z_i)
    for (int i = 0; i < n; i++) {
      y[i] = D[IDX(i,i,n)] * z[i];
    }

    // (D+ωU) x = y
    backward_solve(DU, y, x, n);

    // Grava coluna j em M^{-1}
    for (int i = 0; i < n; i++) {
      (*M)[IDX(i,j,n)] = x[i];
    }
  }

  free(DL); free(DU);
  free(ej); free(z); free(y); free(x);

  *tempo = timestamp() - *tempo;
}

// =====================================================
// calcResiduoSL
// -----------------------------------------------------
// Calcula o resíduo r = b - A x e retorna ||r||_2 (norma Euclidiana).
// Para acelerar, quando informado k (banda), podemos limitar a soma
// por linha à vizinhança [i - k/2, i + k/2].
// =====================================================
// =====================================================
// calcResiduoSL (corrigido: sempre denso para SPD)
// =====================================================
real_t calcResiduoSL(real_t *A, real_t *b, real_t *X, int n, int k, real_t *tempo) {
  *tempo = timestamp();

  real_t *y = (real_t *) calloc((size_t)n, sizeof(real_t));
  if (!y) {
    fprintf(stderr, "Erro ao alocar vetor temporário em calcResiduoSL\n");
    exit(2);
  }

  // y = A x
  if (k > 1) {
    matvec_band(A, X, y, n, k);  // para A k-diagonal original
  } else {
    matvec_dense(A, X, y, n);    // para A SPD (já densa)
  }

  // r = b - y
  real_t nrm2 = 0.0;
  for (int i = 0; i < n; i++) {
    real_t ri = b[i] - y[i];
    nrm2 += ri * ri;
  }
  nrm2 = sqrt(nrm2);

  free(y);
  *tempo = timestamp() - *tempo;
  return nrm2;
}


// =====================================================
// Funções auxiliares internas
// =====================================================

// y = A x (matriz cheia densa)
static void matvec_dense(const real_t *A, const real_t *x, real_t *y, int n) {
  for (int i = 0; i < n; i++) {
    real_t acc = 0.0;
    for (int j = 0; j < n; j++) {
      acc += A[IDX(i,j,n)] * x[j];
    }
    y[i] = acc;
  }
}

// y = A x, aproveitando a estrutura k-diagonal (banda)
// A supõe-se com semi-largura p = k/2; fora da banda é zero.
static void matvec_band(const real_t *A, const real_t *x, real_t *y, int n, int k) {
  int p = k / 2;
  for (int i = 0; i < n; i++) {
    int j0 = (i - p < 0) ? 0 : (i - p);
    int j1 = (i + p >= n) ? (n - 1) : (i + p);
    real_t acc = 0.0;
    for (int j = j0; j <= j1; j++) {
      acc += A[IDX(i,j,n)] * x[j];
    }
    y[i] = acc;
  }
}

// Resolve L x = b (L triangular inferior com diagonal não nula)
static void forward_solve(const real_t *L, const real_t *b, real_t *x, int n) {
  for (int i = 0; i < n; i++) {
    real_t sum = 0.0;
    for (int j = 0; j < i; j++) {
      sum += L[IDX(i,j,n)] * x[j];
    }
    real_t lii = L[IDX(i,i,n)];
    if (fabs(lii) < DLIMIT) {
      fprintf(stderr, "forward_solve: diagonal ~0 na linha %d\n", i);
      exit(2);
    }
    x[i] = (b[i] - sum) / lii;
  }
}

// Resolve U x = b (U triangular superior com diagonal não nula)
static void backward_solve(const real_t *U, const real_t *b, real_t *x, int n) {
  for (int i = n - 1; i >= 0; i--) {
    real_t sum = 0.0;
    for (int j = i + 1; j < n; j++) {
      sum += U[IDX(i,j,n)] * x[j];
    }
    real_t uii = U[IDX(i,i,n)];
    if (fabs(uii) < DLIMIT) {
      fprintf(stderr, "backward_solve: diagonal ~0 na linha %d\n", i);
      exit(2);
    }
    x[i] = (b[i] - sum) / uii;
  }
}
