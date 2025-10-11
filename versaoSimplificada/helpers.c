/* helpers.c — utilidades numéricas usadas pelo solver
 *
 * Implementações efetivamente usadas pelo CG sem pré-condicionador:
 *   - vet_produto, vet_norma2, vet_copy, vet_preenche, vet_axpy, vet_scale
 *   - band_matvec (matriz k-diagonal × vetor)
 *
 * Stubs (comentados) para quando você ativar SGS/SSOR/Jacobi no PCG:
 *   - forward_sweep_DL, backward_sweep_DU, extract_diag_and_invd
 */

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "helpers.h"
#include "sislin.h"

/* ---------- Operações vetoriais básicas ---------- */

// Retorna o produto interno de um vetor A por um vetor B
real_t vet_produto(int n, const real_t *x, const real_t *y)
{
    real_t s = 0.0;
    for (int i = 0; i < n; ++i)
        s += x[i] * y[i];
    return s;
}

real_t vet_norma2(int n, const real_t *x)
{
    /* ||x||2 = sqrt(x^T x) */
    return sqrt(vet_produto(n, x, x));
}

// Copia os valores de um vetor para outro vetor
void vet_copy(int n, const real_t *fonte, real_t *dst)
{
    memcpy(dst, fonte, (size_t)n * sizeof(real_t));
}

// Seta os valores de um vetor com o valor val
void vet_preenche(int n, real_t val, real_t *x)
{
    for (int i = 0; i < n; ++i)
        x[i] = val;
}

// Soma y := y + alpha*x (combinação linear de vetores)
void vet_axpy(int n, real_t alpha, const real_t *x, real_t *y)
{
    for (int i = 0; i < n; ++i)
        y[i] += alpha * x[i];
}

// Multiplica todo o vetor x pelo escalar alpha (escala o vetor)
void vet_escala(int n, real_t alpha, real_t *x)
{
    /* x <- alpha*x */
    for (int i = 0; i < n; ++i)
        x[i] *= alpha;
}

/* ---------- Mat-vec otimizado para k-diagonal ---------- */

// calcula o produto matriz–vetor z=Ax explorando que A é k-diagona
void band_matvet(const real_t *A, int n, int k,
                 const real_t *x, real_t *z)
{
    /* A é k-diagonal: k = 2*p + 1, onde p é o raio da banda
       Para cada linha i, só existem coeficientes não nulos
       em j ∈ [max(0,i-p), min(n-1,i+p)].
    */
    const int p = k / 2;
    for (int i = 0; i < n; ++i)
    {
        real_t s = 0.0;
        int j0 = (i - p > 0) ? (i - p) : 0;
        int j1 = (i + p < n - 1) ? (i + p) : (n - 1);
        for (int j = j0; j <= j1; ++j)
        {
            s += A[IDX(i, j, n)] * x[j];
        }
        z[i] = s;
    }
}

/* ---------- Operações vetoriais básicas ---------- */

int extrai_diag_e_invD(const real_t *A, int n, int k, real_t *D,
                       real_t *invD, real_t eps)
{
    (void)k; /* não é necessário para extrair a diagonal */
    for (int i = 0; i < n; ++i)
    {
        real_t d = A[IDX(i, i, n)];
        D[i] = d;
        if (fabs(d) < eps)
            return 1; /* diagonal problemática */
        invD[i] = 1.0 / d;
    }
    return 0;
}

/* Aplicação do Jacobi por iteração */
void aplica_jacobi(int n, const real_t *invD, const real_t *r, real_t *y)
{
    for (int i = 0; i < n; ++i)
        y[i] = invD[i] * r[i];
}

/* Setup empacotado: aloca D/invD e preenche*/
int setup_jacobi(const real_t *A, int n, int k,
                 real_t **D_out, real_t **invD_out, real_t eps)
{
    real_t *D = (real_t *)malloc((size_t)n * sizeof(real_t));
    real_t *invD = (real_t *)malloc((size_t)n * sizeof(real_t));
    if (!D || !invD)
    {
        free(D);
        free(invD);
        return 2;
    }

    int rc = extrai_diag_e_invD(A, n, k, D, invD, eps);
    if (rc != 0)
    {
        free(D);
        free(invD);
        return rc;
    }

    *D_out = D;
    *invD_out = invD;
    return 0;
}

/* ---------- Varreduras para SGS/SSOR (PREPARO PARA O PRÓXIMO PASSO) ---------- */
/* As funções abaixo NÃO são usadas pelo CG “NONE”.
   Eu deixei implementações completas (útil para SGS/SSOR),
   mas você pode mantê-las aqui já prontas para quando migrar ao PCG.
*/

void forward_sweep_DL(const real_t *A, int n, int k,
                      real_t diagScale,
                      const real_t *rhs, real_t *t)
{
    /* Resolve (D*diagScale + L) t = rhs via varredura forward.
       Na prática: cada t[i] depende de t[0..i-1]. */
    const int p = k / 2;
    for (int i = 0; i < n; ++i)
    {
        real_t s = rhs[i];
        int j0 = (i - p > 0) ? (i - p) : 0;
        for (int j = j0; j < i; ++j)
        {
            s -= A[IDX(i, j, n)] * t[j]; /* subtrai L*t */
        }
        real_t ddiag = diagScale * A[IDX(i, i, n)];
        t[i] = s / ddiag;
    }
}

void backward_sweep_DU(const real_t *A, int n, int k,
                       real_t diagScale,
                       const real_t *rhs, real_t *y)
{
    /* Resolve (D*diagScale + U) y = rhs via varredura backward.
       Cada y[i] depende de y[i+1..n-1]. */
    const int p = k / 2;
    for (int ii = 0; ii < n; ++ii)
    {
        int i = (n - 1) - ii;
        real_t s = rhs[i];
        int j1 = (i + p < n - 1) ? (i + p) : (n - 1);
        for (int j = i + 1; j <= j1; ++j)
        {
            s -= A[IDX(i, j, n)] * y[j]; /* subtrai U*y */
        }
        real_t ddiag = diagScale * A[IDX(i, i, n)];
        y[i] = s / ddiag;
    }
}

