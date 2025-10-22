#include "helpers.h"

/* ---------- Operações vetoriais básicas ---------- */

// Retorna o produto interno de um vetor A por um vetor B
real_t vet_produto(int n, const real_t *x, const real_t *y)
{
    real_t s = 0.0;
    for (int i = 0; i < n; ++i)
        s += x[i] * y[i];
    return s;
}

// Calcula a norma 2 (euclidiana) de um vetor (raiz da soma dos quadrados dos elementos)
real_t vet_norma2(int n, const real_t *x)
{
    /* ||x||2 = sqrt(x^T x) */
    return sqrt(vet_produto(n, x, x));
}

// Copia os valores de um vetor para outro vetor
void vet_copy(int n, const real_t *fonte, real_t *dst)
{
    for (int i = 0; i < n; ++i)
        dst[i] = fonte[i];
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
    for (int i = 0; i < n; ++i)
        x[i] *= alpha;
}

/* ---------- Matriz densa ---------- */

// Multiplica matriz densa A (n×n) por vetor x: y = A*x
void matvet_densa(const real_t *A, const real_t *x, real_t *y, int n)
{
    for (int i = 0; i < n; ++i)
    {
        real_t s = 0.0;
        const real_t *Ai = &A[(size_t)i * (size_t)n]; /* linha i de A */
        for (int j = 0; j < n; ++j)
            s += Ai[j] * x[j];
        y[i] = s;
    }
}

/* ---------- PC de Jacobi ---------- */

// Extrai diagonal e seu inverso (checando zeros)
int extrai_diag_e_invD(const real_t *A, int n, int k, real_t *D,
                       real_t *invD, real_t eps)
{
    (void)k;
    for (int i = 0; i < n; ++i)
    {
        real_t d = A[IDX(i, i, n)];
        D[i] = d;
        if (fabs(d) < eps)
            return 1;
        invD[i] = 1.0 / d;
    }
    return 0;
}

// Aplica Jacobi: y = D⁻¹ r
void aplica_jacobi(int n, const real_t *invD, const real_t *r, real_t *y)
{
    for (int i = 0; i < n; ++i)
        y[i] = invD[i] * r[i];
}

/* ---------- Varreduras para SGS/SSOR ---------- */
// Baseados em (CARVALHO, VIEIRA, 2022)

// Resolve o sistema triangular inferior na forma (D*escala_diag + L) * y = vet_dir
void varredura_progressiva_DL(const real_t *A, int n, int k,
                              real_t escala_diag,
                              const real_t *vet_dir, real_t *y)
{
    const int p = k / 2;
    for (int i = 0; i < n; ++i)
    {
        real_t s = vet_dir[i];
        int j0 = (i - p > 0) ? (i - p) : 0; /* início da banda na linha i */
        for (int j = j0; j < i; ++j)
        {
            s -= A[IDX(i, j, n)] * y[j]; /* subtrai contribuição de L*y */
        }
        real_t ddiag = escala_diag * A[IDX(i, i, n)]; /* diagonal “D escalada” */
        y[i] = s / ddiag;
    }
}

// Resolve o sistema triangular superior na forma (D*escala_diag + U) * y = vet_dir
void varredura_regressiva_DU(const real_t *A, int n, int k,
                             real_t escala_diag,
                             const real_t *vet_dir, real_t *y)
{
    const int p = k / 2;
    for (int ii = 0; ii < n; ++ii)
    {
        int i = (n - 1) - ii;
        real_t s = vet_dir[i];
        int j1 = (i + p < n - 1) ? (i + p) : (n - 1); /* fim da banda na linha i */
        for (int j = i + 1; j <= j1; ++j)
        {
            s -= A[IDX(i, j, n)] * y[j]; /* subtrai contribuição de U*y */
        }
        real_t ddiag = escala_diag * A[IDX(i, i, n)]; /* diagonal “D escalada” */
        y[i] = s / ddiag;
    }
}

/* ---------- Calculo do Resíduo ---------- */

// Calcula ||r||2 com r = b - A x (usando A densa n×n)
real_t residuo_l2(const real_t *A, const real_t *b, const real_t *x, int n)
{
    real_t *Ax = (real_t *)calloc((size_t)n, sizeof(real_t));
    if (!Ax)
        return NAN;
    matvet_densa(A, x, Ax, n); /* Ax */
    real_t s2 = 0.0;
    for (int i = 0; i < n; ++i)
    {
        real_t ri = b[i] - Ax[i];
        s2 += ri * ri;
    }
    free(Ax);
    return sqrt(s2);
}