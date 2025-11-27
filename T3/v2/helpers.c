#include "helpers.h"

/* ---------- Operações vetoriais básicas ---------- */

// Retorna o produto interno de um vetor A por um vetor B
real_t vet_produto(int n, const real_t *x, const real_t *y)
{
    real_t s = 0.0;
    int i;

    // Unroll de 4 elementos
    for (i = 0; i < n - n % 4; i += 4)
    {
        s += x[i] * y[i] + x[i + 1] * y[i + 1] + x[i + 2] * y[i + 2] + x[i + 3] * y[i + 3];
    }

    // Resto (se n não múltiplo de 4)
    for (; i < n; ++i)
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
    int i;

    // Unroll de 4 elementos
    for (i = 0; i < n - n % 4; i += 4)
    {
        dst[i] = fonte[i];
        dst[i + 1] = fonte[i + 1];
        dst[i + 2] = fonte[i + 2];
        dst[i + 3] = fonte[i + 3];
    }

    for (; i < n; ++i)
        dst[i] = fonte[i];
}

// Seta os valores de um vetor com o valor val
void vet_preenche(int n, real_t val, real_t *x)
{
    for (int i = 0; i < n; ++i)
        x[i] = val;
}

void vet_axpy(int n, real_t alpha, const real_t *x, real_t *y)
{
    int i;

    // Unroll de 4 elementos
    for (i = 0; i < n - n % 4; i += 4)
    {
        y[i] += alpha * x[i];
        y[i + 1] += alpha * x[i + 1];
        y[i + 2] += alpha * x[i + 2];
        y[i + 3] += alpha * x[i + 3];
    }

    for (; i < n; ++i)
        y[i] += alpha * x[i];
}

// Multiplica todo o vetor x pelo escalar alpha (escala o vetor)
void vet_escala(int n, real_t alpha, real_t *x)
{
    int i;

    // Unroll de 4 elementos
    for (i = 0; i < n - n % 4; i += 4)
    {
        x[i] *= alpha;
        x[i + 1] *= alpha;
        x[i + 2] *= alpha;
        x[i + 3] *= alpha;
    }

    for (; i < n; ++i)
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
void varredura_progressiva_DL(const matdiag_t *A, int n, int k,
                              real_t escala_diag,
                              const real_t *vet_dir, real_t *y)
{
    /* Encontra a posição da diagonal principal (offset == 0) */
    int pos_diag_principal = -1;
    for (int d = 0; d < k; ++d)
        if (A->offsets[d] == 0)
        {
            pos_diag_principal = d;
            break;
        }

    if (pos_diag_principal == -1)
    {
        fprintf(stderr, "[varredura_progressiva_DL] erro: matriz sem diagonal principal!\n");
        return;
    }

    /* Resolve (D*escala_diag + L) * y = vet_dir
       Percorre as diagonais da matriz armazenada em bandas para subtrair L*y */
    for (int i = 0; i < n; ++i)
    {
        real_t s = vet_dir[i];

        /* percorre todas as diagonais e subtrai as contribuições de colunas j < i (L) */
        for (int d = 0; d < k; ++d)
        {
            int j = i + A->offsets[d];
            if (j >= 0 && j < i)
            {
                s -= A->diag[d][i] * y[j];
            }
        }

        real_t ddiag = escala_diag * A->diag[pos_diag_principal][i]; /* diagonal “D escalada” */
        y[i] = s / ddiag;
    }
}

// Resolve o sistema triangular superior na forma (D*escala_diag + U) * y = vet_dir
void varredura_regressiva_DU(const matdiag_t *A, int n, int k,
                             real_t escala_diag,
                             const real_t *vet_dir, real_t *y)
{
    /* Encontra a posição da diagonal principal (offset == 0) */
    int pos_diag_principal = -1;
    for (int d = 0; d < k; ++d)
        if (A->offsets[d] == 0)
        {
            pos_diag_principal = d;
            break;
        }

    if (pos_diag_principal == -1)
    {
        fprintf(stderr, "[varredura_regressiva_DU] erro: matriz sem diagonal principal!\n");
        return;
    }

    /* Resolve (D*escala_diag + U) * y = vet_dir
       Varredura regressiva: percorre linhas de baixo para cima e subtrai U*y */
    for (int ii = 0; ii < n; ++ii)
    {
        int i = (n - 1) - ii;
        real_t s = vet_dir[i];

        /* percorre todas as diagonais e subtrai as contribuições de colunas j > i (U) */
        for (int d = 0; d < k; ++d)
        {
            int j = i + A->offsets[d];
            if (j > i && j < n)
            {
                s -= A->diag[d][i] * y[j];
            }
        }

        real_t ddiag = escala_diag * A->diag[pos_diag_principal][i]; /* diagonal “D escalada” */
        y[i] = s / ddiag;
    }
}

void matvet_diagonais(const matdiag_t *A, const real_t *x, real_t *y)
{
    int n = A->n, k = A->k;

    for (int i = 0; i < n; ++i)
        y[i] = 0.0;

    for (int d = 0; d < k; ++d)
    {
        int offset = A->offsets[d];
        real_t *diag = A->diag[d];

        if (offset >= 0)
        {
            // i de 0 até n - offset - 1, pois x[i + offset] deve ser < n
            int limite = n - offset;
            for (int i = 0; i < limite; ++i)
                y[i] += diag[i] * x[i + offset];
        }
        else
        {
            // i de -offset até n-1, pois x[i + offset] >= 0
            int start = -offset;
            for (int i = start; i < n; ++i)
                y[i] += diag[i] * x[i + offset];
        }
    }
}

/* ---------- Calculo do Resíduo ---------- */

real_t residuo_l2_v2(const matdiag_t *A, const real_t *b, const real_t *x)
{
    if (!A || !b || !x)
        return NAN;

    int n = A->n;
    real_t *Ax = (real_t *)calloc((size_t)n, sizeof(real_t));
    if (!Ax)
        return NAN;

    matvet_diagonais(A, x, Ax);

    real_t s2 = 0.0;
    int i = 0;

    for (; i <= n - 4; i += 4)
    {
        real_t r0 = b[i]     - Ax[i];
        real_t r1 = b[i + 1] - Ax[i + 1];
        real_t r2 = b[i + 2] - Ax[i + 2];
        real_t r3 = b[i + 3] - Ax[i + 3];

        s2 += r0 * r0 + r1 * r1 + r2 * r2 + r3 * r3;
    }

    for (; i < n; ++i)
    {
        real_t ri = b[i] - Ax[i];
        s2 += ri * ri;
    }

    free(Ax);
    return sqrt(s2);
}

int extrai_diag_e_invD_diag(const matdiag_t *A, real_t *D, real_t *invD, real_t eps)
{
    int pos_diag_principal = -1;
    for (int d = 0; d < A->k; ++d)
        if (A->offsets[d] == 0)
            pos_diag_principal = d;

    if (pos_diag_principal == -1)
        return -1;

    real_t *diag = A->diag[pos_diag_principal];
    int n = A->n;

    for (int i = 0; i < n; ++i)
    {
        D[i] = diag[i];
        if (fabs(D[i]) < eps)
            return -2;
        invD[i] = 1.0 / D[i];
    }

    return 0;
}

void liberaMatDiag(matdiag_t *A)
{
    for (int d = 0; d < A->k; ++d)
        free(A->diag[d]);
    free(A->diag);
    free(A->offsets);
}
