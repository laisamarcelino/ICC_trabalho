/* cg_none.c — Método dos Gradientes Conjugados (sem pré-condicionador)
 *
 * Esta função segue fielmente o Algoritmo 3.5 do material do professor.
 * Cada bloco do código indica, em comentários, a(s) linha(s) do pseudocódigo.
 *
 * CG sem PC + Jacobi dentro da MESMA função
 *
 * Mantém o mapeamento para o Algoritmo 3.5 (comentários [linha X]).
 * Diferenças do PCG (Alg. 3.6) estão marcadas como:  // JACOBI/PCG: ...
 *
 * Critério de parada exigido no trabalho:
 *   ||Δx||_∞ = max_i |x_i^(k+1) - x_i^(k)| < eps_inf
 *
 * Requisitos:
 *   - A é uma matriz n×n simétrica definida positiva (SPD), k-diagonal.
 *   - b é o vetor de termos independentes.
 *   - x chega como vetor de saída (pode ser zerado aqui).
 *
 * Dependências:
 *   - helpers.h: vet_produto, vet_axpy, vet_copy, vet_norma2, band_matvet etc.
 *   - sislin.h: define real_t e o macro IDX(i,j,n) se usado em helpers.c.
 */

#include "pcgc.h"

int cg_solve_none(const real_t *A, const real_t *b, real_t *x,
                  int n, int k, int maxit, real_t eps_inf,
                  pcg_precond_t M /* NEW: precondicionador: NONE ou JACOBI */)
{
    /* buffers comuns: r (resíduo), v (direção), z (A v) */
    real_t *r = (real_t *)calloc((size_t)n, sizeof(real_t));
    real_t *v = (real_t *)calloc((size_t)n, sizeof(real_t));
    real_t *z = (real_t *)calloc((size_t)n, sizeof(real_t));
    if (!r || !v || !z)
    {
        fprintf(stderr, "[cg/pcg] Falha de alocacao de memoria (r/v/z)\n");
        free(r);
        free(v);
        free(z);
        return -1;
    }

    /* ---------- (PREP) Dados do pré-condicionador (extensível) ---------- */
    /* Para JACOBI precisamos de: y = M^{-1} r e invD = 1/diag(A) */
    real_t *y = NULL;               /* JACOBI/PCG: armazena M^{-1}r */
    real_t *D = NULL, *invD = NULL; /* JACOBI/PCG: diagonal e seu inverso */

    if (M == PCG_PRECOND_JACOBI)
    {
        y = (real_t *)calloc((size_t)n, sizeof(real_t));
        D = (real_t *)malloc((size_t)n * sizeof(real_t));
        invD = (real_t *)malloc((size_t)n * sizeof(real_t));
        if (!y || !D || !invD)
        {
            fprintf(stderr, "[pcg] Falha de alocacao (Jacobi)\n");
            free(r);
            free(v);
            free(z);
            free(y);
            free(D);
            free(invD);
            return -2;
        }
        /* prepara Jacobi uma vez: D e invD (falha graciosa se diagonal ~0) */
        if (extrai_diag_e_invD(A, n, k, D, invD, 1e-30) != 0)
        {
            fprintf(stderr, "[pcg] Jacobi inviavel: diagonal nula/pequena\n");
            free(r);
            free(v);
            free(z);
            free(y);
            free(D);
            free(invD);
            return -3;
        }
    }
    /* (No futuro: trocar este bloco por pcg_setup_precond(...) + switch) */

    /* ------------------------------------------------------------------
    [linha 1]  x^(0) = 0,  r = b,  aux = r^T r,  v = b
    ------------------------------------------------------------------ */
    memset(x, 0, (size_t)n * sizeof(real_t)); /* x^(0) = 0 */
    vet_copy(n, b, r);                        /* r = b     */

    real_t aux = 0.0;
    if (M == PCG_PRECOND_NONE)
    {
        aux = vet_produto(n, r, r); /* aux = r^T r (Alg. 3.5) */
        vet_copy(n, r, v);          /* v = r (Alg. 3.5) */
    }
    else
    {                                 /* JACOBI/PCG: diferenças em relação ao Algoritmo 3.5 */
        aplica_jacobi(n, invD, r, y); /* JACOBI/PCG: y = M^{-1} r */
        aux = vet_produto(n, r, y);   /* JACOBI/PCG: aux = r^T y */
        vet_copy(n, y, v);            /* JACOBI/PCG: v = y */
    }

    /* ------------------------------------------------------------------
    [linha 2]  Para k = 0 : max, faça
    ------------------------------------------------------------------ */
    for (int it = 0; it < maxit; ++it)
    {
        /* [linha 3]  z = A v   (multiplicação banda × vetor) */
        band_matvet(A, n, k, v, z);

        /* [linha 4]  s = aux / (v^T z)  */
        real_t vTz = vet_produto(n, v, z);

        /* robustez numérica: se denom ~0 e “energia” do resíduo já é pequena, finalize */
        if (fabs(vTz) < 1e-30)
        {
            /* NONE: sqrt(r^T r);  JACOBI/PCG: sqrt(r^T y) (aux) */
            real_t res_energy = sqrt(fabs(aux));
            if (res_energy < 1e-12)
            {
                free(r);
                free(v);
                free(z);
                free(y);
                free(D);
                free(invD);
                return it; /* converge sem atualizar x nesta iteração */
            }
            fprintf(stderr, "[cg/pcg] v^T z ~ 0 (divisao instavel)\n");
            break;
        }

        real_t s = aux / vTz;

        /* [linha 5]  x^(k+1) = x^(k) + s v
        (e computa ||Δx||_∞ = max_i |s*v[i]| para o critério do TRABALHO) */
        real_t dx_max = 0.0;
        for (int i = 0; i < n; ++i)
        {
            real_t dx_i = s * v[i];
            x[i] += dx_i;
            real_t adx = fabs(dx_i);
            if (adx > dx_max)
                dx_max = adx;
        }

        /* [linha 6]  r = r − s z  (AXPY com alfa = -s) */
        vet_axpy(n, -s, z, r);

        /* *** Critério do TRABALHO ***: ||Δx||_∞ < eps */
        if (dx_max < eps_inf)
        {
            free(r);
            free(v);
            free(z);
            free(y);
            free(D);
            free(invD);
            return it + 1;
        }

        /* [linha 7]–[linha 11]
        Algoritmo 3.5 (CG “puro”) vs Algoritmo 3.6 (PCG Jacobi):
        - NONE:
         [7]  aux1 = r^T r
         [10] m = aux1/aux; aux = aux1
         [11] v = r + m v
        - JACOBI/PCG:
         [7]  y = M^{-1} r
         [8]  aux1 = r^T y
         [10] m = aux1/aux; aux = aux1
         [11] v = y + m v
        */
        if (M == PCG_PRECOND_NONE)
        {
            real_t aux1 = vet_produto(n, r, r); /* [linha 7] aux1 = r^T r */
            real_t m = aux1 / aux;              /* [linha 10] m = aux1/aux */
            aux = aux1;                         /*               aux = aux1 */
            for (int i = 0; i < n; ++i)         /* [linha 11] v = r + m v */
                v[i] = r[i] + m * v[i];
        }
        else
        {                                       /* JACOBI/PCG: diferenças em relação ao Algoritmo 3.5 */
            aplica_jacobi(n, invD, r, y);       /* [linha 7]  y = M^{-1} r */
            real_t aux1 = vet_produto(n, r, y); /* [linha 8]  aux1 = r^T y */
            real_t m = aux1 / aux;              /* [linha 10] m = aux1/aux */
            aux = aux1;                         /*               aux = aux1 */
            for (int i = 0; i < n; ++i)         /* [linha 11] v = y + m v */
                v[i] = y[i] + m * v[i];
        }
    }

    /* atingiu limite de iterações */
    free(r);
    free(v);
    free(z);
    free(y);
    free(D);
    free(invD);
    return maxit; /* não atingiu tolerância por Δx∞ dentro de maxit */
}
