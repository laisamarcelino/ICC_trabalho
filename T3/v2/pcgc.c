/* Método dos Gradientes Conjugados */

#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#ifdef LIKWID_PERFMON
#include <likwid.h>
#endif

#include "pcgc.h"

/* ------------------------ Setup do pré-condicionador ------------------------ */

int pcg_setup(const matdiag_t *A, int n, pcg_precond_t M, real_t omega, pcg_contexto_t *contexto)
{
    // Zera toda a struct para facilitar limpeza em caso de erro
    memset(contexto, 0, sizeof(*contexto));
    contexto->type = M;

    // SGS: força ω=1.0
    // SSOR: exige 0<ω<2 (do enunciado e teoria)

    // Se é Gauss-Seidel, é forçado ω = 1.0
    if (M == PCG_PRECOND_SGS)
    {
        contexto->omega = 1.0;
    }

    // Se é SSOR, ω é verificado se está 0.0 < ω < 2.0
    else if (M == PCG_PRECOND_SSOR)
    {
        if (!(omega > 0.0 && omega < 2.0))
        {
            fprintf(stderr, "[pcg] omega invalido p/ SSOR: %.6g (esperado 0<omega<2)\n", omega);
            return -5;
        }
        contexto->omega = omega;
    }

    // Senão, é ignorada a atribuição de ω (ω = -1.0 sem pré-condicionador ou ω = 0.0 Jacobi)
    else
    {
        contexto->omega = omega;
    }

    if (M == PCG_PRECOND_NONE)
        return 0;

    // y serve para armazenar M⁻¹r de qualquer PC
    contexto->y = (real_t *)calloc((size_t)n, sizeof(real_t));
    if (!contexto->y)
    {
        pcg_free(contexto);
        return -1;
    }

    // Se for PC de Jacobi, precisamos extrair a diagonal principal e sua inversa
    if (M == PCG_PRECOND_JACOBI)
    {
        contexto->D = (real_t *)malloc((size_t)n * sizeof(real_t));
        contexto->invD = (real_t *)malloc((size_t)n * sizeof(real_t));
        if (!contexto->D || !contexto->invD)
        {
            pcg_free(contexto);
            return -2;
        }

        if (extrai_diag_e_invD_diag(A, contexto->D, contexto->invD, 1e-30) != 0)
        {
            fprintf(stderr, "[pcg] Jacobi inviavel: diagonal nula/pequena\n");
            pcg_free(contexto);
            return -3;
        }
        return 0;
    }

    // SGS/SSOR: precisamos de t e u para armazenar resultados intermediários durante as etapas
    // (varredura sobre o sistema linear modificado)
    contexto->t = (real_t *)calloc((size_t)n, sizeof(real_t));
    contexto->u = (real_t *)calloc((size_t)n, sizeof(real_t));
    if (!contexto->t || !contexto->u)
    {
        pcg_free(contexto);
        return -4;
    }

    return 0;
}

/* ------------------------ Aplicação do pré-condicionador -------------------- */

void pcg_apply(const matdiag_t *A, int n, int k, const pcg_contexto_t *contexto, const real_t *r, real_t *y)
{
    switch (contexto->type)
    {
    // Se não houver PC, apenas copia o resíduo para y = M⁻¹r, que resulta em r
    case PCG_PRECOND_NONE:
        vet_copy(n, r, y);
        break;

    // Se o PC for de Jacobi, faz y = D⁻¹r
    case PCG_PRECOND_JACOBI:
        aplica_jacobi(n, contexto->invD, r, y);
        break;

    // Nos dois métodos abaixo, não transformamos M explicitamente
    // Caso for GS, faz varreduras progressiva e regressiva:
    // (D + L) t = r
    // u = D t
    // (D + U) y = D t

    // Caso for SSOR, faz varreduras progressiva e regressiva (M_ω = (D + ωL) D⁻¹ (D + ωU)):
    // (D/ω + L) t = r
    // u = D t
    // (D/ω + U) y = u

    // Os métodos são semelhantes, o que muda é o peso da diagonal D nas etapas de sweep
    // Para GS, escala_diag = 1.0
    // Para SSOR, escala_diag = 1/ω
    case PCG_PRECOND_SGS:
    case PCG_PRECOND_SSOR:
    {
        real_t escala_diag = 1.0;
        if (contexto->type == PCG_PRECOND_SSOR)
            escala_diag = 1.0 / contexto->omega;

        // (D*escala_diag + L) t = r
        varredura_progressiva_DL(A, n, k, escala_diag, r, contexto->t);

        // Identifica a posição da diagonal principal
        int pos_diag_principal = -1;
        for (int d = 0; d < A->k; ++d)
            if (A->offsets[d] == 0)
            {
                pos_diag_principal = d;
                break;
            }

        if (pos_diag_principal == -1)
        {
            fprintf(stderr, "[pcg_apply] erro: matriz sem diagonal principal!\n");
            return;
        }

        // u = D * t  → usando apenas a diagonal principal armazenada
        real_t *diag_principal = A->diag[pos_diag_principal];
        for (int i = 0; i < n; ++i)
            contexto->u[i] = diag_principal[i] * contexto->t[i];

        // (D*escala_diag + U) y = u
        varredura_regressiva_DU(A, n, k, escala_diag, contexto->u, y);
        break;
    }

    default:
        // Por default, sem PC
        vet_copy(n, r, y);
    }
}

/* ------------------------ Liberação do contexto ----------------------------- */

void pcg_free(pcg_contexto_t *contexto)
{
    free(contexto->y);
    free(contexto->D);
    free(contexto->invD);
    free(contexto->t);
    free(contexto->u);
    memset(contexto, 0, sizeof(*contexto));
}

/* ------------------------ Solver CG/PCG unificado --------------------------- */

// Baseado em CUNHA (2002)
int cg_solve(const matdiag_t *A, const real_t *b, real_t *x,
             int n, int k, int maxit, real_t eps_inf,
             pcg_precond_t M, real_t omega,
             real_t *norma_delta_x_inf_out)
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

    /* contexto do pré-condicionador (Jacobi / GS / SSOR) */
    pcg_contexto_t pc;
    memset(&pc, 0, sizeof(pc));

    if (M != PCG_PRECOND_NONE)
    {
        pcg_setup(A, n, M, omega, &pc);
    }

    /* r = b - A x inicial */
    matvet_diagonais(A, x, r);  /* r = A x */
    for (int i = 0; i < n; ++i) /* r = b - r */
        r[i] = b[i] - r[i];

    /* inicialização específica de CG puro ou PCG */
    real_t aux = 0.0;

    if (M == PCG_PRECOND_NONE)
    {
        aux = vet_produto(n, r, r); // CG puro: aux = r^T r
        vet_copy(n, r, v);          // v = r
    }
    else
    {
        pcg_apply(A, n, k, &pc, r, pc.y); // PCG: y = M⁻¹ r
        aux = vet_produto(n, r, pc.y);    // PCG: aux = r^T y
        vet_copy(n, pc.y, v);             // PCG: v = y
    }

    // Constante para proteger divisões por 0 ao calcular m = aux1/aux
    const real_t eps_den = 1e-30;

#ifdef LIKWID_PERFMON
    LIKWID_MARKER_START("op1");
#endif

    /* ------------------------------------------------------------------
     *  OP1: iteração do método de Gradientes Conjugados (CG/PCG)
     *
     *  Quando M = PCG_PRECOND_JACOBI (ω = 0), este laço implementa:
     *    - Aplicação do pré-condicionador Jacobi:   y = D⁻¹ r
     *    - Atualização da direção conjugada:        v = y + m v
     *    - Produto matriz-vetor esparso:            z = A v
     *    - Atualização da solução:                  x = x + s v
     *    - Atualização do resíduo:                  r = r − s z
     *
     *  É EXATAMENTE esta região que será medida como “op1”
     *  (tempo, banda de memória, cache miss, MFLOP/s) no Trabalho 2.
     * ------------------------------------------------------------------ */
    real_t dx_last = 0.0;

    for (int it = 0; it < maxit; ++it)
    {
        // [linha 3]  z = A v   (multiplicação matriz × vetor)
        matvet_diagonais(A, v, z);

        // [linha 4]  s = aux / (v^T z)
        real_t vTz = vet_produto(n, v, z);

        // robustez numérica: se denom ~0 e “energia” do resíduo já é pequena, finalize
        if (fabs(vTz) < 1e-30)
        {
            // sqrt(r^T r);  PCG: sqrt(r^T y) (aux)
            real_t res_energy = sqrt(fabs(aux));
            if (res_energy < 1e-12)
            {
                if (norma_delta_x_inf_out)
                    *norma_delta_x_inf_out = dx_last; // última Δx∞ conhecida
                pcg_free(&pc);
                free(r);
                free(v);
                free(z);
#ifdef LIKWID_PERFMON
                LIKWID_MARKER_STOP("op1");
#endif
                return it; // converge sem atualizar x nesta iteração
            }
            fprintf(stderr, "[cg/pcg] v^T z ~ 0 (divisao instavel)\n");
            break; // cai no final e libera tudo
        }

        real_t s = aux / vTz;

        // [linha 5]  x^(k+1) = x^(k) + s v
        // (e computa ||Δx||_∞ = max_i |s*v[i]|)
        real_t dx_max = 0.0;
        for (int i = 0; i < n; ++i)
        {
            real_t dx_i = s * v[i];
            x[i] += dx_i;
            real_t adx = fabs(dx_i);
            if (adx > dx_max)
                dx_max = adx;
        }
        dx_last = dx_max; // guarda a ÚLTIMA norma Δx∞ medida

        // [linha 6]  r = r − s z  (AXPY com alfa = -s)
        vet_axpy(n, -s, z, r);

        // ||Δx||_∞ < eps
        if (dx_max < eps_inf)
        {
            if (norma_delta_x_inf_out)
                *norma_delta_x_inf_out = dx_max; // devolve a Δx∞ final exata
            pcg_free(&pc);
            free(r);
            free(v);
            free(z);
#ifdef LIKWID_PERFMON
            LIKWID_MARKER_STOP("op1");
#endif
            return it + 1;
        }

        /* [linha 7]–[linha 11]
           Algoritmo 3.5 (CG “puro”) x PCG (Jacobi/GS/SSOR):
           - NONE:
             [7]  aux1 = r^T r
             [10] m = aux1/aux; aux = aux1
             [11] v = r + m v
           - PCG:
             [7]  y = M⁻¹ r
             [8]  aux1 = r^T y
             [10] m = aux1/aux; aux = aux1
             [11] v = y + m v
        */
        if (M == PCG_PRECOND_NONE)
        {
            real_t aux1 = vet_produto(n, r, r);
            real_t m = 0.0;

            if (fabs(aux) > eps_den)
                m = aux1 / aux;

            aux = aux1;
            for (int i = 0; i < n; ++i)
                v[i] = r[i] + m * v[i];
        }
        else
        {
            pcg_apply(A, n, k, &pc, r, pc.y); // y = M⁻¹ r
            real_t aux1 = vet_produto(n, r, pc.y);
            real_t m = 0.0;

            if (fabs(aux) > eps_den)
                m = aux1 / aux;

            aux = aux1;
            for (int i = 0; i < n; ++i)
                v[i] = pc.y[i] + m * v[i];
        }
    }
#ifdef LIKWID_PERFMON
    LIKWID_MARKER_STOP("op1");
#endif

    // atingiu limite de iterações sem satisfazer ||Δx||_∞
    if (norma_delta_x_inf_out)
        *norma_delta_x_inf_out = dx_last; // reporta a última norma medida

    pcg_free(&pc);
    free(r);
    free(v);
    free(z);

    return maxit;
}
