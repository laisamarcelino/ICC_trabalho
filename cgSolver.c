#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sislin.h"
#include "pcgc.h"
#include "utils.h"
#include "helpers.h"

// Mapeia ω lido em M (tipo de PCG) e normaliza omega_out se necessário
static int escolhe_precond(real_t w_in, pcg_precond_t *M_out, real_t *omega_out)
{
    if (w_in == -1.0)
    { // sem pré-condicionador 
        *M_out = PCG_PRECOND_NONE;
        *omega_out = 0.0;
        return 0;
    }
    if (w_in == 0.0)
    { // Jacobi 
        *M_out = PCG_PRECOND_JACOBI;
        *omega_out = 0.0;
        return 0;
    }
    if (w_in == 1.0)
    { // Gauss-Seidel 
        *M_out = PCG_PRECOND_SGS;
        *omega_out = 1.0;
        return 0;
    }
    if (w_in > 1.0 && w_in < 2.0)
    { // SSOR 
        *M_out = PCG_PRECOND_SSOR;
        *omega_out = w_in; // manter ω (0<ω<2)
        return 0;
    }
    return -1; // inválido
}

int main(void)
{
    int n, k, maxit;
    real_t w, eps;

    /* ----------------- Leitura e validação dos parâmetros ----------------- */
    if (scanf("%d", &n) != 1 || n <= 10)
    {
        fprintf(stderr, "Erro: n deve ser > 10\n");
        return 1;
    }
    if (scanf("%d", &k) != 1 || k <= 1 || (k % 2 == 0))
    {
        fprintf(stderr, "Erro: k deve ser > 1 e ímpar\n");
        return 1;
    }
    if (scanf("%lf", &w) != 1)
    {
        fprintf(stderr, "Erro ao ler ω\n");
        return 1;
    }
    if (scanf("%d", &maxit) != 1 || maxit <= 0)
    {
        fprintf(stderr, "Erro: maxit deve ser > 0\n");
        return 1;
    }
    if (scanf("%lf", &eps) != 1 || eps <= 0.0)
    {
        fprintf(stderr, "Erro: ε deve ser > 0\n");
        return 1;
    }

    /* ----------------- Escolha do pré-condicionador ----------------- */
    pcg_precond_t M = PCG_PRECOND_NONE;
    real_t omega = 0.0;
    if (escolhe_precond(w, &M, &omega) != 0)
    {
        fprintf(stderr, "Erro: Valor de ω inválido (%.6g). Valores aceitos:\n", w);
        fprintf(stderr, "  ω = -1.0 : sem pré-condicionador (NONE)\n");
        fprintf(stderr, "  ω =  0.0 : Jacobi\n");
        fprintf(stderr, "  ω =  1.0 : Gauss-Seidel (SGS)\n");
        fprintf(stderr, "  1.0 < ω < 2.0 : SSOR\n");
        return 1;
    }

    /* Semente do gerador aleatório (conforme enunciado) */
    srandom(20252);

    /* ----------------- Gerar sistema k-diagonal A, b ----------------- */
    real_t *A = NULL, *b = NULL;
    criaKDiagonal(n, k, &A, &b);

    /* ----------------- Transformar para SPD: ASP = AᵗA, bsp = Aᵗ b ----------------- */
    real_t *ASP = NULL;
    real_t *bsp = (real_t *)calloc((size_t)n, sizeof(real_t));
    if (!bsp)
    {
        fprintf(stderr, "Erro: memoria insuficiente para bsp\n");
        free(A);
        free(b);
        return 1;
    }
    rtime_t t_spd = 0.0;
    genSimetricaPositiva(A, b, n, k, &ASP, bsp, &t_spd);

    /* ----------------- Medir tempo do setup do pré-condicionador ----------------- */
    rtime_t t_pc_setup = 0.0;
    {
        pcg_contexto_t tmp = {0};
        rtime_t t0 = timestamp();
        int rc = pcg_setup(ASP, n, k, M, omega, &tmp);
        t_pc_setup = timestamp() - t0;
        if (rc != 0)
        {
            fprintf(stderr, "Erro no setup do pre-condicionador (rc=%d)\n", rc);
            pcg_free(&tmp);
            free(A);
            free(b);
            free(ASP);
            free(bsp);
            return 1;
        }
        pcg_free(&tmp);
    }

    /* ----------------- Resolver por CG/PCG (mede tempo total solve) ----------------- */
    real_t *x = (real_t *)calloc((size_t)n, sizeof(real_t));
    if (!x)
    {
        fprintf(stderr, "Erro: memoria insuficiente para x\n");
        free(A);
        free(b);
        free(ASP);
        free(bsp);
        return 1;
    }

    rtime_t t_solve = 0.0;
    int iters = -1;
    real_t norma_delta_x_inf = NAN; // Recebe ||Δx||_∞ da última iteração
    {
        rtime_t t0 = timestamp();
        iters = cg_solve(ASP, bsp, x, n, k, maxit, eps, M, omega,
                         &norma_delta_x_inf);
        t_solve = timestamp() - t0;
    }

    if (iters < 0)
    {
        fprintf(stderr, "Falha no solver CG/PCG (codigo=%d)\n", iters);
        free(A);
        free(b);
        free(ASP);
        free(bsp);
        free(x);
        return 1;
    }

    /* ----------------- Calcular resíduo final e seus tempos ---------- */
    rtime_t t_res = 0.0;
    real_t res_norm = 0.0;
    {
        rtime_t t0 = timestamp();
        res_norm = residuo_l2(ASP, bsp, x, n); /* ||b - A x||2 */
        t_res = timestamp() - t0;
    }

    /* ----------------- Computar tempos exigidos na saída ------------- */
    rtime_t tempo_pc = t_spd + t_pc_setup;
    rtime_t tempo_iter = 0.0;
    if (iters > 0)
    {
        rtime_t iter_total_aprox = t_solve - t_pc_setup;
        if (iter_total_aprox < 0.0)
            iter_total_aprox = 0.0;
        tempo_iter = iter_total_aprox / (rtime_t)iters;
    }

    /* ----------------- Saída no formato do enunciado ----------------- */
    // Linha 1: n
    printf("%d\n", n);

    // Linha 2: x_1 x_2 ... x_n (com precisão pedida)
    for (int i = 0; i < n; ++i)
        printf("%.16g%s", x[i], (i + 1 < n ? " " : "\n"));

    // Linha 3: norma (||Δx||_∞ na última iteração) — agora valor EXATO devolvido pelo solver
    printf("%.8g\n", norma_delta_x_inf);

    // Linha 4: residuo = ||r||2 = ||b - A x||2
    printf("%.8g\n", res_norm);

    // Linha 5: tempo_pc (ms) = preparar SL + setup PC
    printf("%.8g\n", tempo_pc);

    // Linha 6: tempo_iter (ms) = tempo médio por iteração (aprox)
    printf("%.8g\n", tempo_iter);

    // Linha 7: tempo_residuo (ms)
    printf("%.8g\n", t_res);

    // Liberação de memória
    free(A);
    free(b);
    free(ASP);
    free(bsp);
    free(x);
    return 0;
}
