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
    // k deve ser 7!
    // k = 7;
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
        fprintf(stderr, "Erro: Valor de ω inválido (%.6g)\n", w);
        return 1;
    }

    srandom(20252);

    /* ----------------- Gera sistema k-diagonal esparso ----------------- */
    matdiag_t A;
    real_t *b = NULL;
    criaKDiagonal_v2(n, k, &A, &b);

    /* ----------------- Setup do pré-condicionador ----------------- */
    rtime_t t_pc_setup = 0.0;
    {
        pcg_contexto_t tmp = {0};
        rtime_t t0 = timestamp();
        int rc = pcg_setup(&A, n, k, M, omega, &tmp);
        t_pc_setup = timestamp() - t0;
        if (rc != 0)
        {
            fprintf(stderr, "Erro no setup do pre-condicionador (rc=%d)\n", rc);
            pcg_free(&tmp);
            liberaMatDiag(&A);
            free(b);
            return 1;
        }
        pcg_free(&tmp);
    }

    /* ----------------- Solver CG/PCG ----------------- */
    real_t *x = (real_t *)calloc((size_t)n, sizeof(real_t));
    if (!x)
    {
        fprintf(stderr, "Erro: memoria insuficiente para x\n");
        liberaMatDiag(&A);
        free(b);
        return 1;
    }

    rtime_t t_solve = 0.0;
    int iters = -1;
    real_t norma_delta_x_inf = NAN;

    {
        rtime_t t0 = timestamp();
        iters = cg_solve(&A, b, x, n, k, maxit, eps, M, omega, &norma_delta_x_inf);
        t_solve = timestamp() - t0;
    }

    if (iters < 0)
    {
        fprintf(stderr, "Falha no solver CG/PCG (codigo=%d)\n", iters);
        liberaMatDiag(&A);
        free(b);
        free(x);
        return 1;
    }

    /* ----------------- Cálculo do resíduo final ----------------- */
    rtime_t t_res = 0.0;
    real_t res_norm = 0.0;
    {
        rtime_t t0 = timestamp();
        // Recalcula r = b - A*x em formato esparso
        real_t *Ax = (real_t *)calloc(n, sizeof(real_t));
        matvet_diagonais(&A, x, Ax);
        for (int i = 0; i < n; ++i)
        {
            real_t ri = b[i] - Ax[i];
            res_norm += ri * ri;
        }
        res_norm = sqrt(res_norm);
        free(Ax);
        t_res = timestamp() - t0;
    }

    /* ----------------- Saída final ----------------- */
    printf("%d\n", n);
    for (int i = 0; i < n; ++i)
        printf("%.16g%s", x[i], (i + 1 < n ? " " : "\n"));
    printf("%.8g\n", norma_delta_x_inf);
    printf("%.8g\n", res_norm);
    printf("%.8g\n", t_pc_setup);
    printf("%.8g\n", (t_solve - t_pc_setup) / (rtime_t)iters);
    printf("%.8g\n", t_res);

    /* ----------------- Liberação de memória ----------------- */
    liberaMatDiag(&A);
    free(b);
    free(x);
    return 0;
}
