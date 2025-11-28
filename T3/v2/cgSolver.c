#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef LIKWID_PERFMON
#include <likwid.h>
#endif

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
#ifdef LIKWID_PERFMON
    LIKWID_MARKER_INIT;
#endif
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
    if (scanf("%d", &k) != 1 || k != 7)
    {
        fprintf(stderr, "Erro: k deve ser 7\n");
        return 1;
    }
    if (scanf("%lf", &w) != 1)
    {
        fprintf(stderr, "Erro ao ler ω\n");
        return 1;
    }
    if (scanf("%d", &maxit) != 1 || maxit <= 0 || maxit > 25)
    {
        fprintf(stderr, "Erro: maxit deve estar entre 0 e 25\n");
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
    matdiag_t A = {0};
    real_t *b = NULL;
    criaKDiagonal_v2(n, k, &A, &b);

    /* ----------------- Transformar para SPD ----------------- */
    matdiag_t ASP = {0};
    real_t *bsp = NULL;
    rtime_t t_spd = 0.0;
    genSimetricaPositiva_diag(&A, b, &ASP, &bsp, &t_spd);
    liberaMatDiag(&A);
    free(b);

    /* ----------------- Setup do pré-condicionador ----------------- */
    rtime_t t_pc_setup = 0.0;
    {
        pcg_contexto_t tmp = {0};
        rtime_t t0 = timestamp();
        int rc = pcg_setup(&ASP, n, M, omega, &tmp);
        t_pc_setup = timestamp() - t0;
        if (rc != 0)
        {
            fprintf(stderr, "Erro no setup do pre-condicionador (rc=%d)\n", rc);
            pcg_free(&tmp);
            liberaMatDiag(&ASP);
            free(bsp);
            return 1;
        }
        pcg_free(&tmp);
    }

    /* ----------------- Solver CG/PCG ----------------- */
    real_t *x = (real_t *)calloc((size_t)n, sizeof(real_t));
    if (!x)
    {
        fprintf(stderr, "Erro: memoria insuficiente para x\n");
        liberaMatDiag(&ASP);
        free(bsp);
        return 1;
    }

    // --- NOVO: medir tempo de cada iteração (op1) ---
    rtime_t t_op1_total = 0.0;
    int iters = -1;
    real_t norma_delta_x_inf = NAN;

    {
        // Setup para medir tempo de cada iteração
        rtime_t t0 = timestamp();
        // Chamar cg_solve modificado para rodar exatamente maxit iterações
        // e medir tempo total das iterações (op1)
        iters = cg_solve(&ASP, bsp, x, n, ASP.k, maxit, 0.0, M, omega, &norma_delta_x_inf);
        t_op1_total = timestamp() - t0;
    }

    if (iters < 0)
    {
        fprintf(stderr, "Falha no solver CG/PCG (codigo=%d)\n", iters);
        liberaMatDiag(&ASP);
        free(bsp);
        free(x);
        return 1;
    }

    /* ----------------- Cálculo do resíduo final (op2) ----------------- */
    rtime_t t_op2 = 0.0;
    real_t res_norm = 0.0;
    {
        rtime_t t0 = timestamp();
#ifdef LIKWID_PERFMON
        LIKWID_MARKER_START("op2");
#endif
        res_norm = residuo_l2_v2(&ASP, bsp, x);
#ifdef LIKWID_PERFMON
        LIKWID_MARKER_STOP("op2");
#endif
        t_op2 = timestamp() - t0;
    }

    // --- NOVO: imprimir tempos médios das operações ---
    rtime_t tempo_op1_medio = (iters > 0) ? (t_op1_total / (rtime_t)iters) : 0.0;

    printf("N = %d\n", n);
    printf("Tempo médio op1 (CG iteração): %.8g ms\n", tempo_op1_medio);
    printf("Tempo op2 (resíduo): %.8g ms\n", t_op2);
    
    /* ----------------- Cálculo dos tempos de saída ----------------- */
    rtime_t tempo_iter = 0.0;
    if (iters > 0)
    {
        tempo_iter = t_op1_total / (rtime_t)iters; /* [FIX-02] tempo médio usa apenas o solve */
    }
    else
    {
        tempo_iter = 0.0;
    }

    /* ----------------- Saída final -----------------
     * Formato:
     *   n
     *   x[0] ... x[n-1]
     *   ||Δx||_∞ (última iteração)
     *   ||r||_2  (OP2)
     *   t_spd + t_pc_setup     (geração SPD + setup pré-condicionador)
     *   tempo_iter             (tempo médio por iteração da OP1)
     *   t_res                  (tempo de OP2: residuo_l2_v2)
     * ------------------------------------------------ */
    printf("%d\n", n);
    for (int i = 0; i < n; ++i)
        printf("%.16g%s", x[i], (i + 1 < n ? " " : "\n"));
    printf("%.8g\n", norma_delta_x_inf);
    printf("%.8g\n", res_norm);
    printf("%.8g\n", t_spd + t_pc_setup);
    printf("%.8g\n", tempo_iter); /* [FIX-02] usa valor calculado e válido */
    printf("%.8g\n", t_op2);

    /* ----------------- Liberação de memória ----------------- */
    liberaMatDiag(&ASP);
    free(bsp);
    free(x);
#ifdef LIKWID_PERFMON
    LIKWID_MARKER_CLOSE;
#endif
    return 0;
}
