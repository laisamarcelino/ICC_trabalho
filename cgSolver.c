#include <stdio.h>
#include "sislin.h"
#include "utils.h"
#include "pcgc.h"

int main()
{

    /*-------- teste criaKDiagonal e genSimetricaPositiva -----------*/
    int n = 5, k = 3;
    real_t *A = calloc(n * n, sizeof(real_t));
    real_t *b = calloc(n, sizeof(real_t));
    real_t *ASP = NULL;
    real_t *bsp = calloc(n, sizeof(real_t));
    rtime_t t;

    criaKDiagonal(n, k, &A, &b);
    genSimetricaPositiva(A, b, n, k, &ASP, bsp, &t);

    // Imprime matriz ASP
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%.3f ", ASP[i * n + j]);
        }
        printf("\n");
    }

    return 0;
}