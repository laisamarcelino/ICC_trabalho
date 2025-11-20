#ifndef __SISLIN_H__
#define __SISLIN_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utils.h"

typedef struct
{
  int n;         // ordem do sistema
  int k;         // número de diagonais armazenadas
  int *offsets;  // deslocamento de cada diagonal em relação à principal (-p..+p)
  real_t **diag; // vetores com as diagonais (cada uma com n elementos)
} matdiag_t;

/* Ajustadas as assinaturas para usar real_t / rtime_t (consistente com sislin.c) */
void criaKDiagonal(int n, int k, real_t **A, real_t **B);
void criaKDiagonal_v2(int n, int k, matdiag_t *A, real_t **B);
void genSimetricaPositiva(real_t *A, real_t *b, int n, int k,
                          real_t **ASP, real_t *bsp, rtime_t *tempo);

#endif // __SISLIN_H__

