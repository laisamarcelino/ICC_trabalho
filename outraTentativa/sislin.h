#ifndef __SISLIN_H__
#define __SISLIN_H__

#include <stdio.h>
#include <stdlib.h> /* for exit e random/srandom */
#include <string.h>
#include <math.h>

#include "utils.h"

void criaKDiagonal(int n, int k, real_t **A, real_t **B);

void genSimetricaPositiva(real_t *A, real_t *b, int n, int k, real_t **ASP, real_t *bsp, real_t *tempo);
void geraDLU (real_t *A, int n, int k, real_t **D, real_t **L, real_t **U, real_t *tempo);
void geraPreCond(real_t *D, real_t *L, real_t *U, real_t w, int n, int k, real_t **M, real_t *tempo);
real_t calcResiduoSL (real_t *A, real_t *b, real_t *X, int n, int k, real_t *tempo);

#endif // __SISLIN_H__

