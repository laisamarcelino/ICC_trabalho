#ifndef __SISLIN_H__
#define __SISLIN_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"
#include "helpers.h"

void criaKDiagonal(int n, int k, double **A, double **B);
void genSimetricaPositiva(double *A, double *b, int n, int k, double **ASP, double *bsp, double *tempo);

#endif // __SISLIN_H__

