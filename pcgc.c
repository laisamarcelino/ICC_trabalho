#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "pcgc.h"
#include "utils.h"
#include "helpers.h"
// ---------- O que é o método dos Gradientes Conjugados (CG) ----------

// O objetivo do método (iterativo) dos Gradientes Conjugados é encontrar x, no qual temos Ax=b 
// ou seja, conhecemos o valor da matriz A (daremos) e o valor dos coeficientes b (daremos)

// Objetivo: queremos encontrar x que satisfaz a equação Ax=b.

// Antes de explicar o método, devemos saber que o método só funciona bem se a matriz A
// for simétrica definida positiva. Simétrica é quando ela é igual sua transposta (A = A^t).
// Definida positiva é quando o produto zᵗ*A*z > 0.

// z é um vetor qualquer (não nulo) que serve para testar o comportamento da matriz A.
// A operação zᵗ*(A*z) é um produto interno entre zᵗ e Az, por isso usamos a função auxiliar
// "dot", que representa o produto interno

// Nós precisamos desse teste (zᵗ*A*z > 0) pois ele garante que A nunca transforma um vetor (não nulo)
// em uma direção oposta ao próprio vetor, ou seja, que A é "positiva". Isso quer dizer que o algoritmo
// é numericamente estável e converge.

// Nós usamos esse zᵗ e z pois essa combinação transforma uma matriz em um ÚNICO número escalar.
// Isso permite ver com mais facilidade as propriedades de A.

// Vou continuar a explicação

// O método dos Gradientes Conjugados é iterativo. O gradiente (erro da iteração atual) dá a 
// DIREÇÃO em que o erro diminui mais rápido.

// Mas não podemos ir sempre pelo gradiente! Pois ele "faz zig-zag", não dando bons resultados sempre.

// Por isso, precisamos de direções independentes e "bem escolhidas", que são as direções conjugadas (pk)!

// Cada direção pk aponta para um novo caminho e corrige uma parte do erro. Ela é escolhida de modo
// que não interfere nas correções anteriores (iterações anteriores que vão convergindo para o método)

// ---------- CG sem pré-condicionador (ω = -1) ----------

int cg_no_prec(const real_t *A, const real_t *b, real_t *x,
               int n, int maxit, real_t eps,
               rtime_t *t_iter, rtime_t *t_res,
               real_t *res_norm_out)
{
  if (!A || !b || !x || n <= 0 || maxit <= 0 || eps < 0.0) {
    fprintf(stderr, "cg_no_prec: parâmetros inválidos\n");
    return -1;
  }

  // Alocação de vetores auxiliares
  real_t *r  = (real_t*) malloc((size_t)n * sizeof(real_t));
  real_t *p  = (real_t*) malloc((size_t)n * sizeof(real_t));
  real_t *Ap = (real_t*) malloc((size_t)n * sizeof(real_t));
  real_t *x_old = (real_t*) malloc((size_t)n * sizeof(real_t));
  if (!r || !p || !Ap || !x_old) {
    fprintf(stderr, "cg_no_prec: falha de alocação\n");
    free(r); free(p); free(Ap); free(x_old);
    return -1;
  }

  // Chute inicial do vetor (z, que foi dito acima)
  vec_set_zero(x, n);

  // Primeiro resíduo que será gerado
  vec_copy(r, b, n);

  // p0 = r0 (primeira direção da busca - pk, mas como é a primeira é p0 - sem pré-condicionador)

  // Sabemos que o pré-condicionador M é usado para modificar o resíduo: p0 = M⁻¹r0

  // A matriz M é uma aproximação de A que serve para facilitar sua solução. Quando
  // não temos pré condicionador, é equivalente a dizer que não queremos fazer nenhuma transformação, 
  // ou seja, temos M = I pois a matriz identidade não altera o resultado quando multiplicamos ela por
  // um vetor, logo, p0 = I⁻¹r0 -> p0 = r0

  // Lembrando que usando o métodos CG com pré-condicionadores, quer dizer multiplicar uma matriz
  // M (auxiliar) a um sistema, pra tornar ele mais fácil de resolver. Com a notação M⁻¹Ax = M⁻¹b.
  
  // Porém, ao invés de aplicarmos diretamente M⁻¹Ax = M⁻¹b, a teoria nos diz que precisamos aplicar
  // isso nas etapas de atualização do método, ou seja, nos cálculos de resíduo (rk) e direção (pk)
  // que é aqui mesmo.
  vec_copy(p, r, n);

  // Calcula o produto interno do resíduo com ele mesmo.
  // Isso mete o tamanho do erro, e indica o quão longe (passo) está a solução atual (xk) da exata
  // Sem esse rho, não conseguimos criar as outras direções conjugadas
  real_t rho = dot(r, r, n);
  // Se rho for muito pequeno, x0 já satisfaz Ax=b então método nem precisa iterar
  if (rho <= DBL_EPSILON) {
    if (t_iter) 
      *t_iter = 0.0;

    if (t_res) {
      *t_res = timestamp();
      *t_res = timestamp() - *t_res;
    }
    if (res_norm_out) 
      *res_norm_out = sqrt(rho);
    free(r); free(p); free(Ap); free(x_old);
    return 0;
  }

  int it = 0;
  rtime_t t_loop = timestamp();

  // Caso contrário, temos o início das iterações
  // Aqui é o início das iterações do método CG sem pré-condicionadores
  for (it = 1; it <= maxit; ++it) {
    // Vamos calcular o tamanho do passo na direção de busca pk
    // Isso é bom para mover x na direção de busca pk para reduzir o resíduo
    // matvec_dense é a multiplicação de uma matriz por um vetor
    // Estou multiplicando a matriz A pelo vetor pk
    // O resultado é armazenado no vetor Ap
    matvec_dense(A, p, Ap, n);

    // Aqui estamos medindo a curvatura da função de busca na direção pk
    // Se o "denominador" for muito peuqueno ou negativo, o passo ak não faria sentido
    real_t denom = dot(p, Ap, n); // p^T A p
    if (!(denom > 0.0) || !isfinite(denom)) {
      fprintf(stderr, "cg_no_prec: breakdown numérico (p^T A p = %g)\n", denom);
      free(r); free(p); free(Ap); free(x_old);
      return -1;
    }
    // Com esse valor, calculamos o tamanho do passo
    // Serve para sabermos a quantidade que devemos avançar na direção de pk
    // ak é o passo e nos diz o quanto devemos andar
    real_t alpha = rho / denom;

    // x_{k+1} = x_k + alpha p_k
    // Salvamos o antigo x (para medir a diferença)
    vec_copy(x_old, x, n);
    // Aplicamos o passo na direção de pk e geramos um novo x
    axpy(x, alpha, p, n);
    // Após isso, nós conseguimos avançar na iteração do método CG.
    // Caminhamos na "direção de busca ótima"

    // r_{k+1} = r_k - alpha Ap
    // Aqui, atualizamos o resíduo (diminuindo ele, até que fique pequeno o suficiente para parar)
    // Resíduo = resíduo - passo * Ap
    for (int i = 0; i < n; ++i) 
      r[i] -= alpha * Ap[i];

    // critério do trabalho: ||x_new - x_old||_inf < eps
    // Aqui está o criério de parada da iteração, que decide se o método já convergiu para a solução
    // O cálculo dele está na função norm_inf_diff e vamos comparar ela com a tolerância que queremos (eps)
    real_t aprox_err = norm_inf_diff(x, x_old, n);
    if (aprox_err < eps) 
      break;

    // Nós calculamos o novo valor de rho depois de atualizar o resíduo
    // Isso serve para medir o novo tamanho do erro e usado para calcular a próxima direção de busca (β_k)
    // β_k é um escalar que decide quanto da direção antiga pk deve entrar na nova direção
    real_t rho_new = dot(r, r, n); // (r_{k+1})^T (r_{k+1})
    if (!isfinite(rho_new)) {
      fprintf(stderr, "cg_no_prec: rho_new inválido\n");
      free(r); free(p); free(Ap); free(x_old);
      return -1;
    }

    // Cálcula o escalar da nova direção de busca (β)
    real_t beta = rho_new / rho;

    // A nova direção é a combinação do novo resíduo com a direção de busca anterior pk * β_k
    // p_{k+1} = r_{k+1} + beta p_k
    for (int i = 0; i < n; ++i) 
      p[i] = r[i] + beta * p[i];

    // Atualiza de fato rho
    rho = rho_new;
  }

  // Cálculo do tempo métdio por iteração
  t_loop = timestamp() - t_loop;
  if (t_iter) {
    if (it > 0)
      *t_iter = t_loop / (real_t)it;
    else
      *t_iter = 0.0;
  }

  // Cálculo do resíduo final e tempo do resíduo final
  if (t_res || res_norm_out) {
    if (t_res) *t_res = timestamp();

    // Calcula o resíduo final
    // r = b - A x
    matvec_dense(A, x, Ap, n);
    for (int i = 0; i < n; ++i) 
      r[i] = b[i] - Ap[i];
    // Calcula a norma L2 do resíduo
    if (res_norm_out) 
      *res_norm_out = norm2(r, n);
    // Calcula o tempo gasto para calcular o resíduo final
    if (t_res) 
      *t_res = timestamp() - *t_res;
  }

  // Liberação de memória
  free(r); free(p); free(Ap); free(x_old);
  // Retorna o número de iterações realizadas
  return it;
}

// ---------- CG com pré-condicionador Jacobi (ω = 0.0) ----------

// Função auxiliar para resolver M z = r, onde M = D (diagonal de A)
// M é o pré-condicionador (M = D)
// z é o resíduo pré-condicionado atual (D⁻¹ × r_atual)
// r é o resíduo atual: r = b - A*x

/*
Exemplo:
M[0,0] * z[0] = r[0]  →  z[0] = r[0] / M[0,0] = r[0] / A[0,0]
M[1,1] * z[1] = r[1]  →  z[1] = r[1] / M[1,1] = r[1] / A[1,1]  
M[2,2] * z[2] = r[2]  →  z[2] = r[2] / M[2,2] = r[2] / A[2,2]
*/

// Essa função (jacobi_solve) faz a aplicação do pré-condicionador de Jacobi
// Ela pega o elemento da diagonal principal da Matriz A (Aii)
// Se Aii for muito pequeno, então define z[i] = 0 para evitar divisão por zero
// Caso contraŕio, faz o cálculo de z (z = r / Aii)
// Lembrando que z é um vetor qualquer para testar o comportamento da matriz A
// Também lembrando que o resíduo é r = b - Ax

// O pré-condicionador de Jacobi , pela teoria, usa apenas a diagonal de A, por isso é dessa maneira
// Essa divisão (z = r / Aii) serve para corrigir (ou normalizar) cada componente do resíduo de acordo com a escala da diagonal de A
// Cada z[i] terá os valores normalizados da diagonal principal (Aii) da matriz A
void jacobi_solve(const real_t *A, const real_t *r, real_t *z, int n) {
  for (int i = 0; i < n; ++i) {
    real_t aii = A[i*(size_t)n + i]; // elemento diagonal A[i][i]
    if (fabs(aii) <= DBL_EPSILON) {
      fprintf(stderr, "jacobi_solve: elemento diagonal A[%d][%d] = %g muito pequeno\n", i, i, aii);
      z[i] = 0.0; // Se for muito próximo de zero, é igual a zero
    } else {
      z[i] = r[i] / aii; // resíduo pré-condicionado = resíduo atual / A[i][i]
    }
  }
}

int cg_jacobi_prec(const real_t *A, const real_t *b, real_t *x,
                   int n, int maxit, real_t eps,
                   rtime_t *t_iter, rtime_t *t_res,
                   real_t *res_norm_out)
{
  if (!A || !b || !x || n <= 0 || maxit <= 0 || eps < 0.0) {
    fprintf(stderr, "cg_jacobi_prec: parâmetros inválidos\n");
    return -1;
  }

  // aloca vetores auxiliares
  real_t *r  = (real_t*) malloc((size_t)n * sizeof(real_t));
  real_t *z  = (real_t*) malloc((size_t)n * sizeof(real_t));
  real_t *p  = (real_t*) malloc((size_t)n * sizeof(real_t));
  real_t *Ap = (real_t*) malloc((size_t)n * sizeof(real_t));
  real_t *x_old = (real_t*) malloc((size_t)n * sizeof(real_t));
  if (!r || !z || !p || !Ap || !x_old) {
    fprintf(stderr, "cg_jacobi_prec: falha de alocação\n");
    free(r); free(z); free(p); free(Ap); free(x_old);
    return -1;
  }

  // O resto é exatamente igual a função sem pré-condicionador, mas agora temos z com 
  // valores normalizados (corrigidos) pelo pré-condicionador
  // x0 = 0
  vec_set_zero(x, n);

  // r0 = b - A*x0 = b
  vec_copy(r, b, n);

  // z0 = M^{-1} * r0 (resolve M z = r, onde M = D)
  jacobi_solve(A, r, z, n);

  // p0 = z0 (direção inicial)
  vec_copy(p, z, n);

  real_t rho = dot(r, z, n); // r^T z (produto interno modificado)
  if (rho <= DBL_EPSILON) {
    // já convergiu em x0
    if (t_iter) *t_iter = 0.0;
    // resíduo final (≈ ||b||2, pois x=0)
    if (t_res) {
      *t_res = timestamp();
      *t_res = timestamp() - *t_res;
    }
    if (res_norm_out) *res_norm_out = norm2(r, n);
    free(r); free(z); free(p); free(Ap); free(x_old);
    return 0;
  }

  int it = 0;
  rtime_t t_loop = timestamp();

  for (it = 1; it <= maxit; ++it) {
    // Ap = A * p
    matvec_dense(A, p, Ap, n);

    real_t denom = dot(p, Ap, n); // p^T A p
    if (!(denom > 0.0) || !isfinite(denom)) {
      fprintf(stderr, "cg_jacobi_prec: breakdown numérico (p^T A p = %g)\n", denom);
      free(r); free(z); free(p); free(Ap); free(x_old);
      return -1;
    }

    real_t alpha = rho / denom;

    // x_{k+1} = x_k + alpha p_k
    vec_copy(x_old, x, n);
    axpy(x, alpha, p, n);

    // r_{k+1} = r_k - alpha Ap
    for (int i = 0; i < n; ++i) r[i] -= alpha * Ap[i];

    // critério do trabalho: ||x_new - x_old||_inf < eps
    real_t aprox_err = norm_inf_diff(x, x_old, n);
    if (aprox_err < eps) break;

    // z_{k+1} = M^{-1} * r_{k+1}
    jacobi_solve(A, r, z, n);

    real_t rho_new = dot(r, z, n); // r_{k+1}^T z_{k+1}
    if (!isfinite(rho_new)) {
      fprintf(stderr, "cg_jacobi_prec: rho_new inválido\n");
      free(r); free(z); free(p); free(Ap); free(x_old);
      return -1;
    }

    real_t beta = rho_new / rho;

    // p_{k+1} = z_{k+1} + beta p_k
    for (int i = 0; i < n; ++i) p[i] = z[i] + beta * p[i];

    rho = rho_new;
  }

  t_loop = timestamp() - t_loop;
  if (t_iter) *t_iter = (it > 0 ? t_loop / (real_t)it : 0.0);

  // resíduo final e tempo_residuo
  if (t_res || res_norm_out) {
    if (t_res) *t_res = timestamp();

    // r = b - A x
    matvec_dense(A, x, Ap, n); // reutiliza Ap como A*x
    for (int i = 0; i < n; ++i) r[i] = b[i] - Ap[i];

    if (res_norm_out) *res_norm_out = norm2(r, n);
    if (t_res) *t_res = timestamp() - *t_res;
  }

  free(r); free(z); free(p); free(Ap); free(x_old);
  return it;
}

// ---------- CG com pré-condicionador de Gauss-Seidel (ω = 1.0) ----------

// Essa função (gauss_seidel_solve) faz a aplicação do pré-condicionador de Gauss-Seidel
// ... explicarei depois que entender ...
// ...

// Função auxiliar: resolve (D+L)z = r por forward substitution
void gauss_seidel_solve(const real_t *A, const real_t *r, real_t *z, int n) {
  for (int i = 0; i < n; ++i) {
    real_t sum = 0.0;
    for (int j = 0; j < i; ++j) {
      sum += A[i*(size_t)n + j] * z[j];
    }
    real_t dii = A[i*(size_t)n + i];
    if (fabs(dii) <= DBL_EPSILON) {
      fprintf(stderr, "gauss_seidel_solve: elemento diagonal A[%d][%d] = %g muito pequeno\n", i, i, dii);
      z[i] = 0.0;
    } else {
      z[i] = (r[i] - sum) / dii;
    }
  }
}

int cg_gs_prec(const real_t *A, const real_t *b, real_t *x,
               int n, int maxit, real_t eps,
               rtime_t *t_iter, rtime_t *t_res,
               real_t *res_norm_out)
{
  if (!A || !b || !x || n <= 0 || maxit <= 0 || eps < 0.0) {
    fprintf(stderr, "cg_gs_prec: parâmetros inválidos\n");
    return -1;
  }

  real_t *r  = (real_t*) malloc((size_t)n * sizeof(real_t));
  real_t *z  = (real_t*) malloc((size_t)n * sizeof(real_t));
  real_t *p  = (real_t*) malloc((size_t)n * sizeof(real_t));
  real_t *Ap = (real_t*) malloc((size_t)n * sizeof(real_t));
  real_t *x_old = (real_t*) malloc((size_t)n * sizeof(real_t));
  if (!r || !z || !p || !Ap || !x_old) {
    fprintf(stderr, "cg_gs_prec: falha de alocação\n");
    free(r); free(z); free(p); free(Ap); free(x_old);
    return -1;
  }

  vec_set_zero(x, n);
  vec_copy(r, b, n);

  gauss_seidel_solve(A, r, z, n);
  vec_copy(p, z, n);

  real_t rho = dot(r, z, n);
  if (rho <= DBL_EPSILON) {
    if (t_iter) *t_iter = 0.0;
    if (t_res) {
      *t_res = timestamp();
      *t_res = timestamp() - *t_res;
    }
    if (res_norm_out) *res_norm_out = norm2(r, n);
    free(r); free(z); free(p); free(Ap); free(x_old);
    return 0;
  }

  int it = 0;
  rtime_t t_loop = timestamp();

  for (it = 1; it <= maxit; ++it) {
    matvec_dense(A, p, Ap, n);

    real_t denom = dot(p, Ap, n);
    if (!(denom > 0.0) || !isfinite(denom)) {
      fprintf(stderr, "cg_gs_prec: breakdown numérico (p^T A p = %g)\n", denom);
      free(r); free(z); free(p); free(Ap); free(x_old);
      return -1;
    }

    real_t alpha = rho / denom;

    vec_copy(x_old, x, n);
    axpy(x, alpha, p, n);

    for (int i = 0; i < n; ++i) r[i] -= alpha * Ap[i];

    real_t aprox_err = norm_inf_diff(x, x_old, n);
    if (aprox_err < eps) break;

    gauss_seidel_solve(A, r, z, n);

    real_t rho_new = dot(r, z, n);
    if (!isfinite(rho_new)) {
      fprintf(stderr, "cg_gs_prec: rho_new inválido\n");
      free(r); free(z); free(p); free(Ap); free(x_old);
      return -1;
    }

    real_t beta = rho_new / rho;
    for (int i = 0; i < n; ++i) p[i] = z[i] + beta * p[i];

    rho = rho_new;
  }

  t_loop = timestamp() - t_loop;
  if (t_iter) *t_iter = (it > 0 ? t_loop / (real_t)it : 0.0);

  if (t_res || res_norm_out) {
    if (t_res) *t_res = timestamp();
    matvec_dense(A, x, Ap, n);
    for (int i = 0; i < n; ++i) r[i] = b[i] - Ap[i];
    if (res_norm_out) *res_norm_out = norm2(r, n);
    if (t_res) *t_res = timestamp() - *t_res;
  }

  free(r); free(z); free(p); free(Ap); free(x_old);
  return it;
}

// ---------- CG com pré-condicionador SSOR (1.0 < ω < 2.0) ----------

// Essa função (ssor_prec_solve) faz a aplicação do pré-condicionador SSOR
// ... explicarei depois que entender ...
// ...

// Função auxiliar: aplica o pré-condicionador SSOR (z = M^{-1} r)
void ssor_prec_solve(const real_t *A, const real_t *r, real_t *z, int n, real_t omega) {
  // Monta D, L, U
  real_t *D = (real_t*) calloc((size_t)n * n, sizeof(real_t));
  real_t *L = (real_t*) calloc((size_t)n * n, sizeof(real_t));
  real_t *U = (real_t*) calloc((size_t)n * n, sizeof(real_t));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      real_t aij = A[i*n + j];
      if (i == j) D[i*n + j] = aij;
      else if (j < i) L[i*n + j] = aij;
      else U[i*n + j] = aij;
    }
  }
  // Monta (D + omega*L) e (D + omega*U)
  real_t *DL = (real_t*) calloc((size_t)n * n, sizeof(real_t));
  real_t *DU = (real_t*) calloc((size_t)n * n, sizeof(real_t));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      DL[i*n + j] = (i == j ? D[i*n + i] : 0.0) + omega * L[i*n + j];
      DU[i*n + j] = (i == j ? D[i*n + i] : 0.0) + omega * U[i*n + j];
    }
  }
  // Forward: (D + omega*L) y = r
  real_t *y = (real_t*) calloc((size_t)n, sizeof(real_t));
  for (int i = 0; i < n; ++i) {
    real_t sum = 0.0;
    for (int j = 0; j < i; ++j)
      sum += DL[i*n + j] * y[j];
    y[i] = (r[i] - sum) / DL[i*n + i];
  }
  // D y
  real_t *Dy = (real_t*) calloc((size_t)n, sizeof(real_t));
  for (int i = 0; i < n; ++i)
    Dy[i] = D[i*n + i] * y[i];
  // Backward: (D + omega*U) z = D y
  for (int i = n-1; i >= 0; --i) {
    real_t sum = 0.0;
    for (int j = i+1; j < n; ++j)
      sum += DU[i*n + j] * z[j];
    z[i] = (Dy[i] - sum) / DU[i*n + i];
  }
  free(D); free(L); free(U); free(DL); free(DU); free(y); free(Dy);
}

int cg_ssor_prec(const real_t *A, const real_t *b, real_t *x,
                 int n, int maxit, real_t eps, real_t omega,
                 rtime_t *t_iter, rtime_t *t_res,
                 real_t *res_norm_out)
{
  if (!A || !b || !x || n <= 0 || maxit <= 0 || eps < 0.0 || omega <= 1.0 || omega >= 2.0) {
    fprintf(stderr, "cg_ssor_prec: parâmetros inválidos\n");
    return -1;
  }

  real_t *r  = (real_t*) malloc((size_t)n * sizeof(real_t));
  real_t *z  = (real_t*) malloc((size_t)n * sizeof(real_t));
  real_t *p  = (real_t*) malloc((size_t)n * sizeof(real_t));
  real_t *Ap = (real_t*) malloc((size_t)n * sizeof(real_t));
  real_t *x_old = (real_t*) malloc((size_t)n * sizeof(real_t));
  if (!r || !z || !p || !Ap || !x_old) {
    fprintf(stderr, "cg_ssor_prec: falha de alocação\n");
    free(r); free(z); free(p); free(Ap); free(x_old);
    return -1;
  }

  vec_set_zero(x, n);
  vec_copy(r, b, n);

  ssor_prec_solve(A, r, z, n, omega);
  vec_copy(p, z, n);

  real_t rho = dot(r, z, n);
  if (rho <= DBL_EPSILON) {
    if (t_iter) *t_iter = 0.0;
    if (t_res) {
      *t_res = timestamp();
      *t_res = timestamp() - *t_res;
    }
    if (res_norm_out) *res_norm_out = norm2(r, n);
    free(r); free(z); free(p); free(Ap); free(x_old);
    return 0;
  }

  int it = 0;
  rtime_t t_loop = timestamp();

  for (it = 1; it <= maxit; ++it) {
    matvec_dense(A, p, Ap, n);

    real_t denom = dot(p, Ap, n);
    if (!(denom > 0.0) || !isfinite(denom)) {
      fprintf(stderr, "cg_ssor_prec: breakdown numérico (p^T A p = %g)\n", denom);
      free(r); free(z); free(p); free(Ap); free(x_old);
      return -1;
    }

    real_t alpha = rho / denom;

    vec_copy(x_old, x, n);
    axpy(x, alpha, p, n);

    for (int i = 0; i < n; ++i) r[i] -= alpha * Ap[i];

    real_t aprox_err = norm_inf_diff(x, x_old, n);
    if (aprox_err < eps) break;

    ssor_prec_solve(A, r, z, n, omega);

    real_t rho_new = dot(r, z, n);
    if (!isfinite(rho_new)) {
      fprintf(stderr, "cg_ssor_prec: rho_new inválido\n");
      free(r); free(z); free(p); free(Ap); free(x_old);
      return -1;
    }

    real_t beta = rho_new / rho;
    for (int i = 0; i < n; ++i) p[i] = z[i] + beta * p[i];

    rho = rho_new;
  }

  t_loop = timestamp() - t_loop;
  if (t_iter) *t_iter = (it > 0 ? t_loop / (real_t)it : 0.0);

  if (t_res || res_norm_out) {
    if (t_res) *t_res = timestamp();
    matvec_dense(A, x, Ap, n);
    for (int i = 0; i < n; ++i) r[i] = b[i] - Ap[i];
    if (res_norm_out) *res_norm_out = norm2(r, n);
    if (t_res) *t_res = timestamp() - *t_res;
  }

  free(r); free(z); free(p); free(Ap); free(x_old);
  return it;
}