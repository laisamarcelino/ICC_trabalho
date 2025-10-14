# Trabalho 01 - Solução de Sistemas Lineares Esparsos com Pré-condicionadores

## Autoria

- Rafael Ribeiro Kluge - GRR20244439
- Laisa Marcelino Santos Rodrigues - GRR20243783

## Estrutura do Projeto

O projeto implementa um solucionador de sistemas lineares utilizando o método dos Gradientes Conjugados (CG) e suas variantes pré-condicionadas (PCG), incluindo Jacobi, Gauss-Seidel (SGS) e SSOR. O código está modularizado nos seguintes arquivos:

### 1. helpers.h / helpers.c

Funções utilitárias para operações com vetores e matrizes.

- vet_produto: Produto interno de dois vetores.
- vet_norma2: Norma euclidiana de um vetor.
- vet_copy: Cópia de um vetor para outro.
- vet_preenche: Preenche um vetor com um valor constante.
- vet_axpy: Soma linear de vetores (y := y + alpha*x).
- vet_escala: Multiplica vetor por escalar.
- matvet_densa: Multiplica matriz densa A (n×n) por vetor x: y = A*x
- extrai_diag_e_invD: Extrai diagonal e seu inverso de uma matriz.
- aplica_jacobi: Aplica pré-condicionador de Jacobi.
- forward_sweep_DL: Varredura forward para SGS/SSOR.
- backward_sweep_DU: Varredura backward para SGS/SSOR.

### 2. pcgc.h / pcgc.c

Implementação do método dos Gradientes Conjugados e seus pré-condicionadores.

- pcg_precond_t: Define o conjunto dos tipos de pré-condicionador.
- pcg_contexto_t: Estrutura de contexto para buffers e parâmetros do pré-condicionador.
- pcg_setup: Inicializa buffers e parâmetros do pré-condicionador.
- pcg_apply: Aplica o pré-condicionador ao resíduo.
- pcg_free: Libera memória do contexto do pré-condicionador.
- cg_solve: Solver unificado para CG puro e PCG.

### 3. cgSolver.c

Programa principal para leitura de parâmetros, montagem do sistema, chamada do solver e saída dos resultados.

- residuo_l2: Calcula a norma do resíduo final.
- escolhe_precond: Define o valor de ω para o tipo de pré-condicionador.
- main: Fluxo principal do programa, incluindo leitura de entrada, preparação do sistema, chamada do solver e impressão dos resultados.

### 4. sislin.h / sislin.c

Módulo responsável pela geração e manipulação dos sistemas lineares k-diagonais e pela transformação para sistemas simétricos definidos positivos (SPD).

- criaKDiagonal: Gera uma matriz k-diagonal e o vetor de termos independentes.
- genSimetricaPositiva: Transforma o sistema original em um sistema simétrico definido positivo (AᵗA, Aᵗb).
- generateRandomA: Função interna que gera os coeficientes dos elementos da matriz k-diagonal.
- generateRandomB: Função interna que gera os termos independentes do sistema linear.

### 5. utils.h / utils.c

Funções utilitárias para medição de tempo e manipulação de strings.

- timestamp: Retorna o tempo atual em milissegundos, útil para medir tempos de execução de trechos do código.
- markerName: Gera uma string em um formato específico, útil para marcadores de ferramentas como LIKWID.
- Essas funções foram fornecidas no enunciado do trabalho.

## Observações

- Cada função está documentada com comentários no código-fonte, explicando seu propósito, parâmetros e funcionamento.
- O projeto segue as recomendações do enunciado e a preparação do sistema para garantir simetria e definição positiva.