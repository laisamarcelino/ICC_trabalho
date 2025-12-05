# Trabalho 02 - Otimização de Desempenho para Solução de Sistemas Lineares Esparsos com Pré-condicionadores

## Autoria

- Rafael Ribeiro Kluge - GRR20244439
- Laisa Marcelino Santos Rodrigues - GRR20243783

## Estrutura do Projeto

O projeto é a versão 2 que implementa um solucionador de sistemas lineares utilizando o método dos Gradientes Conjugados (CG) e suas variantes pré-condicionadas (PCG), incluindo Jacobi, Gauss-Seidel (SGS) e SSOR, de maneira otimizada em relação a versão 1. O código está modularizado nos seguintes arquivos:

### 1. helpers.h / helpers.c

Funções utilitárias para operações com vetores e matrizes.

- vet_produto: Produto interno de dois vetores.
- vet_norma2: Norma euclidiana de um vetor.
- vet_copy: Cópia de um vetor para outro.
- vet_preenche: Preenche um vetor com um valor constante.
- vet_axpy: Soma linear de vetores (y := y + alpha*x).
- vet_escala: Multiplica vetor por escalar.
- matvet_diagonais: Multiplica matriz k-diagonal (armazenada em matdiag_t) por vetor: y = A*x.
- extrai_diag_e_invD_diag: Extrai a diagonal principal de uma matdiag_t e calcula seu inverso.
- aplica_jacobi: Aplica pré-condicionador de Jacobi (y = D^{-1} r).
- varredura_progressiva_DL / varredura_regressiva_DU: Varreduras forward/backward usadas por SGS/SSOR.
- residuo_l2_v2: Calcula a norma L2 do resíduo r = b - A x usando matvet_diagonais.
- liberaMatDiag: Libera memória de uma estrutura matdiag_t.

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

- escolhe_precond: Define o valor de ω para o tipo de pré-condicionador.
- main: Fluxo principal do programa, incluindo leitura de entrada, preparação do sistema, chamada do solver e impressão dos resultados.

### 4. sislin.h / sislin.c

Módulo responsável pela geração e manipulação dos sistemas lineares k-diagonais e pela transformação para sistemas simétricos definidos positivos (SPD).

- criaKDiagonal_v2: Gera uma matriz k-diagonal armazenada em matdiag_t e o vetor de termos independentes.
- genSimetricaPositiva_diag: Constrói AᵗA e Aᵗb diretamente na estrutura matdiag_t.
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

## Referências Teóricas

- CUNHA, M.CRISTINA; Métodos Numéricos, 2002, Ed. UNICAMP, Cap. 3, seções 3.3, 3.4 e 3.5.

- CARVALHO, Rafael Aleixo de; VIEIRA, Felipe. Métodos Iterativos para Problemas de Quadrados Mínimos Lineares. São Carlos: Sociedade Brasileira de Matemática Aplicada e Computacional (SBMAC), 2022. 200 p. (Notas em Matemática Aplicada, v. 93). ISBN 978-65-86388-09-1.