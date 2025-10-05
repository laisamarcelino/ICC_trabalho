# Plano em etapas — Trabalho 1

## Objetivo: deixar o projeto compilável e reproduzível localmente.

### 1. Entendimento rápido e preparação do repositório

Ações:

Colocar todos os arquivos (os que você mostrou) no diretório do projeto.

Verificar e corrigir inconsistências de tipos/headers (notei sislin.h inclui sislin.h novamente — isso é erro).

Garantir o Makefile gerará cgSolver no diretório raiz do pacote conforme a especificação (regra all produz cgSolver).

Confirmar presença dos arquivos listados em DISTFILES: *.c *.h Makefile LEIAME.

Arquivos a revisar/editar:

Makefile — ajustar MODULES, SRCS se necessário; garantir all: produz cgSolver.

sislin.h — remover inclusão recursiva e corrigir protótipos para usar real_t e rtime_t do utils.h (consistência com sislin.c).

LEIAME — criar/esboçar (autoria, compilação, uso).

Validação:

Rodar make e obter executável cgSolver (mesmo que incompleto).

make clean funciona.

### 2. Corrigir/normalizar headers e tipos

Objetivo: evitar erros de compilação por headers conflitantes e tipos inconsistentes.

Tarefas:

Em sislin.h:

Remover #include "sislin.h" recursivo.

Incluir utils.h.

Usar real_t e rtime_t nas assinaturas (atualmente há mistura double e real_t).

Assegurar protótipos batem com sislin.c (parâmetros e tipos).

Em sislin.c:

Uniformizar tipos (real_t, rtime_t) e protótipos conforme sislin.h.

Incluir utils.h (já tem).

Validação:

gcc -c de cada .c não deve falhar por tipos/headers.

### 3. Implementar geração da matriz k-diagonal e vetor b

Objetivo: implementar criaKDiagonal usando generateRandomA e generateRandomB, armazenando em formato compacto (vetor denso ou representação por bandas).

Decisões de representação (escolha que impacta eficiência):

Opção A (simples): armazenar A como vetor denso n*n — mais simples para implementar e depurar.

Opção B (recomendado para esparsidade): armazenar por bandas/k-diagonal em n*k (cada linha armazena k valores) — atende requisitos e facilita multiplicação eficiente.

O que implementar:

criaKDiagonal(int n, int k, real_t **A, real_t **B):

Alocar A no formato escolhido (p.ex. real_t *A = calloc(n*k, sizeof(real_t))) e preencher somente as k diagonais com generateRandomA(i,j,k).

Para B, alocar vetor real_t *B = calloc(n,...) e preencher com generateRandomB(k).

Preencher A de forma que A represente matriz simétrica (ou gerar e depois transformar) — lembrar que generateRandomA não garante SPD.

Documentar formato no comentário (importante para LEIAME).

Validação:

Implementar função de debug que imprime as primeiras linhas/diagonais.

Rodar teste: gerar n=20 k=3 e imprimir A e B.

### 4. Gerar matriz simétrica e positiva definida (SPD)

Objetivo: transformar A gerada em uma matriz SPD apropriada para o CG.

Tarefas:

Implementar genSimetricaPositiva(real_t *A, real_t *b, int n, int k, real_t **ASP, real_t *bsp, rtime_t *tempo):

Estratégia simples e eficaz: construir ASP = A^T * A + α I (produto que garante SPD). Porém produto denso tem custo; como é para teste e n pode ser moderado, ok.

Alternativa para k-diagonal: garantir diagonal dominante — por exemplo, para cada linha, somar valores absolutos das off-diagonais e acrescentar fator na diagonal para garantir estrita diagonal dominante (makes SPD likely). Especifique no código.

Implementar *tempo = timestamp(); conforme assinatura (medir apenas a transformação).

Produzir bsp (vetor b transformado consistentemente — se usar ASP = A^T*A, então bsp = A^T*b).

Validação:

Verificar ASP simétrica (comparar ASP[i,j] e ASP[j,i]).

Teste simples: checar que x^T * ASP * x > 0 para alguns vetores aleatórios x (teste de positividade).

Observação: o enunciado diz “Observe que o sistema linear resultante não atende — precisa ser transformado antes da aplicação” — a transformação via A^T A ou forcando diagonal dominante é aceitável; documente a escolha.

### 5. Extrair D, L, U (função geraDLU)

Objetivo: extrair D (diagonal), L (strict lower with zeros on diagonal) e U (strict upper with zeros) a partir de ASP armazenada.

Tarefas:

Implementar geraDLU (real_t *A, int n, int k, real_t **D, real_t **L, real_t **U, rtime_t *tempo):

Fazer iteração sobre a estrutura de armazenamento de ASP e preencher:

D: vetor length n com diagonais.

L: mesma estrutura de bandas (ou compacto) contendo apenas elementos i>j.

U: idem para i<j.

Medir tempo em *tempo como pedido.

Validação:

Recombinar L + D + U e comparar com ASP (diferença pequena).

### 6. Gerar pré-condicionador M (função geraPreCond)

Objetivo: construir/inverter (ou preparar aplicação de) M para os casos necessários (ω=-1 identidade; ω=0 Jacobi; ω=1 Gauss-Seidel; 1<ω<2 SSOR). Lembrar que cálculo de M⁻¹ completo pode ser custoso — o enunciado pede M⁻¹ (ou função que resolve M y = r).

Tarefas:

Implementar geraPreCond(real_t *D, real_t *L, real_t *U, real_t w, int n, int k, real_t **M, rtime_t *tempo):

Caso ω = -1: M = I (representação n*vetor com 1s) — trivial.

Caso ω = 0.0 (Jacobi): M = D — armazenar diagonal inversa se você preferir aplicar M⁻¹ facilmente (ideal armazenar M_inv).

Caso ω = 1.0 (Gauss-Seidel) e SSOR (ω>1.0): como bônus, implementar aproximação:

Opção simples: armazenar factors to solve M y = r via forward/back substitution if M is triangular or using ILU0 / incomplete Cholesky (opcional).

Se não fizer ILU/IC, documente que geraPreCond prepara estrutura para solução, mas a aplicação de M^{-1} será feita por resolver sistemas triangular aproximados. Para avaliação básica, você pode marcar Gauss-Seidel / SSOR como “bônus” e implementar apenas Jacobi + sem pré-condicionador.

Medir tempo de pré-condicionamento em *tempo (não incluir geração da matriz de entrada).

Design recomendação:

Em vez de calcular M⁻¹ explícito (que pode densificar), armazene uma rotina apply_precond(M_struct, r, y) que resolve M y = r eficientemente (Jacobi: y = D^{-1} r; Identity: y=r).

Validação:

Teste apply_precond com vetores simples e comparar com solução direta quando possível.

### 7. Implementar Gradientes Conjugados pré-condicionado

Objetivo: implementar o algoritmo de PCG (preconditioned conjugate gradients) com critérios de parada e tempos conforme enunciado.

Tarefas:

Criar módulo pcgc.c (parece já citado no Makefile) ou completar se existir:

Função pcg_solve(ASP, bsp, n, k, w, maxit, eps, x, tempo_iter_avg, tempo_residuo) ou similar.

Regras importantes:

Estimativa inicial x0 = 0.

Critério de parada: max(|xi - xi-1|) < ε e/ou ||r||_2 convergência; mas enunciado pede norma máxima em x para ε.

Calcular resíduo final: r = b - A x e ||r||_2.

Medir tempos:

tempo_pc: já medido na etapa de pré-condicionamento (geraPreCond).

tempo_iter: tempo médio por iteração (medir tempo total de iterações / número de iterações).

tempo_residuo: tempo para calcular norma euclidiana do resíduo ao final (usar timestamp()).

Tratamento de falhas:

Se o método divergir, imprimir mensagem em stderr e encerrar com código != 0.

Detectar NaN/Inf em vetores e abortar com mensagem.

Implementação:

Implementar versão sem pré-condicionador (M = I).

Implementar versão com Jacobi (M = D) — aplicar M_inv multiplicando.

Opcional: implementar Gauss-Seidel / SSOR (bônus), ou preparar estrutura para isso.

Validação:

Testes com n pequeno (ex.: n=50, k=3) e verificar convergência.

Comparar solução com lib (se disponível) ou com eliminação direta para pequenas dimensões.

### 8. Implementar cálculo do resíduo e métricas

Objetivo: função calcResiduoSL e cálculo da norma máxima entre iterações.

Tarefas:

Implementar real_t calcResiduoSL (real_t *A, real_t *b, real_t *X, int n, int k, rtime_t *tempo) para calcular r = b - A x e retornar ||r||_2.

Implementar função auxiliar para norma_max_diff(x_old, x_new, n) (infinito-norma do x - x_old).

Medir tempo_residuo usando timestamp().

Validação:

Testar com vetor x conhecido e comparar resultado contra cálculo manual ou multiplicação densa.

### 9. Entrada e saída conforme especificação

Objetivo: ler stdin os 5 valores (n k w maxit ε) e imprimir a saída no formato exato pedido.

Tarefas:

Implementar parsing robusto de stdin (ler com scanf, validar ranges: n>10, k>1 && k odd, ω possibilities).

Validar maxit e ε.

Saída (stdout) exatamente:

--------------------- versao 2 ----------------------------

## 📅 Cronograma de Desenvolvimento

### ✅ 05/10 (Sábado) — Planejamento
- Leitura completa do enunciado.
- Análise dos arquivos já fornecidos.
- Definição do escopo mínimo (ω = -1, 0.0) e escopo bônus (ω = 1.0, >1.0).
- Estudo do método dos Gradientes Conjugados e pré-condicionadores.

### ✅ 06–07/10 (Dom–Seg) — Geração de matriz k-diagonal
- Implementar `criaKDiagonal()` (`sislin.c`)  
  _Responsável: Rafael_  
- Implementar `genSimetricaPositiva()`  
  _Responsável: Laisa_  
- Testes de geração com `srandom(20252)`

### ✅ 08–09/10 (Ter–Qua) — Gradientes Conjugados (sem pré-condicionador)
- Implementar método CG com ω = -1  
  _Responsável: Laisa_  
- Cálculo de erro e critério de parada com ε  
  _Responsável: Rafael_

### ✅ 10–11/10 (Qui–Sex) — Pré-condicionador Jacobi
- Gerar matriz M = D  
  _Responsável: Rafael_  
- Ajustar CG para uso de M⁻¹b e M⁻¹r  
  _Responsável: Laisa_

### ⏳ 12–13/10 (Sáb–Dom) — Gauss-Seidel (opcional)
- Implementar `geraDLU()`  
  _Responsável: Laisa_  
- Implementar `geraPreCond()` com ω = 1.0  
  _Responsável: Rafael_

### ⏳ 14/10 (Terça) — SSOR (opcional)
- Generalizar `geraPreCond()` para ω > 1.0  
  _Responsável: Rafael_  
- Integrar SSOR ao solver  
  _Responsável: Laisa_

### ⏳ 15/10 (Quarta) — Medição de tempos
- Medir:
  - `tempo_pc`
  - `tempo_iter`
  - `tempo_residuo`  
  _Responsável: Laisa_
- Testar desempenho para diferentes parâmetros  
  _Responsável: Rafael_

### ⏳ 16/10 (Quinta) — Tratamento de erros
- Mensagens em `stderr` e encerramento com `exit(1)`  
  _Responsável: Laisa_  
- Testes de não convergência e falhas numéricas  
  _Responsável: Rafael_

### ⏳ 17/10 (Sexta) — Finalização
- Escrever arquivo `LEIAME` com:
  - Autores
  - RA
  - Descrição das funções e estruturas  
  _Responsável: Rafael_
- Verificar `Makefile` (`all`, `clean`, `purge`, `dist`)  
  _Responsável: Laisa_

### ⏳ 18/10 (Sábado) — Revisão e Entrega
- Testes finais (entrada, saída, erro)  
  _Responsável: Laisa e Rafael_  
- Compactação com `tar/gzip` no formato `login1-login2.tgz`
- Submissão no Moodle até 23:59

---


