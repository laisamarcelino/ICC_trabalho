Lista de tarefas

Configuração do ambiente:
    Criação do repositório com versão 1 (v1) e versão 2 (v2). ✅
    Compilador GCC e flags: -O3 -march=native -mavx -fopt-info-vec. ✅
    Biblioteca LIKWID instalada e testada (likwid-topology, likwid-perfctr).

Otimizações:
    (OP1) Na operação de iteração do método de Gradiente Conjugado com Pré-condicionador de Jacobi (op1);
        - Substituir matrizes completas por formato de armazenamento de diagonais (band storage). ✅
        - Armazenar apenas 7 diagonais (k=7), cada uma em um vetor. ✅
        - Evitar armazenar diagonais nulas. ✅
        - Aplicar vetorização manual ou automática nas operações internas ✅
            - Otimizações feitas na prova -> aplicar em helpers.c ~✅ (UNROLL APENAS)

    (OP2) Na operação de cálculo do resíduo (op2): R = b - Ax; 
        - Reescrever multiplicação Ax considerando a estrutura esparsa (somente diagonais válidas).✅
        - Reduzir acessos à memória (usar dados contíguos e alinhamento AVX).✅
        - Garantir que a operação b - Ax use vetorização. ✅
            - Otimizações feitas na prova -> aplicar em helpers.c ✅

Testes de desempenho:
    Inserir marcadores LIKWID_MARKER_START/STOP em:
    op1 — iteração do CG;
    op2 — cálculo do resíduo.

    Rodar testes com os tamanhos:
    N = {32, 64, 128, 256, 512, 1000, 2000, 4000, 8000, 9000, 10000, 20000}.
    Cada teste: 25 iterações, k=7.

    Medir e registrar:
        Tempo médio (timestamp() em ms);
        Bandwidth (MEM/L3);
        Cache miss L2/L1;
        MFLOP/s (FLOPS_DP e FLOPS_AVX).

    Gerar gráficos (log-log):
        Tempo × N (para op1 e op2);
        Bandwidth × N;
        Cache miss × N;
        MFLOP/s × N.

