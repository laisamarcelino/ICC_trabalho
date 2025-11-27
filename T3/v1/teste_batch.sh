#!/bin/bash

# Lista dos tamanhos de matriz a testar
Ns=(32 64 128 256 512 1000 2000 4000 8000 9000 10000 20000)

# Parâmetros fixos
K=7
W=-1
MAXIT=25
EPS=1e-6

# Arquivo de saída
OUTFILE="resultados_v2.csv"

# Cabeçalho do CSV
echo "N,tempo_op1_ms,tempo_op2_ms" > "$OUTFILE"

for N in "${Ns[@]}"; do
    # Executa o programa e captura as linhas de tempo
    OUTPUT=$(echo -e "$N\n$K\n$W\n$MAXIT\n$EPS" | ./cgSolver)
    # Extrai os tempos das linhas correspondentes
    TEMPO_OP1=$(echo "$OUTPUT" | grep "Tempo médio op1" | awk '{print $(NF-1)}')
    TEMPO_OP2=$(echo "$OUTPUT" | grep "Tempo op2" | awk '{print $(NF-1)}')
    # Salva no CSV
    echo "$N,$TEMPO_OP1,$TEMPO_OP2" >> "$OUTFILE"
    echo "Teste N=$N concluído."
done

echo "Todos os testes finalizados. Resultados em $OUTFILE"
