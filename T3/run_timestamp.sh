#!/bin/bash

if [[ $# -ne 1 ]]; then
  echo "Uso: $0 v1|v2"
  exit 1
fi

VER="$1"
if [[ "$VER" != "v1" && "$VER" != "v2" ]]; then
  echo "Versão inválida: $VER (use v1 ou v2)"
  exit 1
fi

# Lista dos tamanhos de matriz a testar
Ns=(32 64 128 256 512 1000 2000 4000 8000 9000 10000 20000)

# Parâmetros fixos
K=7
W=0
MAXIT=25
EPS=1e-6

# Caminho do executável e saída
DIR="./${VER}"
EXE="${DIR}/cgSolver"
RESULTS_DIR="${DIR}/results"
OUTFILE="${RESULTS_DIR}/${VER}_timestamps.csv"

# Compila antes de rodar os testes
cd "$DIR"
make clean
make
cd - > /dev/null

# Garante que o diretório de resultados existe
mkdir -p "$RESULTS_DIR"

# Cabeçalho do CSV
echo "N,tempo_op1_ms,tempo_op2_ms" > "$OUTFILE"

for N in "${Ns[@]}"; do
    # Executa o programa e captura as linhas de tempo
    OUTPUT=$(echo -e "$N\n$K\n$W\n$MAXIT\n$EPS" | "$EXE")
    # Extrai os tempos das linhas correspondentes
    TEMPO_OP1=$(echo "$OUTPUT" | grep "Tempo médio op1" | awk '{print $(NF-1)}')
    TEMPO_OP2=$(echo "$OUTPUT" | grep "Tempo op2" | awk '{print $(NF-1)}')
    # Salva no CSV
    echo "$N,$TEMPO_OP1,$TEMPO_OP2" >> "$OUTFILE"
    echo "Teste N=$N concluído."
done

echo "Todos os testes finalizados. Resultados em $OUTFILE"
