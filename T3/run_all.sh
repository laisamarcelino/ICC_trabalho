#!/usr/bin/env bash
# Script mestre otimizado para rodar todos os testes (timestamp, memória, cache miss, FLOPS) para v1 ou v2
# Sempre usando LIKWID
# Uso: ./run_all.sh v1|v2

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Uso: $0 v1|v2"
  exit 1
fi

VER="$1"
if [[ "$VER" != "v1" && "$VER" != "v2" ]]; then
  echo "Versão inválida: $VER (use v1 ou v2)"
  exit 1
fi

# Lista de tamanhos de matriz
NS=(32 64 128 256 512 1000 2000 4000 8000 9000 10000 20000)
K=7
W=0
MAXIT=25
EPS=1e-6

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VER_DIR="${ROOT}/${VER}"
RESULTS_DIR="${VER_DIR}/results"
mkdir -p "${RESULTS_DIR}"

# ================================
# Compilação única com LIKWID
# ================================
echo "[${VER}] Compilando com LIKWID..."
cd "$VER_DIR"
USE_LIKWID=1 make clean
USE_LIKWID=1 make
cd - > /dev/null
EXEC="${VER_DIR}/cgSolver"

# ================================
# 1) Testes de timestamp
# ================================
echo "[${VER}] Rodando testes de tempo (timestamp)..."
OUT_TS="${RESULTS_DIR}/${VER}_timestamps.csv"
echo "N,tempo_op1_ms,tempo_op2_ms" > "$OUT_TS"

for N in "${NS[@]}"; do
    OUTPUT=$(echo -e "$N\n$K\n$W\n$MAXIT\n$EPS" | "$EXEC")
    TEMPO_OP1=$(echo "$OUTPUT" | grep "Tempo médio op1" | awk '{print $(NF-1)}')
    TEMPO_OP2=$(echo "$OUTPUT" | grep "Tempo op2" | awk '{print $(NF-1)}')
    echo "$N,$TEMPO_OP1,$TEMPO_OP2" >> "$OUT_TS"
    echo "Timestamp: N=$N concluído."
done

# ================================
# 2) Testes de memória (L3)
# ================================
echo "[${VER}] Rodando testes de memória (L3)..."
LOG_MEM="${RESULTS_DIR}/${VER}_L3.log"
: > "$LOG_MEM"

for N in "${NS[@]}"; do
    echo "### versao=${VER} grupo=L3 N=${N}" | tee -a "$LOG_MEM"
    likwid-perfctr -C 0 -g L3 -m "$EXEC" <<EOF >> "$LOG_MEM"
${N} ${K} ${W} ${MAXIT} ${EPS}
EOF
    echo >> "$LOG_MEM"
done

# ================================
# 3) Testes de cache miss (L2CACHE)
# ================================
echo "[${VER}] Rodando testes de cache miss (L2CACHE)..."
LOG_CM="${RESULTS_DIR}/${VER}_L2CACHE_cachemiss.log"
: > "$LOG_CM"

for N in "${NS[@]}"; do
    echo "### versao=${VER} grupo=L2CACHE N=${N}" | tee -a "$LOG_CM"
    likwid-perfctr -C 0 -g L2CACHE -m "$EXEC" <<EOF >> "$LOG_CM"
${N} ${K} ${W} ${MAXIT} ${EPS}
EOF
    echo >> "$LOG_CM"
done

# ================================
# 4) Testes de FLOPS (FLOPS_DP)
# ================================
echo "[${VER}] Rodando testes de FLOPS_DP..."
LOG_DP="${RESULTS_DIR}/${VER}_flops.log"
: > "$LOG_DP"

for N in "${NS[@]}"; do
    echo "[${VER}] N=${N} (FLOPS_DP)"
    {
        echo
        echo "=============================="
        echo "N=${N}"
        echo "=============================="
        likwid-perfctr -C 0 -g FLOPS_DP -m "$EXEC" <<EOF
${N} ${K} ${W} ${MAXIT} ${EPS}
EOF
    } >> "$LOG_DP" 2>&1 || {
        echo "Erro ao medir FLOPS_DP para N=${N}" >> "$LOG_DP"
    }
done

echo "[${VER}] Todos os testes finalizados."
echo "Resultados:"
echo " - Timestamp: $OUT_TS"
echo " - Memória (L3): $LOG_MEM"
echo " - Cache miss (L2CACHE): $LOG_CM"
echo " - FLOPS_DP: $LOG_DP"
