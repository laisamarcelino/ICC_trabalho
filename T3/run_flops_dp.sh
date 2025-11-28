#!/usr/bin/env bash
# Mede MFLOP/s usando LIKWID (FLOPS_DP) para v1 ou v2
# e salva o log completo de cada execução.
#
# Uso:
#   ./run_flops_dp.sh v1
#   ./run_flops_dp.sh v2

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

NS=(32 64 128 256 512 1000 2000 4000 8000 9000 10000 20000)
K=7
W=0
MAXIT=25
EPS=1e-6

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VER_DIR="${ROOT}/${VER}"
RESULTS_DIR="${VER_DIR}/results"
mkdir -p "${RESULTS_DIR}"

LOG_DP="${RESULTS_DIR}/${VER}_flops_dp.log"
: > "${LOG_DP}"

echo "[${VER}] Compilando com LIKWID..."
cd "${VER_DIR}"
make clean > /dev/null
USE_LIKWID=1 make > /dev/null
cd - > /dev/null

for N in "${NS[@]}"; do
  echo "[${VER}] N=${N} (FLOPS_DP)"
  {
    echo
    echo "=============================="
    echo "N=${N}"
    echo "=============================="
    likwid-perfctr -C 0 -g FLOPS_DP -m "${VER_DIR}/cgSolver" <<EOF
${N} ${K} ${W} ${MAXIT} ${EPS}
EOF
  } >> "${LOG_DP}" 2>&1 || {
    echo "Erro ao medir FLOPS_DP para N=${N}" >> "${LOG_DP}"
  }
done

echo
echo "[${VER}] Testes de FLOPS_DP finalizados."
echo "Log completo salvo em: ${LOG_DP}"