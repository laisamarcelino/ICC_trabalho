#!/usr/bin/env bash
# Mede MFLOP/s usando LIKWID (FLOPS_AVX) para v1 ou v2
# e salva o log completo de cada execução.
#
# Uso:
#   ./run_flops_avx.sh v1
#   ./run_flops_avx.sh v2

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

LOG_AVX="${RESULTS_DIR}/${VER}_flops_avx.log"
: > "${LOG_AVX}"

echo "[${VER}] Compilando com LIKWID..."
cd "${VER_DIR}"
make clean > /dev/null
USE_LIKWID=1 make > /dev/null
cd - > /dev/null

for N in "${NS[@]}"; do
  echo "[${VER}] N=${N} (FLOPS_AVX)"
  {
    echo
    echo "=============================="
    echo "N=${N}"
    echo "=============================="
    likwid-perfctr -C 0 -g FLOPS_AVX -m "${VER_DIR}/cgSolver" <<EOF
${N} ${K} ${W} ${MAXIT} ${EPS}
EOF
  } >> "${LOG_AVX}" 2>&1 || {
    echo "Erro ao medir FLOPS_AVX para N=${N}" >> "${LOG_AVX}"
  }
done

echo
echo "[${VER}] Testes de FLOPS_AVX finalizados."
echo "Log completo salvo em: ${LOG_AVX}"