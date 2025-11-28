#!/usr/bin/env bash
# Mede MFLOP/s (FLOPS_AVX) para v1 ou v2 nos N do enunciado (usando LIKWID Marker API).

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

LOG_FILE="${RESULTS_DIR}/${VER}_FLOPS_AVX_flops.log"

echo "[${VER}] Compilando para grupo FLOPS_AVX..."
cd "${VER_DIR}"
make clean
USE_LIKWID=1 make
cd - > /dev/null

if [[ ! -x "${VER_DIR}/cgSolver" ]]; then
  echo "Erro: ${VER_DIR}/cgSolver não encontrado ou sem permissão de execução."
  exit 2
fi

echo "[${VER}] Rodando testes de MFLOP/s (grupo FLOPS_AVX) com Marker API..."
: > "${LOG_FILE}"

for N in "${NS[@]}"; do
  echo "### versao=${VER} grupo=FLOPS_AVX N=${N}" | tee -a "${LOG_FILE}"

  (
    cd "${VER_DIR}"
    echo "Executando: likwid-perfctr -C 0 -g FLOPS_AVX -m ./cgSolver"
    likwid-perfctr -C 0 -g FLOPS_AVX -m ./cgSolver <<EOF
${N} ${K} ${W} ${MAXIT} ${EPS}
EOF
  ) >> "${LOG_FILE}" 2>&1

  echo >> "${LOG_FILE}"
done

echo "[${VER}] Concluído. Resultados em: ${LOG_FILE}"