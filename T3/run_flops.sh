#!/usr/bin/env bash
# Mede MFLOP/s (FLOPS_DP e FLOPS_AVX) para v1 ou v2 nos N do enunciado.

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

# Tenta grupos em ordem de preferência
GROUPS=("FLOPS_DP" "FLOPS_AVX" "FLOPS" "DP" "AVX")
GROUP_FOUND=""

for G in "${GROUPS[@]}"; do
  if likwid-perfctr -g list | grep -q "$G"; then
    GROUP_FOUND="$G"
    break
  fi
done

if [[ -z "$GROUP_FOUND" ]]; then
  echo "Nenhum grupo FLOPS disponível no LIKWID para esta arquitetura."
  exit 2
fi

LOG_FILE="${RESULTS_DIR}/${VER}_${GROUP_FOUND}_flops.log"
echo "[${VER}] Usando grupo: ${GROUP_FOUND}"

cd "${VER_DIR}"
make clean
USE_LIKWID=1 make
cd - > /dev/null

echo "[${VER}] Rodando testes de MFLOP/s (grupo ${GROUP_FOUND})..."
: > "${LOG_FILE}"

for N in "${NS[@]}"; do
  echo "### versao=${VER} grupo=${GROUP_FOUND} N=${N}" | tee -a "${LOG_FILE}"

  likwid-perfctr -C 0 -g "${GROUP_FOUND}" -m ./cgSolver <<EOF >> "${LOG_FILE}"
${N} ${K} ${W} ${MAXIT} ${EPS}
EOF

  echo >> "${LOG_FILE}"
done

echo "[${VER}] Concluído. Resultados em: ${LOG_FILE}"
