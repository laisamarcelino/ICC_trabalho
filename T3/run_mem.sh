#!/usr/bin/env bash
# Mede banda de memória (via grupo L3, substituindo MEM) para v1 ou v2
# nos N do enunciado.
#
# Uso:
#   ./run_mem.sh v1
#   ./run_mem.sh v2
#
# Parâmetros fixos:
#   N = 32 64 128 256 512 1000 2000 4000 8000 9000 10000 20000
#   k = 7
#   w = 0 (Jacobi)
#   maxit = 25
#   eps = 1e-6
#
# Resultados:
#   T3/v1/results/v1_L3.log
#   T3/v2/results/v2_L3.log
#
# OBS: Grupo MEM não existe nesta instalação; usamos L3 conforme enunciado.

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

# Lista de N do enunciado
NS=(32 64 128 256 512 1000 2000 4000 8000 9000 10000 20000)

K=7
W=0
MAXIT=25
EPS=1e-6

# Grupo de desempenho a usar (L3 no lugar de MEM)
GROUP="L3"

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VER_DIR="${ROOT}/${VER}"
RESULTS_DIR="${VER_DIR}/results"
LOG_FILE="${RESULTS_DIR}/${VER}_${GROUP}.log"

mkdir -p "${RESULTS_DIR}"

echo "[${VER}] Diretório raiz: ${ROOT}"
echo "[${VER}] Diretório da versão: ${VER_DIR}"
echo "[${VER}] Grupo LIKWID: ${GROUP}"
echo "[${VER}] Compilando com LIKWID..."

cd "${VER_DIR}"
make clean
USE_LIKWID=1 make

echo "[${VER}] Rodando testes de 'banda de memória' (grupo ${GROUP})..."
: > "${LOG_FILE}"

for N in "${NS[@]}"; do
  echo "### versao=${VER} grupo=${GROUP} N=${N}" | tee -a "${LOG_FILE}"

  likwid-perfctr -C 0 -g "${GROUP}" -m ./cgSolver <<EOF >> "${LOG_FILE}"
${N} ${K} ${W} ${MAXIT} ${EPS}
EOF

  echo >> "${LOG_FILE}"
done

echo "[${VER}] Concluído. Resultados em: ${LOG_FILE}"
