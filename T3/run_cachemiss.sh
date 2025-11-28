#!/usr/bin/env bash
# Mede "Data cache miss ratio" (via grupo L2CACHE) para v1 ou v2
# nos N do enunciado e extrai o valor do cache miss ratio do LIKWID.

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

# Grupo de desempenho a usar (L2CACHE, L1CACHE ou CACHE)
GROUP="L2CACHE"

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VER_DIR="${ROOT}/${VER}"
RESULTS_DIR="${VER_DIR}/results"
LOG_FILE="${RESULTS_DIR}/${VER}_${GROUP}_cachemiss.log"

mkdir -p "${RESULTS_DIR}"

echo "[${VER}] Diretório raiz: ${ROOT}"
echo "[${VER}] Diretório da versão: ${VER_DIR}"
echo "[${VER}] Grupo LIKWID: ${GROUP}"
echo "[${VER}] Compilando com LIKWID..."

cd "${VER_DIR}"
make clean
USE_LIKWID=1 make

echo "[${VER}] Rodando testes de 'cache miss ratio' (grupo ${GROUP})..."
: > "${LOG_FILE}"

for N in "${NS[@]}"; do
  echo "### versao=${VER} grupo=${GROUP} N=${N}" | tee -a "${LOG_FILE}"

  # Executa e salva saída completa
  likwid-perfctr -C 0 -g "${GROUP}" -m ./cgSolver <<EOF >> "${LOG_FILE}"
${N} ${K} ${W} ${MAXIT} ${EPS}
EOF

  # Removido: extração e print do ratio
  echo >> "${LOG_FILE}"
done

echo "[${VER}] Concluído. Resultados em: ${LOG_FILE}"
