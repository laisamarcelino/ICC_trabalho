#!/usr/bin/env bash
set -euo pipefail

# --- Config: arquivos de saída no MESMO diretório do script ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPORT_TXT="${SCRIPT_DIR}/cgsolver_test_report.txt"
REPORT_CSV="${SCRIPT_DIR}/cgsolver_test_report.csv"

CGSOLVER="${1:-./cgSolver}"
if [[ ! -x "$CGSOLVER" ]]; then
  echo "ERRO: executável não encontrado ou sem permissão: '$CGSOLVER'" >&2
  exit 2
fi

# --- Helpers ---
bold() { printf "\033[1m%s\033[0m" "$*"; }
ok()   { printf "  ✅ %s\n" "$*"; }
warn() { printf "  ⚠️  %s\n" "$*"; }
err()  { printf "  ❌ %s\n" "$*"; }

validate_output() {
  local out_file="$1" n_expected="$2"
  local nlines; nlines=$(wc -l < "$out_file" || echo 0)
  if [[ "$nlines" -ne 7 ]]; then
    err "Saída possui $nlines linhas (esperado 7)."
    return 1
  fi
  local L1; L1=$(sed -n '1p' "$out_file")
  if ! [[ "$L1" =~ ^[0-9]+$ ]]; then
    err "Linha 1 não é inteiro: '$L1'"; return 1;
  fi
  local L2; L2=$(sed -n '2p' "$out_file")
  local count; count=$(awk '{print NF}' <<< "$L2")
  if (( count != n_expected )); then
    warn "Vetor x tem $count elementos (esperado $n_expected)."
  fi
  for i in 3 4 5 6 7; do
    local Li; Li=$(sed -n "${i}p" "$out_file")
    if ! [[ "$Li" =~ ^[+\-]?[0-9]*\.?[0-9]+([eE][+\-]?[0-9]+)?$ ]]; then
      err "Linha $i não parece float: '$Li'"; return 1;
    fi
  done
  ok "Formato OK (7 linhas, números válidos)."
  return 0
}

run_case() {
  local NAME="$1" n="$2" k="$3" w="$4" maxit="$5" eps="$6"

  echo "" | tee -a "$REPORT_TXT" >/dev/null
  echo "===== CASE: $NAME =====" | tee -a "$REPORT_TXT" >/dev/null
  echo "params: n=$n k=$k ω=$w maxit=$maxit ε=$eps" | tee -a "$REPORT_TXT" >/dev/null

  local tmp_out tmp_err
  tmp_out="$(mktemp)"; tmp_err="$(mktemp)"
  {
    printf "%d\n%d\n%s\n%d\n%.16g\n" "$n" "$k" "$w" "$maxit" "$eps" \
      | "$CGSOLVER" 1> "$tmp_out" 2> "$tmp_err"
  } || true
  local rc=$?

  echo "--- exit_code: $rc"     | tee -a "$REPORT_TXT" >/dev/null
  echo "--- stdout:"            | tee -a "$REPORT_TXT" >/dev/null
  sed 's/^[[:space:]]*$/<empty>/' "$tmp_out" | tee -a "$REPORT_TXT" >/dev/null
  echo "--- stderr:"            | tee -a "$REPORT_TXT" >/dev/null
  sed 's/^[[:space:]]*$/<empty>/' "$tmp_err" | tee -a "$REPORT_TXT" >/dev/null

  if [[ $rc -eq 0 ]]; then
    if validate_output "$tmp_out" "$n"; then
      local norma resid tpc titer tres
      norma=$(sed -n '3p' "$tmp_out")
      resid=$(sed -n '4p' "$tmp_out")
      tpc=$(sed -n '5p' "$tmp_out")
      titer=$(sed -n '6p' "$tmp_out")
      tres=$(sed -n '7p' "$tmp_out")
      echo "OK;$NAME;$n;$k;$w;$maxit;$eps;$norma;$resid;$tpc;$titer;$tres" >> "$REPORT_CSV"
    else
      echo "BAD_FORMAT;$NAME;$n;$k;$w;$maxit;$eps;;;;;" >> "$REPORT_CSV"
    fi
  else
    echo "ERROR_RC_$rc;$NAME;$n;$k;$w;$maxit;$eps;;;;;" >> "$REPORT_CSV"
  fi

  rm -f "$tmp_out" "$tmp_err"
}

# --- Início do relatório ---
: > "$REPORT_TXT"
echo "status;name;n;k;omega;maxit;eps;norma;resid;tempo_pc;tempo_iter;tempo_residuo" > "$REPORT_CSV"
echo "Relatório de Testes cgSolver"            | tee -a "$REPORT_TXT" >/dev/null
echo "Executável: $CGSOLVER"                   | tee -a "$REPORT_TXT" >/dev/null
echo "Data: $(date '+%Y-%m-%d %H:%M:%S %Z')"   | tee -a "$REPORT_TXT" >/dev/null

# --- Casos válidos ---
run_case "no-precond"     200 5  -1   2000 1e-8
run_case "jacobi"         200 5   0.0 2000 1e-8
run_case "gauss-seidel"   200 5   1.0 2000 1e-8
run_case "ssor-1.5"       200 5   1.5 2000 1e-8
run_case "jacobi-n800-k9" 800 9   0.0 3000 1e-6

# --- Casos inválidos (espera RC != 0) ---
run_case "erro-n-pequeno"     10  5 0.0 1000 1e-8
run_case "erro-k-par"        200  6 0.0 2000 1e-8
run_case "erro-k-min"        200  1 0.0 2000 1e-8
run_case "erro-omega-bizarro" 200 5 2.5 2000 1e-8

echo "" | tee -a "$REPORT_TXT" >/dev/null
echo "Arquivos gerados:"      | tee -a "$REPORT_TXT" >/dev/null
echo " - $REPORT_TXT"         | tee -a "$REPORT_TXT" >/dev/null
echo " - $REPORT_CSV"         | tee -a "$REPORT_TXT" >/dev/null
bold "Concluído. Veja '$REPORT_TXT' e '$REPORT_CSV'."
