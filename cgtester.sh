#!/usr/bin/env bash
# cgtester.sh — Testes automatizados para o Trabalho 1 (CG + Pré-condicionadores)
# Ambiente-alvo: Arch Linux
#
# Uso:
#   ./cgtester.sh caminho/para/login1-login2.tgz
#   ./cgtester.sh caminho/para/diretorio/login1-login2
#
# O script:
#  - Garante dependências básicas
#  - Extrai o .tgz (se fornecido) para um tmpdir
#  - Executa `make clean && make all` e confere cgSolver
#  - Roda testes obrigatórios (-1 e 0.0) e opcionais (1.0 e 1.5)
#  - Injeta entradas via stdin (n, k, ω, maxit, ε)
#  - Valida formato de saída (linhas e contagens)
#  - Checa tratamento de erros (n<=10, k par, k<=1, ω fora do escopo)
#  - Resume PASS/FAIL/SKIP; sai com código !=0 se algum obrigatório falhar
#
# Requisitos:
#  - base-devel (make, gcc), tar, gzip, coreutils, awk, grep, sed, timeout
#  - Em Arch: sudo pacman -S --needed base-devel tar gzip gawk coreutils grep sed

set -Eeuo pipefail

# ---------- Utilidades ----------
require_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "Falta dependência: $1" >&2; exit 1; }; }
for c in tar gzip make awk grep sed timeout head tail wc tr; do require_cmd "$c"; done

readonly GREEN="$(tput setaf 2 || true)"
readonly YELLOW="$(tput setaf 3 || true)"
readonly RED="$(tput setaf 1 || true)"
readonly CYAN="$(tput setaf 6 || true)"
readonly BOLD="$(tput bold || true)"
readonly RESET="$(tput sgr0 || true)"

log()  { echo -e "${CYAN}${*}${RESET}"; }
ok()   { echo -e "${GREEN}✔${RESET} $*"; }
warn() { echo -e "${YELLOW}⚠${RESET} $*"; }
err()  { echo -e "${RED}✘${RESET} $*"; }

usage() {
  cat <<EOF
${BOLD}Uso:${RESET} $0 <caminho para .tgz OU diretório do projeto>
Exemplos:
  $0 ./login1-login2.tgz
  $0 ~/trabalhos/login1-login2
EOF
}

[[ $# -eq 1 ]] || { usage; exit 2; }

INPUT_PATH="$(realpath "$1")"
WORKDIR="$(mktemp -d -t cgtester-XXXXXX)"
ARTIFACTS="$WORKDIR/artifacts"
mkdir -p "$ARTIFACTS"
cleanup() { [[ -d "$WORKDIR" ]] && rm -rf "$WORKDIR"; }
trap cleanup EXIT

export LC_ALL=C
export LC_NUMERIC=C

# ---------- Encontrar/Preparar projeto ----------
PROJECT_DIR=""
if [[ -d "$INPUT_PATH" ]]; then
  PROJECT_DIR="$INPUT_PATH"
elif [[ -f "$INPUT_PATH" ]]; then
  case "$INPUT_PATH" in
    *.tgz|*.tar.gz)
      log "Extraindo pacote em $WORKDIR ..."
      tar -xzf "$INPUT_PATH" -C "$WORKDIR"
      # Descobrir diretório raiz criado pela extração
      CANDIDATES=("$WORKDIR"/*)
      [[ -d "${CANDIDATES[0]}" ]] || { err "Não foi possível localizar diretório após extração."; exit 1; }
      PROJECT_DIR="${CANDIDATES[0]}"
      ;;
    *)
      err "Arquivo não suportado: $INPUT_PATH (use .tgz/.tar.gz ou um diretório)."
      exit 1
      ;;
  esac
else
  err "Caminho inválido: $INPUT_PATH"
  exit 1
fi

log "Diretório do projeto: $PROJECT_DIR"
[[ -f "$PROJECT_DIR/Makefile" ]] || { err "Makefile não encontrado em $PROJECT_DIR"; exit 1; }

# Heurística: avisar se não houver LEIAME
[[ -f "$PROJECT_DIR/LEIAME" ]] || warn "Arquivo LEIAME não encontrado (apenas aviso)."

# Opcional: checar presença de srandom(20252) (apenas aviso)
if ! grep -R --include='*.c' -n 'srandom( *20252 *);' "$PROJECT_DIR" >/dev/null 2>&1; then
  warn "Não encontrei 'srandom(20252);' nos fontes (.c). Verifique a especificação."
fi

# ---------- Compilação ----------
pushd "$PROJECT_DIR" >/dev/null
log "Compilando (make clean && make all) ..."
if ! make clean >/dev/null 2>&1; then warn "make clean retornou código !=0 (seguindo assim mesmo)"; fi
if ! make -j"$(nproc)" all; then err "Falha na compilação (make all)."; exit 1; fi
[[ -x "./cgSolver" ]] || { err "Executável cgSolver não encontrado após make all."; exit 1; }
ok "Compilação OK e cgSolver presente."

# ---------- Runner/Validações ----------
PASS=0; FAIL=0; SKIP=0

is_number() {
  # aceita inteiro/real/científico (ex.: -1, 0.5, 1e-3, -2.3E+10)
  printf '%s' "$1" | grep -Eq '^[+-]?([0-9]+(\.[0-9]*)?|\.[0-9]+)([eE][+-]?[0-9]+)?$'
}

validate_stdout() {
  # Esperado:
  # 1) linha 1: n
  # 2) linha 2: x_1 ... x_n (n números)
  # 3) linha 3: norma (1 número)
  # 4) linha 4: resíduo (1 número)
  # 5) linha 5: tempo_pc (1 número)
  # 6) linha 6: tempo_iter (1 número)
  # 7) linha 7: tempo_residuo (1 número)
  local out="$1"
  local expected_n="$2"
  local eps="$3"

  local lines cnt
  cnt="$(wc -l < "$out" | tr -d ' ')"
  if [[ "$cnt" -lt 7 ]]; then
    err "Saída possui menos de 7 linhas (=$cnt)."
    return 1
  fi

  local L1 L2 L3 L4 L5 L6 L7
  L1="$(sed -n '1p' "$out" | tr -d '\r')"
  L2="$(sed -n '2p' "$out" | tr -d '\r')"
  L3="$(sed -n '3p' "$out" | tr -d '\r')"
  L4="$(sed -n '4p' "$out" | tr -d '\r')"
  L5="$(sed -n '5p' "$out" | tr -d '\r')"
  L6="$(sed -n '6p' "$out" | tr -d '\r')"
  L7="$(sed -n '7p' "$out" | tr -d '\r')"

  # L1: n
  if ! [[ "$L1" =~ ^[0-9]+$ ]]; then
    err "Linha 1 não é um inteiro (n): '$L1'"
    return 1
  fi
  if [[ "$L1" -ne "$expected_n" ]]; then
    warn "Linha 1 (n=$L1) diverge do esperado ($expected_n)."
  fi

  # L2: vetor x com n números
  local num_tokens
  num_tokens="$(awk '{print NF}' <<<"$L2")"
  if [[ "$num_tokens" -ne "$expected_n" ]]; then
    err "Linha 2 não possui $expected_n números (tem $num_tokens)."
    return 1
  fi
  # checar se todos são números válidos
  for tok in $L2; do
    is_number "$tok" || { err "Elemento de x não numérico: '$tok'"; return 1; }
  done

  # L3..L7: números
  for label in "norma:$L3" "residuo:$L4" "tempo_pc:$L5" "tempo_iter:$L6" "tempo_residuo:$L7"; do
    name="${label%%:*}"; val="${label#*:}"
    is_number "$val" || { err "Linha '$name' não numérica: '$val'"; return 1; }
  done

  # Heurística (não reprova): norma não muito maior que eps
  awk -v nrm="$L3" -v e="$eps" 'BEGIN{ if (nrm > 10*e) exit 1; else exit 0; }' \
    || warn "norma (= $L3) > 10*ε (= $(awk -v e="$eps" 'BEGIN{printf "%.6g", 10*e}')) — pode indicar não convergência."

  # tempos não-negativos
  awk -v a="$L5" -v b="$L6" -v c="$L7" 'BEGIN{ if (a<0||b<0||c<0) exit 1; }' \
    && true || { err "Algum tempo é negativo."; return 1; }

  return 0
}

run_case() {
  local name="$1" n="$2" k="$3" w="$4" maxit="$5" eps="$6" required="$7"
  local out="$ARTIFACTS/${name}.out"
  local errf="$ARTIFACTS/${name}.err"
  local rc=0

  log "Rodando ${BOLD}${name}${RESET}  (n=$n, k=$k, ω=$w, maxit=$maxit, ε=$eps)"
  if ! printf "%s\n" "$n" "$k" "$w" "$maxit" "$eps" \
      | timeout 20s ./cgSolver 1>"$out" 2>"$errf"; then
    rc=$?
  fi

  # Tratamento de opcionais (Gauss-Seidel/SSOR) não implementados
  if [[ "$required" = "no" && "$rc" -ne 0 ]]; then
    if grep -Eiq "não implementad|nao implementad|opcional|unsupported|not implemented" "$errf"; then
      warn "Pré-condicionador opcional não implementado — marcando como SKIP."
      ((SKIP++)) || true
      return 0
    fi
  fi

  if [[ "$rc" -ne 0 ]]; then
    err "${name}: execução retornou código $rc"
    echo "---- stderr ----"; sed 's/^/  /' "$errf" | head -n 20; echo "--------------"
    ((FAIL++)) || true
    return 1
  fi

  # Se stderr tiver conteúdo para caso "válido", tratar como aviso (pode ser logging)
  if [[ -s "$errf" ]]; then
    warn "${name}: houve saída em stderr (ver ${errf})."
  fi

  if validate_stdout "$out" "$n" "$eps"; then
    ok "${name}: PASS"
    ((PASS++)) || true
    return 0
  else
    err "${name}: FAIL (formato/validação)"
    echo "---- stdout (7 primeiras linhas) ----"; head -n 7 "$out" | sed 's/^/  /'
    echo "---- stderr (até 20 linhas) ----"; head -n 20 "$errf" | sed 's/^/  /'
    ((FAIL++)) || true
    return 1
  fi
}

run_error_case() {
  local name="$1" n="$2" k="$3" w="$4" maxit="$5" eps="$6"
  local errf="$ARTIFACTS/${name}.err"
  local out="$ARTIFACTS/${name}.out"
  local rc=0
  log "Rodando (erro esperado) ${BOLD}${name}${RESET}  (n=$n, k=$k, ω=$w)"
  if printf "%s\n" "$n" "$k" "$w" "$maxit" "$eps" \
      | timeout 10s ./cgSolver 1>"$out" 2>"$errf"; then
    rc=0
  else
    rc=$?
  fi

  if [[ "$rc" -eq 0 ]]; then
    err "${name}: esperado erro (exit!=0), mas retornou 0"
    echo "---- stdout ----"; head -n 10 "$out" | sed 's/^/  /'
    ((FAIL++)) || true
    return 1
  fi

  if [[ ! -s "$errf" ]]; then
    err "${name}: esperado mensagem em stderr, mas stderr está vazio."
    ((FAIL++)) || true
    return 1
  fi

  ok "${name}: PASS (erro tratado corretamente)"
  ((PASS++)) || true
  return 0
}

# ---------- Testes ----------
# Casos "válidos" (obrigatórios)
run_case "T1_sem_precond"  "50" "5" "-1"   "5000" "1e-8" "yes"
run_case "T2_jacobi"       "60" "7" "0.0"  "5000" "1e-8" "yes"

# Casos "válidos" (opcionais)
run_case "T3_gauss_seidel" "70" "5" "1.0"  "5000" "1e-8" "no"   || true
run_case "T4_ssor"         "80" "9" "1.5"  "5000" "1e-8" "no"   || true

# Casos de ERRO (devem falhar graciosamente)
run_error_case "E1_n_menor_igual_10"  "10" "5" "-1"  "100" "1e-6" || true
run_error_case "E2_k_par"             "50" "4" "-1"  "100" "1e-6" || true
run_error_case "E3_k_menor_igual_1"   "50" "1" "-1"  "100" "1e-6" || true
run_error_case "E4_omega_fora_faixa"  "50" "5" "2.1" "100" "1e-6" || true

# ---------- Resumo ----------
echo
echo "${BOLD}Resumo:${RESET} ${GREEN}PASS=$PASS${RESET}  ${YELLOW}SKIP=$SKIP${RESET}  ${RED}FAIL=$FAIL${RESET}"
echo "  Artefatos/logs: $ARTIFACTS"
echo

# Saída com falha se algum obrigatório falhar
if [[ "$FAIL" -gt 0 ]]; then
  exit 1
else
  exit 0
fi
