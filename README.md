## 📅 Cronograma de Desenvolvimento

### ✅ 05/10 — Planejamento L[x] R[]
- Leitura completa do enunciado. 
- Análise dos arquivos já fornecidos.
- Definição do escopo mínimo (ω = -1, 0.0) e escopo bônus (ω = 1.0, >1.0).
- Estudo do método dos Gradientes Conjugados e pré-condicionadores.

### ✅ 06–07/10 — Geração de matriz k-diagonal L[x] R[]
- Implementar `criaKDiagonal()` (`sislin.c`) 
  _Responsável: Rafael_  
- Implementar `genSimetricaPositiva()`  
  _Responsável: Laisa_  
- Testes de geração com `srandom(20252)`

### ✅ 08–09/10 — Gradientes Conjugados (sem pré-condicionador) L[x] R[]
- Implementar método CG com ω = -1  
  _Responsável: Laisa_  
- Cálculo de erro e critério de parada com ε  
  _Responsável: Rafael_

### ✅ 10–11/10 — Pré-condicionador Jacobi L[x] R[]
- Gerar matriz M = D  
  _Responsável: Rafael_  
- Ajustar CG para uso de M⁻¹b e M⁻¹r  
  _Responsável: Laisa_

### ⏳ 12–13/10 — Gauss-Seidel (opcional) L[x] R[]
- Implementar `geraDLU()`  
  _Responsável: Laisa_  
- Implementar `geraPreCond()` com ω = 1.0  
  _Responsável: Rafael_

### ⏳ 14/10  — SSOR (opcional) L[x] R[]
- Generalizar `geraPreCond()` para ω > 1.0  
  _Responsável: Rafael_  
- Integrar SSOR ao solver  
  _Responsável: Laisa_

### ⏳ 15/10 — Medição de tempos L[x] R[]
- Medir:
  - `tempo_pc`
  - `tempo_iter`
  - `tempo_residuo`  
  _Responsável: Laisa_
- Testar desempenho para diferentes parâmetros  
  _Responsável: Rafael_

### ⏳ 16/10 — Tratamento de erros L[x] R[]
- Mensagens em `stderr` e encerramento com `exit(1)`  
  _Responsável: Laisa_  
- Testes de não convergência e falhas numéricas  
  _Responsável: Rafael_

### ⏳ 17/10 — Finalização L[x] R[]
- Escrever arquivo `LEIAME` com:
  - Autores
  - RA
  - Descrição das funções e estruturas  
  _Responsável: Rafael_
- Verificar `Makefile` (`all`, `clean`, `purge`, `dist`)  
  _Responsável: Laisa_

### ⏳ 18/10 — Revisão e Entrega L[x] R[]
- Testes finais (entrada, saída, erro)  
  _Responsável: Laisa e Rafael_  
- Compactação com `tar/gzip` no formato `login1-login2.tgz`
- Submissão no Moodle até 23:59

---

## 📂 Estrutura Esperada

