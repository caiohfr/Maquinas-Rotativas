# Máquinas Rotativas – Rotor Dynamics Lab (Mestrado)

Este repositório reúne e traduz para **Python** os cinco trabalhos da disciplina
de Máquinas Rotativas do mestrado. Ele funciona como um pequeno laboratório de
dinâmica de rotores, incluindo desde modelos massa–mola até modelagem por
elementos finitos com mancais hidrodinâmicos.

O objetivo é duplo:

1) **Documentar e reproduzir fielmente** os trabalhos acadêmicos.  
2) **Preparar terreno** para uma futura evolução em direção a uma toolbox
   modular, orientada a objetos e com ferramentas avançadas de controle.

---

# 📚 Conteúdo Acadêmico e Referências Bibliográficas

Os modelos implementados se baseiam em formulações clássicas de dinâmica de rotores,
mancais hidrodinâmicos e elementos finitos. Principais referências:

### **Livros fundamentais**
- Childs, Dara. *Turbomachinery Rotordynamics*. Wiley.  
- Lalanne, Christian; Ferraris, Guy. *Rotordynamics Prediction in Engineering*. Wiley.  
- Genta, Giancarlo. *Dynamics of Rotating Systems*. Springer.  
- Rao, J. S. *Rotor Dynamics*. New Age International.  
- Muszynska, Agnes. *Rotordynamics*. CRC Press.  
- Bently, Donald; Hatch, Charles. *Fundamentals of Rotating Machinery Diagnostics*. Bently Pressurized Bearing Co.

### **Artigos / Clássicos**
- Sommerfeld, A. *"Zur Theorie der hydrodynamischen Schmierung"*, 1904.  
- Reynolds, O. *"On the Theory of Lubrication"*, 1886.  
- Lund, J. W. *"A Review of Rotor-Dynamic Analysis"*, 1980.

### **Modelagem FEM**
- Zienkiewicz & Taylor. *The Finite Element Method*.

Essas referências sustentam diretamente o que está implementado nos Trabalhos 3–5.

---

# 📂 Estrutura do Repositório

```
maquinas-rotativas/
  README.md
  requirements.txt

  trabalho1_massa_mola_3gdl/
  trabalho2_laval_rotor/
  trabalho3_rotor_disco_descentrado/
  trabalho4_mancal_hidrodinamico/
  trabalho5_fem_rotor_mancais/
```

Cada pasta possui:

- `README.md` próprio  
- script Python completo  
- (opcional) códigos originais MATLAB/Octave  

---

# ▶️ Como Rodar

## 1. Instale dependências

```bash
pip install -r requirements.txt
```

## 2. Execute qualquer trabalho

```bash
cd trabalho4_mancal_hidrodinamico
python trabalho4_mancal_hidrodinamico.py
```

---

# ⚠️ Nota Importante — Função ε(S*)

Nos Trabalhos 4 e 5 existe uma função fundamental:

\[
arepsilon(S^*) \quad 	ext{(excentricidade adimensional do mancal curto)}
\]

No MATLAB original ela é definida em um arquivo externo (`epsilon.m`), não incluso
no anexo. Portanto, no código Python:

- A estrutura numérica está toda implementada;  
- A função `epsilon_residual()` possui um **placeholder**;  
- Basta inserir a equação exata da disciplina para obter resultados físicos.

---

# 🚀 Roadmap da Expansão Futura  
*(a versão premium do projeto que você pode construir quando quiser)*

Aqui está um plano **realista e profissional** para transformar este repositório
numa **toolbox completa de rotodinâmica**.

---

## **1. Versão 1 — Modularização (curto prazo)**  
**Objetivo:** transformar scripts isolados em módulos reutilizáveis.

- Criar estrutura `rotor_dynamics/`
- Extrair:
  - `LavalModel`
  - `HydroBearingModel`
  - `FEMRotor`
  - `OrbitTools`, `CampbellTools`
- Criar interface simples:
  ```python
  from rotor_dynamics import LavalRotor
  ```

---

## **2. Versão 2 — Orientação a Objetos (médio prazo)**  
**Objetivo:** transformar o código em um framework formal.

Classes principais:

### 🌪️ `LavalRotor`
- estados: θ, uy, uz, ...
- métodos: `simulate()`, `orbit()`

### 🔧 `HydroBearing`
- métodos: `stiffness(w)`, `damping(w)`, `sommerfeld(w)`

### 🧱 `FEMRotor`
- montagem M, C, K, G
- inserção de discos e mancais
- análise modal

---

## **3. Versão 3 — Controle Ativo (diferencial absurdo)**  
Inserir:

- PID anti-órbita  
- AMB simplificado (Active Magnetic Bearing)  
- controle feedforward  
- controle de passagem de ressonância  

Isso adiciona modernidade e vira *portfolio gold*.

---

## **4. Versão 4 — Análises Avançadas**  
- Campbell refinado  
- FRF (magnitude/fase)  
- Bode/Nyquist  
- mapas de estabilidade (root locus com Ω)  
- animação de modos

---

## **5. Versão 5 — Interface Streamlit**  
Interface gráfica com sliders:

- massa, rigidez, viscosidade  
- folga radial  
- velocidade  
- modos  
- força de desbalanceamento  
- órbitas animadas  
- diagramas em tempo real

---

# 🌟 Resultado Esperado

Ao seguir este roadmap, você terá:

### ✔ Uma toolbox própria de dinâmica de rotores  
### ✔ Uma interface intuitiva de engenharia  
### ✔ Um diferencial profissional (real) para indústria automotiva e powertrain  
### ✔ Um laboratório completo para estudo e pesquisa

---

# 📜 Licença  
Definir (MIT, BSD, ou uso acadêmico).

---
