# Trabalho 3 – Análise no Domínio da Frequência de um Rotor com Disco Descentrado

Este repositório contém a implementação em Python do **Trabalho 3** da disciplina  
**IM342 – Análise de Máquinas Rotativas (FEM/UNICAMP)**.

O objetivo é analisar, **no domínio da frequência**, um rotor com disco descentrado em um modelo de **4 graus de liberdade**:

- Translações: \(u_y, u_z\)
- Rotações: \(\varphi_y, \varphi_z\)

O código:

1. Monta as matrizes \(M\), \(K\) e \(G(\Omega)\) do rotor descentrado;
2. Calcula os **autovalores em função da rotação** para:
   - Caso **não amortecido** → frequências naturais (diagrama de Campbell);
   - Caso **amortecido** → frequências naturais amortecidas e fatores de amortecimento \(\zeta(\Omega)\);
3. Calcula a **matriz FRF** \(H(j\omega)\) para:
   - Caso não amortecido;
   - Caso amortecido;
4. Gera os mesmos tipos de gráficos discutidos no relatório:
   - Diagramas de Campbell (full + zoom);
   - Diagrama de \(\zeta\) vs \(\Omega\);
   - FRFs (amplitude em dB e fase em graus) para ambos os casos.

---

## 1. Modelo físico e parâmetros utilizados

### 1.1 Geometria e material

Os dados do enunciado (e do relatório) são:

- Diâmetro do eixo: \(D_e = 0{,}012\ \text{m}\)
- Comprimento do eixo: \(L_e = 0{,}85\ \text{m}\)
- Diâmetro do disco: \(D_d = 0{,}10\ \text{m}\)
- Comprimento do disco: \(L_d = 0{,}01\ \text{m}\)
- Módulo de elasticidade: \(E = 2 \cdot 10^{11}\ \text{Pa}\)
- Densidade: \(\rho = 7850\ \text{kg/m}^3\)
- Divisão do eixo: \(a = \frac{2}{3}L_e,\quad b = \frac{1}{3}L_e\)

Momento de inércia de área do eixo:

\[
I_e = \frac{\pi D_e^4}{64}
\]

Massa do disco e raios:

\[
m_d = \rho \frac{\pi D_d^2}{4} L_d,\quad r = \frac{D_d}{2}
\]

Momentos de inércia do disco:

\[
I_p = \frac{1}{2} m_d r^2
\]
\[
I_d = \frac{1}{4} m_d r^2 + \frac{1}{12} m_d L_d^2
\]

### 1.2 Rigidezes equivalentes

As rigidezes são obtidas como em Resistência dos Materiais, reproduzindo a Tabela do relatório:

\[
K = 3 E I_e \frac{a^3 + b^3}{a^3 b^3}
\]
\[
k = 3 E I_e \frac{L_e (a - b)}{a^2 b^2}
\]
\[
K_r = 3 E I_e \frac{L_e}{a b}
\]

Matriz de rigidez:

\[
K =
\begin{bmatrix}
K & 0 & 0 & -k \\
0 & K & k & 0 \\
0 & k & K_r & 0 \\
-k & 0 & 0 & K_r
\end{bmatrix}
\]

### 1.3 Matrizes \(M\), \(C\) e \(G(\Omega)\)

Matriz de massa (4 GDL: \(u_y, u_z, \varphi_y, \varphi_z\)):

\[
M = \mathrm{diag}(m_d, m_d, I_d, I_d)
\]

Caso não amortecido:

\[
C = 0
\]

Caso amortecido (amortecimento proporcional):

\[
C = \alpha K,\quad \alpha = 2 \times 10^{-4}
\]

Matriz giroscópica (dependente da rotação \(\Omega\) em rad/s):

\[
G(\Omega) =
\begin{bmatrix}
0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 \\
0 & 0 & 0 & I_p \Omega \\
0 & 0 & -I_p \Omega & 0
\end{bmatrix}
\]

---

## 2. Equações de movimento e formulação em espaço de estados

As equações gerais do rotor são:

### 2.1 Caso não amortecido

\[
M \ddot{x} + G(\Omega) \dot{x} + K x = 0
\]

### 2.2 Caso amortecido

\[
M \ddot{x} + [C + G(\Omega)] \dot{x} + K x = 0
\]

com:

\[
x =
\begin{bmatrix}
u_y \\ u_z \\ \varphi_y \\ \varphi_z
\end{bmatrix}
\]

A forma em espaço de estados é:

\[
y(t) =
\begin{bmatrix}
x \\ \dot{x}
\end{bmatrix},\quad
\dot{y}(t) = A(\Omega)\, y(t)
\]

onde:

\[
A(\Omega) =
\begin{bmatrix}
0 & I \\
-M^{-1} K & -M^{-1} (C + G(\Omega))
\end{bmatrix}
\]

Para cada valor de \(\Omega\) (varrendo o espectro), calculam-se os autovalores de \(A(\Omega)\):

- Caso não amortecido → \(\lambda = \pm j \omega_{ni}\)
- Caso amortecido → \(\lambda_i = -\zeta_i \omega_{ni} \pm j \omega_{di}\)

O fator de amortecimento é obtido em função dos autovalores:

\[
\zeta_i = \left[ \left( \frac{\Im(\lambda_i)}{\Re(\lambda_i)} \right)^2 + 1 \right]^{-1/2}
\]

---

## 3. Matriz FRF

A matriz de impedância no domínio de Laplace é:

\[
Z(s) = M s^2 + (C + G(\Omega)) s + K
\]

Para a FRF (substituindo \(s = j\omega\)):

\[
H(j\omega) = Z(j\omega)^{-1}
= \left[-M \omega^2 + j\omega (C + G(\Omega)) + K \right]^{-1}
\]

No código, considera-se **excitação síncrona**:

- \(\omega = \Omega\)
- Para cada velocidade de rotação \(\Omega\), monta-se \(Z(j\Omega)\) e computa-se \(H(j\Omega)\)

Os termos de interesse da matriz FRF são:

- \(H_{11}, H_{22}, H_{33}, H_{44}\) – respostas diretas;
- \(H_{12}, H_{13}, H_{14}, H_{23}, H_{24}, H_{34}\) – respostas cruzadas.

> Importante: aqui os índices 1,2,3,4 referem-se aos GDL \((u_y, u_z, \varphi_y, \varphi_z)\),  
> não aos modos do diagrama de Campbell.

---

## 4. Estrutura do código Python

Arquivo principal: **`trabalho3_rotor_disco_descentrado.py`**

Principais blocos:

- **Definição de parâmetros e matrizes**
  - Geometria, material e constantes \(M\), \(K\), \(C\), \(G(\Omega)\).

- **`G_matrix(Omega_rad)`**
  - Monta \(G(\Omega)\) a partir de \(\Omega\) em rad/s.

- **`state_matrix(M, K, C, G)`**
  - Monta a matriz de estados \(A(\Omega)\).

- **`campbell_diagram(wmax_Hz, dph)`**
  - Calcula autovalores para o **caso não amortecido** (\(C = 0\));
  - Retorna matriz `omega1` com frequências naturais (rad/s) para cada \(\Omega\).

- **`campbell_with_damping(wmax_Hz, dph, alpha)`**
  - Usa \(C = \alpha K\);
  - Retorna `omega2` (frequências amortecidas) e `zeta` (fatores de amortecimento).

- **`compute_frf_vs_speed(wmax_Hz, dph, C)`**
  - Calcula FRFs para excitação síncrona \(\omega = \Omega\);
  - Retorna um dicionário com `H11`, `H22`, `H33`, `H44`, `H12`, `H13`, `H14`, `H23`, `H24`, `H34`.

- **Funções de plot**
  - `plot_campbell_pair(...)` → Diagrama de Campbell (janela full 0–1000 Hz e zoom 0–40 Hz);
  - `plot_zeta(...)` → \(\zeta(\Omega)\) vs \(\Omega\);
  - `plot_frf(...)` → Amplitude (dB) e fase (graus) vs \(\Omega\).

- **`main()`**
  - Define `wmax_Hz = 1000.0`, `dph = 100`, `alpha = 2e-4`;
  - Gera:
    - Campbell sem amortecimento;
    - Campbell com amortecimento;
    - Gráfico de \(\zeta(\Omega)\);
    - FRF sem amortecimento;
    - FRF com amortecimento.

---

## 5. Como executar

### 5.1 Dependências

- Python 3.10+
- `numpy`
- `matplotlib`

Instalação:

```bash
pip install numpy matplotlib
```

### 5.2 Execução

Na pasta do projeto:

```bash
python trabalho3_rotor_disco_descentrado.py
```

Será gerado um conjunto de figuras:

1. **Diagrama de Campbell – não amortecido (full + zoom)**
2. **Diagrama de Campbell – amortecido (full + zoom)**
3. **Fator de amortecimento \(\zeta_i(\Omega)\)**
4. **FRF – caso não amortecido (amplitude em dB + fase)**
5. **FRF – caso amortecido (amplitude em dB + fase)**

---

## 6. Como reproduzir os resultados e conclusões do relatório

### 6.1 Frequências naturais em \(\Omega = 0\)

Ao inspecionar `omega1[:, 0]/(2π)` e `omega2[:, 0]/(2π)` (se desejar, pode imprimir no código), você deve encontrar valores muito próximos aos do relatório:

- **Caso não amortecido:**
  - \(\omega_{30}(0) = \omega_{40}(0) \approx 28{,}74\ \text{Hz}\)
  - \(\omega_{10}(0) = \omega_{20}(0) \approx 458{,}43\ \text{Hz}\)

- **Caso amortecido:**
  - \(\omega_{30}(0) = \omega_{40}(0) \approx 28{,}73\ \text{Hz}\)
  - \(\omega_{10}(0) = \omega_{20}(0) \approx 439{,}00\ \text{Hz}\)

Essas duplicidades em \(\Omega = 0\) refletem a **isotropia** do eixo (propriedades idênticas em y e z).

---

### 6.2 Diagramas de Campbell – precessão e velocidades críticas

Nos gráficos de Campbell (não amortecido e amortecido):

- Em \(\Omega = 0\), as curvas 1–2 e 3–4 partem de frequências coincidentes;
- As curvas \(1 \times \Omega\) e \(2 \times \Omega\) (linhas tracejadas) representam:
  - **1×Ω** → excitações típicas de **desbalanceamento**;
  - **2×Ω** → excitações típicas de **desalinhamento**.

Observe os cruzamentos entre essas linhas e os modos:

- A linha \(2×\Omega\) cruza a frequência 2 (modo backward) em torno de **158,2 Hz**;
- A linha \(1×\Omega\) cruza a frequência 2 em torno de **265,3 Hz**.

Esses cruzamentos representam **velocidades críticas** associadas à mudança de precessão:

- O rotor inicia em **precessão forward** (mesmo sentido da rotação);
- Ao cruzar a frequência 2 (modo backward), ocorre mudança para **precessão backward**;
- Devido à proximidade entre as frequências 3 e 4, as mudanças no cruzamento dessas curvas não são claramente percebidas graficamente (como comentado no texto).

No caso amortecido, o comportamento geral é o mesmo, com pequenas variações numéricas:

- As velocidades críticas das linhas 1×Ω e 2×Ω permanecem próximas de 158,2 Hz e 259–265 Hz;
- As frequências naturais amortecidas são ligeiramente menores, devido ao amortecimento pequeno \(\alpha = 2 \times 10^{-4}\).

---

### 6.3 Diagrama de \(\zeta(\Omega)\) – auto-centragem

No gráfico de \(\zeta_i(\Omega)\):

- **\(\zeta_1\) e \(\zeta_2\)** (modos associados principalmente às **translações** \(u_y, u_z\)) mantêm-se praticamente constantes com o aumento da rotação;
- **\(\zeta_3\) e \(\zeta_4\)** (modos associados principalmente às **rotações** \(\varphi_y, \varphi_z\)) **decrescem com \(\Omega\)**.

Essa queda de \(\zeta_3\) e \(\zeta_4\) está relacionada ao fenômeno de **auto-centragem**:

- Com o aumento da rotação, a órbita de precessão do eixo diminui;
- As forças elásticas e de amortecimento associadas às rotações ficam menores;
- O amortecimento “efetivo” desses modos decai.

---

### 6.4 FRFs – caso não amortecido

Na figura de FRFs **sem amortecimento**:

- Os picos de amplitude (em dB) aparecem próximos de:
  - **~28 Hz** → primeira frequência natural (modos 3 e 4 no Campbell);
  - **~265 Hz** → velocidade crítica associada à curva **1×Ω**;
  - **~458 Hz** → frequência natural associada principalmente aos modos rotacionais \(\omega_1\) e \(\omega_2\).

Observe:

- Os picos de ~458 Hz aparecem apenas nas curvas envolvendo índices 3 e 4 (relacionados a \(\varphi_y, \varphi_z\)), confirmando que:
  - Essa frequência está ligada aos **modos rotacionais**;
- Mudanças bruscas de fase em torno das frequências naturais reforçam a identificação dos modos.

---

### 6.5 FRFs – caso amortecido

Na figura de FRFs **com amortecimento**:

- Os mesmos picos aparecem, mas:
  - A região da velocidade crítica (~258–265 Hz) apresenta **transição mais suave** na amplitude;
  - Os picos são “arredondados” pelo amortecimento.

Isso evidencia as conclusões do relatório:

- O amortecimento (mesmo pequeno) **suaviza a resposta** na região de ressonância;
- A localização das frequências naturais amortecidas permanece muito próxima dos valores não amortecidos, devido ao baixo nível de \(C\).

---

## 7. Síntese – o que o usuário deve conseguir observar

Rodando o script e analisando os gráficos:

1. **Isotropia e duplicidade de modos**:
   - Frequências naturais coincidentes em \(\Omega = 0\) para modos de translação (1–2) e rotação (3–4).

2. **Efeito giroscópico** nos diagramas de Campbell:
   - Separação das frequências 1 e 2 com o aumento de \(\Omega\);
   - Cruzamentos com as linhas 1×Ω e 2×Ω marcando velocidades críticas (~158 Hz e ~265 Hz).

3. **Mudança de precessão**:
   - Associação dos cruzamentos com transição de **forward** para **backward whirl**.

4. **Auto-centragem**:
   - Decaimento de \(\zeta_3\) e \(\zeta_4\) nas curvas de amortecimento;
   - Ligação com a redução da órbita de precessão em altas rotações.

5. **Coerência Campbell × FRF**:
   - Picos de FRF (28 Hz, ~265 Hz, ~458 Hz) alinhados com as frequências do diagrama de Campbell;
   - Picos rotacionais aparecendo apenas em termos com índices 3 e 4.

Seguindo este README, qualquer pessoa com Python instalado consegue:

- Rodar o código;
- Gerar os mesmos tipos de diagramas do relatório;
- E interpretar os gráficos sob a mesma ótica física: precessão, velocidades críticas, auto-centragem e efeito do amortecimento.
