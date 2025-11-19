# Trabalho 1 – Sistema Massa–Mola–Amortecedor de 3 GDL

Este repositório contém a implementação em Python do Trabalho 1 da disciplina  
**IM342 – Análise de Máquinas Rotativas**, cujo objetivo é estudar um sistema
massa–mola–amortecedor com **3 graus de liberdade (3 GDL)** por meio de:

1. Formulação matricial e **análise modal** em espaço de estados;
2. Cálculo das **Funções de Resposta em Frequência (FRF)**;
3. Obtenção da **resposta em frequência de deslocamento**;
4. **Resposta temporal** (livre e forçada) com integração numérica (Runge–Kutta).

O código em Python foi organizado para reproduzir os mesmos resultados e permitir
a mesma interpretação física apresentada no relatório original do trabalho.

---

## 1. Modelo físico e matemático

O sistema é um conjunto de três massas ligadas por molas e amortecedores
viscosos, em deslocamento horizontal:

- Massas: \(m_1 = 2\ 	ext{kg}\), \(m_2 = 5\ 	ext{kg}\), \(m_3 = 1\ 	ext{kg}\);
- Rigidezes:
  - \(k_1 = 4000\ 	ext{N/m}\),
  - \(k_2 = 8000\ 	ext{N/m}\),
  - \(k_3 = 6000\ 	ext{N/m}\);
- Amortecimentos:
  - \(c_1 = 100\ 	ext{N·s/m}\),
  - \(c_2 = 200\ 	ext{N·s/m}\),
  - \(c_3 = 150\ 	extN·s/m}\).

As equações de movimento podem ser escritas na forma matricial:

\[
M \ddot{x} + C \dot{x} + K x = F(t),
\]

com:

- \(x = [x_1, x_2, x_3]^T\) (deslocamentos),
- \(M\), \(C\), \(K\) montadas a partir de \(m_i\), \(c_i\), \(k_i\),
- \(F(t)\) representando a força externa.

No caso forçado, a excitação é harmônica aplicada na massa 2:

\[
F(t) = [0,\; F_0 \sin(\omega t),\; 0]^T, \quad F_0 = 50\ \text{N}.
\]

---

## 2. Estrutura do código

O código principal (por exemplo, `trabalho1_massa_mola_3gdl.py`) está organizado em:

- `build_mck()`: monta as matrizes \(M\), \(C\), \(K\);
- `modal_analysis(...)`: monta a matriz dinâmica \(A\), calcula autovalores/autovetores,
  frequências naturais \(\omega_n\), fatores de amortecimento \(arepsilon\) e modos relevantes;
- `compute_frf(...)`: calcula \(H(\omega)\) a partir da matriz de rigidez dinâmica
  \[
  Z(\omega) = -M\omega^2 + jC\omega + K, \quad H(\omega) = Z(\omega)^{-1};
  \]
- `compute_steady_state_response(H, F)`: calcula
  \[
  X(\omega) = H(\omega)\,F
  \]
  para obter as amplitudes e fases dos deslocamentos de cada GDL;
- `vib_rhs(...)` / `simulate_time_response(...)`: montam e integram a resposta no tempo
  (Runge–Kutta via `solve_ivp`);
- Funções de **plot**:
  - FRFs \(H_{ij}(\omega)\),
  - Amplitude e fase de \(X(\omega)\),
  - Deslocamentos \(x_i(t)\) e velocidades \(\dot{x}_i(t)\).

---

## 3. Como executar

### 3.1 Dependências

- Python 3.10+
- `numpy`
- `scipy`
- `matplotlib`

Instale com:

```bash
pip install numpy scipy matplotlib
```

### 3.2 Execução básica

Na pasta do projeto:

```bash
python trabalho1_massa_mola_3gdl.py
```

O script irá:

1. Montar \(M\), \(C\), \(K\);
2. Calcular parâmetros modais;
3. Calcular e plotar as FRFs \(H_{ij}(\omega)\);
4. Calcular e plotar a resposta em frequência \(X(\omega)\) (amplitude e fase);
5. Integrar a resposta no tempo para uma excitação harmônica e condições iniciais definidas.

---

## 4. Roteiro para reproduzir as ANÁLISES do relatório

Esta seção orienta passo a passo o que observar nos gráficos gerados para que o
usuário chegue às **mesmas conclusões físicas** discutidas no trabalho.

### 4.1 Parâmetros modais (espaço de estados)

1. A partir da matriz dinâmica \(A\), o código obtém os autovalores \(\lambda_i\).
2. A partir de \(\lambda_i\), são calculados:
   - Frequências naturais \(\omega_{n,i}\),
   - Fatores de amortecimento \(arepsilon_i\).

**O que conferir:**

- Você deve encontrar duas frequências naturais principais em torno de:
  - \(\omega_{n,1} pprox 19{,}5\ 	ext{rad/s}\),
  - \(\omega_{n,2} pprox 77{,}5\ 	ext{rad/s}\),
  com fatores de amortecimento diferentes (um mais baixo, outro muito próximo
  de amortecimento crítico).
- Apenas os autovetores associados a **modos subamortecidos** (\(arepsilon_i < 1\))
  são relevantes para análise modal.
- Compare essas frequências com as posições dos picos nas FRFs \(H(\omega)\):
  os picos mais relevantes devem ocorrer próximos a \(\omega_{n,1}\) e \(\omega_{n,2}\).

---

### 4.2 Funções de Resposta em Frequência \(H(\omega)\)

O código varre um intervalo de frequências \(\omega\) e, para cada ponto:

\[
Z(\omega) = -M\omega^2 + jC\omega + K,\quad
H(\omega) = Z(\omega)^{-1}.
\]

São extraídas as magnitudes:

- \(H_{11}, H_{22}, H_{33}\) (respostas nos próprios GDL),
- \(H_{12}, H_{23}, H_{13}\) (acoplamento entre GDLs).

**O que observar (ligado ao texto do relatório):**

- Os maiores picos de \(|H_{ij}(\omega)|\) aparecem **próximo à primeira frequência natural**,
  indicando que o sistema responde fortemente nessa região.
- As funções ligadas à **massa 3 (GDL 3)** tendem a apresentar **maior amplitude**,
  pois essa massa é a mais leve e está associada a uma rigidez/amortecimento menores.
- A massa 2 também mostra amplitudes elevadas porque é exatamente o ponto
  de aplicação da força.
- Observe a reciprocidade:
  \(|H_{12}(\omega)|\) e \(|H_{21}(\omega)|\), \(|H_{13}(\omega)|\) e \(|H_{31}(\omega)|\), etc.,
  devem ser compatíveis (sistema linear satisfaz a reciprocidade de Maxwell).

---

### 4.3 Resposta em frequência de deslocamento \(X(\omega)\)

Com a excitação:

\[
F =
egin{bmatrix}
0 \\
50 \\
0
\end{bmatrix},
\quad
X(\omega) = H(\omega)\,F,
\]

o código calcula, para cada GDL:

- **Amplitude** \(|X_i(\omega)|\) (plotada em **mm**, via \(1000\cdot|X_i|\));
- **Fase** \(ngle X_i(\omega)\) em graus.

**O que observar para reproduzir as conclusões do trabalho:**

1. **Amplitude (mm):**
   - As curvas de \(x_1, x_2, x_3\) apresentam um pico pronunciado próximo à **primeira frequência natural**.
   - A massa 3 (GDL 3) tende a ter a **maior amplitude**, devido à menor massa, rigidez e amortecimento.
   - Para frequências muito acima da primeira natural, as amplitudes caem bastante,
     mesmo próximo à segunda frequência natural, evidenciando o efeito do amortecimento.

2. **Fase (graus):**
   - Abaixo da primeira frequência natural, as três massas tendem a oscilar aproximadamente em fase.
   - Após cruzar a primeira frequência natural, surge uma **defasagem** mais evidente,
     especialmente da massa 2 em relação às massas 1 e 3.
   - Essa defasagem aumenta com a frequência, o que é consistente com o comportamento
     de sistemas com amortecimento.

Essas observações batem com os comentários do relatório, onde se destaca:

- A maior amplitude na massa 3;
- A mudança de fase da massa 2 após a primeira frequência natural.

---

### 4.4 Resposta temporal (livre e forçada)

A resposta no tempo é obtida integrando:

\[
M \ddot{x} + C \dot{x} + K x = F(t),
\]

com a mesma estrutura em espaço de estados usada no trabalho original.

#### 4.4.1 Vibração livre

Para reproduzir os gráficos de vibração livre:

1. No código, defina:
   - \(F_0 = 0\) (sem força externa),
   - Condições iniciais, por exemplo:
     \[
     x(0) = [0{,}01,\ 0{,}005,\ 0{,}005],\quad
     \dot{x}(0) = [0,\ 1,\ 0].
     \]
2. Rode a simulação e plote \(x_i(t)\) e \(\dot{x}_i(t)\).

**O que observar:**

- Decaimento da amplitude ao longo do tempo (efeito do amortecimento).
- Frequência dominante próxima à primeira frequência natural identificada na análise modal.
- Diferenças de amplitude entre as massas, coerentes com os modos de vibração.

#### 4.4.2 Vibração forçada – variação de frequência

Para reproduzir a sequência de gráficos do relatório (ω = 10, 19,5, 30, 50, 77,5 rad/s):

1. Mantenha as mesmas condições iniciais utilizadas na vibração livre;
2. Para cada simulação, escolha um valor de `w_exc` (10, 19.5, 30, 50, 77.5 rad/s)
   e um valor não nulo de `F0` (por exemplo, 50 N);
3. Execute a simulação e plote:

   - Deslocamentos \(x_1(t), x_2(t), x_3(t)\)
   - Velocidades \(\dot{x}_1(t), \dot{x}_2(t), \dot{x}_3(t)\)

**O que verificar para bater com as conclusões:**

- **ω = 10 rad/s (abaixo da 1ª natural):**
  - Resposta relativamente pequena;
  - Sistema não está em ressonância, a energia é mais facilmente “dissipada” pelo amortecimento.

- **ω ≈ 19,5 rad/s (1ª frequência natural):**
  - Maior amplitude de vibração, principalmente na massa 3 e na massa 2;
  - Tempo maior para o sistema estabilizar;
  - Perfil em regime permanente claramente maior que nos outros casos.

- **ω = 30 e 50 rad/s (acima da 1ª natural, abaixo da 2ª):**
  - Amplitudes menores que no caso em ressonância;
  - Fica evidente que, conforme a frequência aumenta, **o amortecimento se torna mais relevante**,
    reduzindo a resposta.

- **ω ≈ 77,5 rad/s (2ª frequência natural):**
  - Apesar de ser uma frequência natural, o amortecimento mais alto neste modo
    faz com que a resposta não atinja amplitudes tão elevadas quanto na primeira natural;
  - Isso está em linha com a conclusão de que **apenas próximo à primeira frequência natural
    há presença significativa de vibrações**.

---

## 5. Como conectar os gráficos às conclusões do trabalho

Ao final, o usuário deve conseguir chegar às mesmas conclusões principais:

1. **Coerência entre domínio do tempo e da frequência**  
   - As amplitudes observadas em \(x(t)\) nas frequências excitadas
     são consistentes com os picos de \(|X(\omega)|\) e \(|H(\omega)|\).

2. **Predominância da primeira frequência natural**  
   - É na vizinhança da primeira frequência natural que surgem as maiores amplitudes;
   - Acima desta, a resposta tende a ser bem menor, mesmo na região da segunda natural.

3. **Papel do amortecimento**  
   - Modos de alta frequência apresentam fatores de amortecimento mais elevados,
     o que explica a menor amplificação na segunda frequência natural;
   - No tempo, isso aparece como uma resposta com menor amplitude e decaimento mais rápido.

4. **Diferença de comportamento entre as massas**  
   - A massa 3 apresenta as maiores amplitudes de deslocamento;
   - A massa 2 é significativamente excitada por ser o ponto de aplicação da força;
   - A curva de fase mostra a defasagem crescente da massa 2 em relação às massas 1 e 3
     à medida que a frequência aumenta.

Seguindo esse roteiro, qualquer usuário que rode o código em Python consegue
não apenas reproduzir os gráficos, mas também **interpretá-los da mesma forma**
que foi feito no relatório original.
