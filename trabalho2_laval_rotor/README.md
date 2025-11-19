# Trabalho 2 – Análise no Domínio do Tempo de um Rotor Laval

Este repositório contém a implementação em Python do **Trabalho 2** da disciplina  
**IM342 – Análise de Máquinas Rotativas (FEM/UNICAMP)**.

O objetivo é estudar, no domínio do tempo, um **rotor de Laval/Jeffcott**:

- Eixo flexível, isotrópico e de massa desprezível;
- Disco rígido, de massa concentrada, montado no meio do vão;
- Centro de massa do disco desbalanceado em relação ao centro geométrico;
- Apoios rígidos nas extremidades;
- Deslocamentos laterais no plano **Y–Z** com três graus de liberdade:  
  \(u_y, u_z\) (translações) e \(\theta\) (rotação).

A partir do equacionamento via Lagrange e da implementação em espaço de estados, o código:

1. Calcula a **resposta temporal** para três torques diferentes (aceleração rápida, lenta e insuficiente para atravessar a 1ª ressonância);
2. Gera **órbitas** dos pontos W (centro geométrico) e G (centro de massa) para velocidades constantes;
3. Compara amplitudes obtidas numericamente com a **fórmula de Kramer**, reproduzindo os resultados e conclusões do relatório.

---

## 1. Modelo físico e parâmetros utilizados

Os dados geométricos e de material são:

- Diâmetro do eixo: \(D_e = 0{,}012\ \text{m}\)
- Comprimento do eixo: \(L_e = 0{,}85\ \text{m}\)
- Diâmetro do disco: \(D_d = 0{,}10\ \text{m}\)
- Comprimento do disco: \(L_d = 0{,}015\ \text{m}\)
- Módulo de elasticidade: \(E = 2 \cdot 10^{11}\ \text{Pa}\)
- Densidade: \(\rho = 7850\ \text{kg/m}^3\)
- Gravidade: \(g = 9{,}81\ \text{m/s}^2\)

A partir destes dados, o script calcula:

- Massa do disco:
  \[
  m_d = \rho \frac{\pi D_d^2}{4} L_d \approx 0{,}92481\ \text{kg}
  \]
- Excentricidade:
  \[
  m_d e = 1{,}5 \times 10^{-3}\ \Rightarrow\ e \approx 0{,}001622\ \text{m}
  \]
- Momento de inércia de área do eixo:
  \[
  I_e = \frac{\pi D_e^4}{64} \approx 1{,}02 \times 10^{-9}\ \text{m}^4
  \]
- Rigidez equivalente do eixo:
  \[
  k = \frac{48 E I_e}{L_e^3} \approx 15911\ \text{N/m}
  \]
- Amortecimento viscoso:
  \[
  c = 10^{-3} k \approx 15{,}911\ \text{N·s/m}
  \]
- Momento de inércia polar do disco:
  \[
  I_d = \frac{1}{2} m_d \left(\frac{D_d}{2}\right)^2 + m_d e^2
  \approx 1{,}156 \times 10^{-3}\ \text{kg·m}^2
  \]

Para a análise modal plana (direções y e z), as matrizes equivalentes são:

- Massa:
  \[
  M = \begin{bmatrix}
  0{,}92481 & 0 \\
  0 & 0{,}92481
  \end{bmatrix}
  \]
- Rigidez:
  \[
  K = \begin{bmatrix}
  15911 & 0 \\
  0 & 15911
  \end{bmatrix}
  \]
- Amortecimento:
  \[
  C = \begin{bmatrix}
  15{,}911 & 0 \\
  0 & 15{,}911
  \end{bmatrix}
  \]

---

## 2. Equações de movimento

### 2.1 Sistema completo (3 GDL: \(u_y, u_z, \theta\))

Usando Lagrange com as considerações do relatório, obtêm-se:

\[
m_d \ddot{u}_y + c \dot{u}_y + k u_y = m_d e \left( \ddot{\theta} \sin\theta + \dot{\theta}^2 \cos\theta \right)
\]

\[
m_d \ddot{u}_z + c \dot{u}_z + k u_z = m_d e \left( -\ddot{\theta} \cos\theta + \dot{\theta}^2 \sin\theta \right) - m_d g
\]

\[
I_d \ddot{\theta} = T + k u_z e \cos\theta - k u_y e \sin\theta - c \dot{u}_y e \sin\theta + c \dot{u}_z e \cos\theta
\]

O vetor de estado usado no código (caso torque) é:

\[
y =
\begin{bmatrix}
\theta \\ u_y \\ u_z \\ \dot{\theta} \\ \dot{u}_y \\ \dot{u}_z
\end{bmatrix},
\quad
\dot{y} = f(t,y,T)
\]

implementado na função `laval_rhs(t, y, T)`.

### 2.2 Velocidade angular constante

Quando a velocidade angular é imposta como constante \(\dot{\theta} = \Omega\), \(\ddot{\theta} = 0\), as equações se desacoplam:

\[
m_d \ddot{u}_y + c \dot{u}_y + k u_y = m_d e \Omega^2 \cos(\Omega t)
\]

\[
m_d \ddot{u}_z + c \dot{u}_z + k u_z = m_d e \Omega^2 \sin(\Omega t) - m_d g
\]

Com o vetor de estado reduzido:

\[
u =
\begin{bmatrix}
u_y \\ u_z \\ \dot{u}_y \\ \dot{u}_z
\end{bmatrix},
\quad
\dot{u} = f(t,u,\Omega)
\]

implementado em `vibwconst_rhs(t, u, w)`.

---

## 3. Estrutura do código Python

Arquivo principal: **`trabalho2_laval_rotor.py`**

Principais blocos:

- **Definição de parâmetros**  
  Calcula `me`, `md`, `e`, `Ie`, `k`, `c`, `Id` exatamente como no enunciado e no relatório.

- **`laval_rhs(t, y, T)`**  
  Implementa o sistema completo com torque aplicado (equações 8, 9 e 10 do relatório), com:
  ```text
  y = [theta, uy, uz, theta_dot, uy_dot, uz_dot]
  ```

- **`vibwconst_rhs(t, u, w)`**  
  Implementa o caso de velocidade angular constante (equações 11 e 12), com:
  ```text
  u = [uy, uz, uy_dot, uz_dot]
  ```

- **Simulações de torque:**
  - `simulate_case_1()` → T = 2.0 N·m, t_end = 1.0 s (cruza rapidamente a 1ª ressonância);
  - `simulate_case_2()` → T = 0.3 N·m, t_end = 4.0 s (cruza lentamente);
  - `simulate_case_3()` → T = 0.25 N·m, t_end = 4.0 s (torque limite / insuficiente).

- **Simulação de órbita em velocidade constante:**
  - `simulate_orbit_constant_speed(w, t_final=2.0, n_points=4000, transient_fraction=0.6)`  
    Integra até `t_final` e descarta os **60% iniciais** para remover o transiente (igual ao trecho `i=6*L/10:L` do código Octave).

- **Funções de plot:**
  - `plot_case(...)` → gráficos \(u_y(t)\) e \(u_z(t)\) em mm, com \(\dot{\theta}(t)\) no eixo secundário;
  - `plot_orbit(...)` → órbitas no plano Y–Z dos pontos W (centro do rotor) e G (centro de massa).

- **Comparação com Kramer:**
  - `compare_amplitudes(uy_orbit, w)` → calcula:
    \[
    r_{\text{teo}} = \frac{m_d e \Omega^2}{\sqrt{(k - m_d \Omega^2)^2 + (c \Omega)^2}}
    \]
    e compara com a amplitude numérica máxima de \(|u_y|\) em regime permanente.

---

## 4. Como executar

### 4.1 Dependências

- Python 3.10+
- `numpy`
- `scipy`
- `matplotlib`

Instalação:

```bash
pip install numpy scipy matplotlib
```

### 4.2 Execução padrão

Na pasta do projeto:

```bash
python trabalho2_laval_rotor.py
```

Isso irá:

1. Rodar os **três casos de torque** (T1, T2, T3);
2. Plotar os históricos \(u_y(t)\), \(u_z(t)\) (em mm) com \(\dot{\theta}(t)\);
3. Simular uma órbita em velocidade constante (por padrão, `w_const = 1000.0 rad/s`);
4. Plotar as órbitas de W e G;
5. Imprimir no terminal a amplitude teórica (Kramer) e a amplitude numérica para esse caso.

---

## 5. Como reproduzir os resultados e conclusões do relatório

### 5.1 Parâmetros modais (opcional)

No relatório, a análise modal via espaço de estados (matriz A) fornece:

- Autovalores:  
  \[
  \lambda = -8{,}60 \pm 130{,}89\,i
  \]
- Frequência natural:
  \[
  \omega_n \approx 131{,}17\ \text{rad/s}
  \]
- Fator de amortecimento:
  \[
  \varepsilon \approx 0{,}065584
  \]

Esses valores podem ser reobtidos em Python montando:

\[
A = \begin{bmatrix}
0 & I \\
-M^{-1}K & -M^{-1}C
\end{bmatrix}
\]

e usando `numpy.linalg.eig(A)`.

**Interpretação (relacionar com o resto):**

- Só há **uma frequência natural relevante** (e mesma propriedade em y e z → isotropia);
- Os picos de vibração e as amplitudes máximas nas órbitas ocorrem na vizinhança desta \(\omega_n\).

---

### 5.2 Casos de torque – resposta temporal

#### 5.2.1 Caso 1 – \(T_1 = 2\ \text{N·m},\ t = 1\ \text{s}\)

- Aceleração angular muito rápida;
- Nos gráficos:
  - \(\dot{\theta}(t)\) cresce quase linear e atravessa a região de \(\omega_n \approx 131\ \text{rad/s}\) rapidamente;
  - \(u_y(t)\) e \(u_z(t)\) em mm mostram apenas um pequeno aumento na região crítica.

**Conclusão esperada (como no relatório):**

> A rotação passa pela frequência natural de maneira praticamente despercebida,  
> ocorre um aumento de amplitude naquela região, mas ela é **significativamente menor**  
> que nos casos de torque mais baixo (aproximadamente metade da amplitude máxima  
> observada nos outros casos).

#### 5.2.2 Caso 2 – \(T_2 = 0{,}3\ \text{N·m},\ t = 4\ \text{s}\)

- Aceleração angular mais lenta, o rotor fica mais tempo perto de \(\omega_n\);
- Nos gráficos:
  - \(\dot{\theta}(t)\) sobe mais devagar;
  - Observa-se uma “barriga” mais nítida nas amplitudes de \(u_y(t)\) e \(u_z(t)\)
    enquanto o rotor permanece na região de ressonância;
  - Após sair da vizinhança crítica, as amplitudes se reduzem e se aproximam
    das do Caso 1.

**Conclusão:**

> Com torque moderado, há tempo suficiente para que a energia de vibração se acumule,  
> resultando em amplitudes bem maiores na região crítica. Esse caso evidencia  
> com clareza o efeito da ressonância em um rotor desbalanceado.

#### 5.2.3 Caso 3 – \(T_3 = 0{,}25\ \text{N·m},\ t = 4\ \text{s}\)

- Torque apenas ligeiramente menor que no Caso 2, mas suficiente para **não**
  permitir que a rotação ultrapasse a 1ª crítica;
- Nos gráficos:
  - A curva de \(\dot{\theta}(t)\) se “achata” e não atinge velocidades acima de \(\omega_n\);
  - As amplitudes de \(u_y(t)\) e \(u_z(t)\) ficam **altas durante praticamente toda a simulação**;
  - Fica evidente que o rotor permanece “preso” na região de ressonância.

**Conclusão de limite de torque (como no relatório):**

> Entre **0,25 N·m e 0,3 N·m** está o limite de torque para que o rotor consiga  
> superar a região de ressonância. Torques abaixo disso resultam em operação  
> praticamente “travada” na frequência crítica, com amplitudes elevadas.

---

### 5.3 Órbitas W e G – auto-centragem e precessão

No relatório, quatro velocidades constantes foram analisadas:  
\(\Omega = 90,\ 130,\ 150,\ 1000\ \text{rad/s}\).

Para reproduzir cada uma no Python:

1. Edite o valor de `w_const` na função `main()`:

   ```python
   w_const = 90.0   # ou 130.0, 150.0, 1000.0
   ```

2. (Opcional) Ajuste `t_final` e `transient_fraction` em `simulate_orbit_constant_speed` para garantir que:
   - haja tempo suficiente para atingir regime permanente,
   - o corte de 60% inicial elimine a maior parte do transiente.

3. Execute o script e observe o gráfico de órbita (W e G) + resultados impressos.

#### 5.3.1 Abaixo da frequência natural – \(\Omega = 90\ \text{rad/s}\)

- Teoria de Kramer:
  - Amplitude: \(r \approx 1{,}4225\ \text{mm}\)
  - Deslocamento estático: \(y_s \approx 0{,}57\ \text{mm}\) (devido ao peso)
- Numérico:
  - \(r_W \approx 1{,}4267\ \text{mm}\)
- Observação visual:
  - A órbita de G é **maior** que a de W → “barriga para fora” no movimento de precessão;
  - Condiz com o comportamento típico **abaixo** da frequência crítica.

#### 5.3.2 Próximo da frequência natural – \(\Omega \approx 130\ \text{rad/s}\)

- Teoria:
  - \(r \approx 12{,}143\ \text{mm}\)
- Numérico:
  - \(r_W \approx 12{,}126\ \text{mm}\)
- Observação visual:
  - As órbitas de W e G são **praticamente coincidentes**;
  - A amplitude é máxima (região de ressonância);
  - Evidencia o fato de que, na crítica, o desbalanceamento e a flexibilidade fazem
    o rotor vibrar como se todo o sistema se movesse junto.

#### 5.3.3 Acima da frequência natural – \(\Omega = 150\ \text{rad/s}\)

- Teoria:
  - \(r \approx 2{,}8143\ \text{mm}\)
- Numérico:
  - \(r_W \approx 2{,}8152\ \text{mm}\)
- Observação visual:
  - A órbita de G se torna **menor** que a de W;
  - A amplitude diminui em relação ao caso na crítica;
  - Confirma o comportamento típico **acima** da velocidade crítica
    (auto-centragem começando a aparecer).

#### 5.3.4 Velocidade muito alta – \(\Omega = 1000\ \text{rad/s}\)

- Teoria:
  - \(r \approx 1{,}6501\ \text{mm}\)
- Numérico:
  - \(r_W \approx 1{,}6504\ \text{mm}\)
- Observação visual:
  - A órbita de G se aproxima de um ponto praticamente em O (apoio) →  
    **fenômeno de auto-centragem** bem característico;
  - A órbita de W ainda possui um raio finito;
  - Ruídos numéricos aparecem com mais intensidade (como descrito no relatório),
    mas a amplitude externa bate muito bem com a fórmula de Kramer.

---

## 6. Síntese das conclusões que o usuário deve conseguir ver

Rodando o código e seguindo as instruções acima, o usuário consegue reproduzir:

1. **Dinâmica da passagem pela 1ª velocidade crítica:**
   - Torques altos → passagem rápida → pequena amplificação;
   - Torques médios → amplificação grande mas finita;
   - Torques abaixo de ~0,3 N·m → rotor “preso” na ressonância.

2. **Consistência numérico × analítico (Kramer):**
   - Amplitudes calculadas pela fórmula de Kramer e obtidas via integração numérica
     são muito próximas para todas as rotações analisadas.

3. **Fenômenos de precessão e auto-centragem:**
   - Abaixo da crítica: órbita de G maior que W (barriga para fora);
   - Na crítica: órbitas praticamente coincidentes;
   - Acima da crítica: órbita de G menor que W;
   - Em velocidades muito altas: G se aproxima do centro → auto-centragem clara.

4. **Influência da gravidade:**
   - O deslocamento médio em Z (\(y_s \approx 0{,}57\ \text{mm}\)) está sempre presente,
     puxando as órbitas para baixo, exatamente como mostrado nos valores numéricos
     do relatório.

Seguindo esse README, qualquer pessoa com o Python instalado consegue não só
**rodar o código**, mas também **interpretar os resultados** da mesma forma que foi feito
no trabalho escrito.
