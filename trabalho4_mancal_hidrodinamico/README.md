# Trabalho 4 – Rotor com Mancais Hidrodinâmicos (IM342)

Este diretório contém a implementação em **Python** do Trabalho 4 de Máquinas Rotativas,  
no qual é analisado um **rotor Laval apoiado em mancais hidrodinâmicos**.

O objetivo é reproduzir, o mais fielmente possível, os resultados do relatório:

- parâmetros do mancal obtidos a partir da **equação de Reynolds** (mancal curto);
- variação da **excentricidade adimensional** com a rotação;
- número de **Sommerfeld modificado** e **ângulo de lócus**;
- coeficientes auxiliares `γik` e `βik`, rigidezes `kik` e amortecimentos `cik`;
- resposta temporal para três **torques de aceleração** (modelo 5 GDL);
- órbitas do eixo dentro do mancal para **rotações constantes** (modelo 4 GDL);
- identificação qualitativa da instabilidade **oil-whirl / oil-whip**.

---

## 1. Estrutura do diretório

- `trabalho4_mancal_hidrodinamico.py`  
  Script principal em Python. Faz **tudo**:
  - define os parâmetros do rotor e do mancal;
  - calcula excentricidade, Sommerfeld, γ, β, k, c;
  - integra as equações de movimento para 5 GDL (torque × tempo);
  - integra as equações de 4 GDL (rotação constante);
  - gera os gráficos equivalentes aos do relatório.

> Não há dependências “secretas”: é só ter Python + NumPy + Matplotlib + SciPy.

---

## 2. Como rodar

No terminal, dentro da pasta do trabalho:

```bash
python trabalho4_mancal_hidrodinamico.py
```

O script vai abrir **várias figuras** em sequência.  
Se quiser salvar, use o botão de salvar de cada janela do Matplotlib.

Se preferir trabalhar em um ambiente interativo (VS Code, Jupyter, etc.),  
basta executar o arquivo como script ou importar as funções e chamar `main()`.

---

## 3. Modelagem do mancal – parte “Reynolds + Sommerfeld”

### 3.1. Dados de entrada

Os parâmetros geométricos e físicos são os mesmos do enunciado:

- diâmetro do mancal `Dm`, comprimento `Lm`, raio `R` e folga radial `δ`;
- viscosidade do óleo `η`;
- propriedades do eixo e do disco (`De`, `Le`, `Dd`, `Ld`, `E`, `ρ`, etc.);
- desbalanceamento `e` construído de forma que `md · e = 1,5×10⁻⁴ kg·m`;
- força axial `F0 = md·g/2`, usada como carga no mancal.

Esses valores estão encapsulados em `build_parameters()`.

### 3.2. Excentricidade adimensional ε(Ω)

O programa considera um vetor de velocidades angulares:

- `omega_Hz` de **0,01 a 50 Hz** (step pequeno),
- `omega = 2π · omega_Hz` em rad/s.

Para cada rotação, é calculada a força de referência:

\[
F_\eta = \frac{\eta L^3 R}{2 \, \delta^2} \, \Omega
\]

Em seguida:

- número de Sommerfeld **modificado**:
  \[
  S^\* = \frac{F_0}{F_\eta}
  \]
- para cada valor de `S*`, é resolvida numericamente a relação analítica:

  \[
  S^\*(\varepsilon) = \frac{\pi}{2}\,\frac{\varepsilon}{(1-\varepsilon^2)^2}
  \sqrt{1-\varepsilon^2 + \left(\frac{4}{\pi}\varepsilon\right)^2}
  \]

  Usando um método de **bisseção** entre `ε = 0` e `ε ≈ 1` para obter `ε(S*)`.

O primeiro conjunto de gráficos que você deve conferir:

1. **Excentricidade ε × velocidade Ω (Hz)**  
   → Deve reproduzir a figura de excentricidade do relatório.

2. **S* e ângulo α × ε**  
   - `S*` em escala log (à esquerda);  
   - `α` (rad ou deg, dependendo da versão) à direita.  
   O ângulo α é dado por:
   \[
   \alpha = \arctan\left(\frac{\pi}{4}\frac{\sqrt{1-\varepsilon^2}}{\varepsilon}\right)
   \]

3. **Lócus do eixo** (se a versão Python tiver o gráfico polar)  
   – excentricidade vs. ângulo α, mostrando a evolução do tipo de lubrificação
   (contato, filme misto, regime hidrodinâmico), como discutido na introdução.

### 3.3. Fatores γik, βik e coeficientes do mancal

Usando as expressões analíticas da solução de Reynolds para mancais curtos,
o código calcula:

- função auxiliar:
  \[
  A(\varepsilon) = 
  \frac{4}{\left[\pi^2 + (16-\pi^2)\varepsilon^2\right]^{3/2}}
  \]
- **fatores γik(ε)** (`gama11`, `gama12`, `gama21`, `gama22`);
- **fatores βik(ε)** (`beta11`, `beta12`, `beta21`, `beta22`).

Em seguida:

- rigidezes:
  \[
  k_{ik} = \gamma_{ik} \frac{F_0}{\delta}
  \]
- amortecimentos:
  \[
  c_{ik} = \frac{\beta_{ik}}{\Omega}\frac{F_0}{\delta}
  \]

Gráficos a observar:

1. `γik × ε` em escala log  
   - procure a variação brusca de `γ₁₂` próxima de ε ≈ 0,7,
     como comentado no relatório (associada à instabilidade fluido-induzida).

2. `βik × ε` em escala log.

3. `kik × Ω`  
   – rigidezes hidrodinâmicas em função da velocidade angular.

4. `cik × Ω`  
   – amortecimentos hidrodinâmicos em função da velocidade angular.

Essas curvas devem ser comparadas com as figuras 9, 10, 11 e 12 do relatório.

---

## 4. Modelo de 5 GDL – torque variável (ode45 → solve_ivp)

Nesta parte o rotor é modelado com **5 graus de liberdade**:

- deslocamentos do disco: `uy`, `uz`;
- deslocamentos do eixo dentro do mancal: `ym`, `zm`;
- rotação `θ` do eixo.

O script monta o espaço de estados de acordo com as equações (18) e (20)  
do relatório e integra com `solve_ivp` (equivalente ao `ode45` do MATLAB).

São simulados **três cenários de torque**:

1. **Caso 1 – T = 0,010 N·m**  
   - torque suficiente para cruzar rapidamente a primeira ressonância;
   - a rotação sobe rápido, pico de deslocamento bem localizado.

2. **Caso 2 – T = 0,0045 N·m**  
   - cruzamento **lento** da região de ressonância;
   - maiores amplitudes de vibração, resposta mais “gorda”.

3. **Caso 3 – T = 0,0040 N·m**  
   - torque **insuficiente** para ultrapassar a ressonância;
   - a rotação “encalha” perto da crítica, com deslocamentos persistentes.

Os gráficos correspondentes são:

- `uy(t)` e `uz(t)` (em mm) vs. tempo,  
- rotação `ω(t)` vs. tempo,  
- deslocamentos do centro do mancal `ym(t)`, `zm(t)`.

Compare com as figuras 13, 14 e 15 do relatório.  
A frequência natural estimada em todos os casos deve estar em torno de **72 rad/s**,  
e o torque necessário para cruzar a ressonância é **menor** que no caso de mancais rígidos (Trabalho 2).

---

## 5. Modelo de 4 GDL – rotação constante e órbitas

Para estudar a instabilidade fluido-induzida, o relatório reduz o modelo para 4 GDL  
assumindo **rotação constante** (`θ(t) = Ωt`):

- deslocamentos `uy`, `uz` do disco;
- deslocamentos `ym`, `zm` do eixo dentro do mancal.

O script Python reproduz isso na função de rotação constante, usando:

- os coeficientes `kik` e `cik` **interpolados** para cada valor de Ω desejado;
- integração em regime permanente e posterior plotagem das órbitas.

São avaliadas as órbitas para:

- Ω = 30 rad/s  
- Ω = 50 rad/s  
- Ω = 72 rad/s (frequência natural)  
- Ω = 120 rad/s  
- Ω = 145 rad/s  
- Ω = 150 rad/s  

Dica ao rodar: o código costuma descartar o início da simulação e plota somente
os **últimos períodos em regime permanente** para filtrar os transientes,  
como feito manualmente no MATLAB (seleção de índices `i1`, `i2` no script original).

Os gráficos de órbita (`uy` × `uz`, em mm) correspondem às figuras 16 a 21:

- aumento de amplitude até perto de Ω ≈ 72 rad/s;
- redução da órbita após a passagem pela crítica;
- comportamento **caótico** próximo a Ω ≈ 150 rad/s,
  com deslocamentos chegando à ordem da folga radial → evidência de **oil-whip**.

---

## 6. Como usar o código para estudo

Sugestões de uso didático (segue o espírito do relatório):

1. **Refazer todas as figuras** do relatório diretamente em Python:
   - rodar o script e conferir se cada figura tem o mesmo padrão de forma e escala;
   - ajustar `xlim`, `ylim` e escalas log/linear quando necessário.

2. **Explorar a influência dos parâmetros**:
   - alterar viscosidade `n_visc`, folga radial `folga_r` ou carga `F0`;
   - observar como isso desloca as curvas de `ε`, `S*`, `kik`, `cik`  
     e como afeta a posição das órbitas em diferentes velocidades.

3. **Testar outros torques** nos casos 5 GDL:
   - ver como torque mais alto ou mais baixo altera o tempo de passagem pela ressonância  
   e a amplitude de vibração.

4. **Investigar oil-whirl / oil-whip**:
   - focar nas órbitas próximas a 2× a frequência natural;
   - comparar com as discussões do relatório sobre a linha 0,5Ω e o valor típico ~0,48Ω.

---

## 7. Possíveis extensões futuras

Algumas ideias para expandir o trabalho (caso você queira levar isso para portfólio):

- implementar outros modelos de mancais (longos, anisotrópicos, desgaste);
- incluir não-linearidades no filme de óleo e comparar com o modelo linear atual;
- criar uma interface gráfica simples (por exemplo em **Streamlit**)  
  para alterar parâmetros e visualizar em tempo real os gráficos de órbita e resposta temporal.

---

Se algo fugir muito das figuras do relatório (escala, forma da curva, etc.),  
o primeiro passo é conferir:

- se a unidade do eixo (Hz × rad/s × rpm) está coerente;
- se os coeficientes `kik` e `cik` estão sendo interpolados na mesma variável (Ω em rad/s);
- e se os transientes foram devidamente descartados nas órbitas.
