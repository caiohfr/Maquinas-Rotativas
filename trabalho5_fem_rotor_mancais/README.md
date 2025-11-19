# 📘 README – Trabalho 5: Análise de um Rotor com Mancais Hidrodinâmicos usando Elementos Finitos (FEM)
### Disciplina: IM342 – Análise de Máquinas Rotativas  
### Aluno: *Caio Henrique Ferreira Rocha*   
### Professor: Tiago Henrique Machado  
### Universidade Estadual de Campinas – UNICAMP  

---

# 📌 Resumo

Este projeto implementa, em **Python**, todo o desenvolvimento teórico e numérico apresentado no **Trabalho 5 de Máquinas Rotativas**. O objetivo é analisar um rotor apoiado sobre mancais hidrodinâmicos por meio do **método dos elementos finitos (FEM)**, incluindo:

- Montagem das matrizes globais de massa, rigidez, amortecimento e giroscópica;  
- Cálculo dos coeficientes dinâmicos dos mancais via Sommerfeld modificado;  
- Diagrama de Campbell do sistema completo;  
- FRFs dos mancais e do disco;  
- Respostas temporais para rotações constantes;  
- Formas modais e órbitas 3D.

O código Python é uma **reprodução fiel** do MATLAB fornecido no PDF, porém reorganizado e documentado.

---

# 🧩 Introdução

Rotores são elementos presentes em turbinas, motores, compressores, máquinas industriais e sistemas automotivos. A análise dinâmica permite identificar:

- Frequências críticas;  
- Modos forward e backward;  
- Instabilidade fluido-induzida (oil-whirl e oil-whip);  
- Interação rotor–mancal;  
- Efeitos giroscópicos;  
- Respostas vibratórias sob desbalanceamento.

O objetivo final é construir o **Diagrama de Campbell**, ferramenta essencial na análise de máquinas rotativas.

---

# 🎯 Objetivos

Seguindo exatamente o PDF do Trabalho 5, este projeto implementa:

1. Matrizes FEM do eixo (K, M, G);  
2. Disco com massa e giroscopia;  
3. Coeficientes do mancal (ε, α, γ, β, k, c);  
4. Montagem de K(Ω), C(Ω), G(Ω);  
5. Solução dos autovalores (Campbell);  
6. FRFs de mancais e disco;  
7. Simulação temporal via espaço de estados;  
8. Modos de vibrar e órbitas 3D.

---

# 🧠 Desenvolvimento Teórico

## 1) Equação Dinâmica

\[
M \ddot{x} + (C + \Omega G)\dot{x} + Kx = F(t)
\]

## 2) FEM – Viga Euler-Bernoulli

Cada nó possui 4 DOFs:

\[
[y, 	heta_z, z, 	heta_y]
\]

Com 6 nós → 24 DOFs totais.

### Matriz de Rigidez (PDF Eq. 5)
\[
K_e = rac{EI}{L^3}
egin{bmatrix}
12 & 6L & \dots
\end{bmatrix}
\]

### Matriz de Massa Consistente (PDF Eq. 6)

### Matriz Giroscópica (PDF Eq. 7)

Todas elas implementadas em:

```python
Ke, Me, Ge, G, Ce, M = montar_matrizes_rotor(params)
```

---

## 3) Disco do Rotor

No nó 3 são adicionadas:

- massa;  
- inércia polar;  
- matriz giroscópica local.

---

## 4) Mancais Hidrodinâmicos

### 4.1 Excentricidade via Sommerfeld Modificado

\[
S^*(arepsilon) =
rac{\pi}{2} rac{arepsilon}{(1-arepsilon^2)^2}
\sqrt{1 - arepsilon^2 + \left(rac{4arepsilon}{\pi}
ight)^2}
\]

Resolvida via **bisseção estável** (mais robusto que o MATLAB):

```python
eps = resolver_epsilon(S)
```

### 4.2 Cálculo de γ, β, k e c

\[
k_{ik} = \gamma_{ik} rac{F_0}{c_r}
\qquad
c_{ik} = rac{eta_{ik} F_0}{c_r \Omega}
\]

Gerando as curvas:

- ε × Ω  
- S* × ε  
- α × ε  
- γik × ε  
- βik × ε  
- kik × Ω  
- cik × Ω  

---

# 📊 Diagrama de Campbell

A matriz de estado:

\[
A = egin{bmatrix}
0 & I \
-M^{-1}K & -M^{-1}(C + \Omega G)
\end{bmatrix}
\]

Os autovalores fornecem:

- Frequências naturais  
- Amortecimentos  
- Modos forward/backward

Plotado como:

```python
plt.plot(omega/(2*np.pi), Wn[:, :24] / (2*np.pi))
```

Inclui as linhas teóricas:  
- ω = Ω  
- ω = 2Ω  

---

# 🔊 FRFs

\[
H(j\omega) =
\left[
-M\omega^2 + j\omega(C+\Omega G) + K

ight]^{-1}
\]

FRFs plotadas para:

- Mancal 1  
- Mancal 2  
- Disco  

---

# ⏱️ Resposta Temporal (w = constante)

Força de desbalanceamento:

\[
F(t) = m_d e \Omega^2
egin{bmatrix}
\cos(\Omega t) \
\sin(\Omega t)
\end{bmatrix}
\]

Equações integradas com `solve_ivp(BDF)`:

```python
t, x = simular_wconst(...)
```

---

# 🌀 Modos de Vibrar / Órbitas 3D

A partir dos deslocamentos:

\[
r = \sqrt{y^2 + z^2}
\]

Plotam-se as curvas orbitais para cada nó do rotor.

---

# 📂 Estrutura do Repositório

```
📁 trabalho5_rotor_mancais/
│
├── trabalho5_fem_rotor_mancais.py
├── README.md
├── figs/
│   ├── eccentricidade.png
│   ├── gamma_beta.png
│   ├── campbell.png
│   ├── frf_mancal1.png
│   ├── frf_mancal2.png
│   ├── frf_disco.png
│   ├── respostas_temporais.png
│   └── orbitas_3d.png
└── requirements.txt
```

---

# 📦 Requirements

```
numpy
scipy
matplotlib
```

---

# ▶️ Como Executar

```bash
python trabalho5_fem_rotor_mancais.py
```

---

# 📚 Referências (do seu PDF)

1. Notas de aula IM342 – Máquinas Rotativas.  
2. Kramer, E. *Dynamics of Rotors and Foundations*.  
3. Childs, D. *Turbomachinery Rotordynamics*.  
4. Machado, T. H., *UNiCAMP – Hidrodinâmica de Mancais*.  
5. Código MATLAB original do Trabalho 5.

---

# ✨ Observação Final

Este README foi adaptado **diretamente do seu PDF oficial**, mantendo rigor teórico e fidelidade acadêmica, e estruturado para uso em **portfólio profissional** (GitHub, entrevistas, CV técnico).
