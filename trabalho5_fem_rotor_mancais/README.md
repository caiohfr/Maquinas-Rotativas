# Trabalho 5 – Rotor em Elementos Finitos com Mancais Hidrodinâmicos

Este projeto implementa, em Python, o **Trabalho 5** da disciplina
IM342A – Análise de Máquinas Rotativas, cujo objetivo é modelar um rotor
flexível por **elementos finitos de viga** apoiado em **mancais
hidrodinâmicos**, e analisar seu comportamento dinâmico.

## Modelo físico

- Eixo modelado como viga de Euler–Bernoulli:
  - Discretização em `nelm = 5` elementos;
  - 4 GDL por nó (`uy, θy, uz, θz`) → 24 GDL totais;
- Disco rígido montado em um dos nós intermediários (`s1 = 9:12`);
- Matriz de massa, rigidez e giroscópica globais:
  - `M = Me + Md`
  - `K = Ke`
  - `G = Ge + Gd`
- Mancais hidrodinâmicos nas extremidades (nós 1 e 6):
  - Propriedades dependentes do **número de Sommerfeld modificado** `S*`;
  - Excentricidade adimensional `ε(S*)`;
  - Fatores γ\_{ij} e β\_{ij};
  - Rigidez `k\_{ij}` e amortecimento `c\_{ij}`.

As equações de movimento globais no domínio do tempo são:

\[
M \ddot{q} + (C + S(\Omega)) \dot{q} + K(\Omega) q = F(t)
\]

onde:

- \(q \in \mathbb{R}^{24}\) é o vetor de deslocamentos globais;
- \(M\) é a matriz de massa global;
- \(K(\Omega)\) inclui a rigidez dos mancais (k\_{ij});
- \(C(\Omega)\) inclui o amortecimento dos mancais (c\_{ij}) e do eixo;
- \(S(\Omega) = \Omega G + C(\Omega)\) é a matriz usada na formulação em espaço de estados para o Campbell;
- \(F(t)\) aplica forças de desbalanceamento no nó do disco.

## Funcionalidades implementadas

O script Python:

1. Monta as matrizes locais `Ke_local`, `Me_local`, `Ge_local` e as globais `Ke`, `Me`, `Ge`;
2. Adiciona as contribuições de disco:
   - `Md` (massa e momento de inércia do disco);
   - `Gd` (matriz giroscópica do disco);
3. Calcula:
   - Número de Sommerfeld modificado `S*`;
   - Excentricidade adimensional `ε(S*)` via minimização numérica;
   - Fatores γ\_{ij} e β\_{ij};
   - Fatores de rigidez `k\_{ij}` e amortecimento `c\_{ij}` dos mancais;
4. Para cada velocidade de rotação:
   - Monta `K{i}` e `C{i}` com os termos dos mancais nas extremidades;
   - Calcula a matriz de estado `A(Ω)` e os autovalores → **diagrama de Campbell**;
5. Simula a resposta temporal do rotor para rotações específicas usando o modelo de 2ª ordem em **espaço de estados**:
   - Equivalente ao `ode15s(@wcons, ...)` do MATLAB original;
   - Força de desbalanceamento aplicada nos GDL do disco.

## Dependências

- Python 3.10+
- `numpy`
- `scipy`
- `matplotlib`

Instalação:

```bash
pip install numpy scipy matplotlib
