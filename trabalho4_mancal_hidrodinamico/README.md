# Trabalho 4 – Rotor Laval com Mancais Hidrodinâmicos

Este projeto implementa, em Python, a análise numérica de um rotor Laval
apoiado em um mancal hidrodinâmico, baseada no código original em MATLAB/Octave.

O modelo considera:

- Eixo flexível simplesmente apoiado;
- Disco rígido no meio do vão;
- Mancal hidrodinâmico de mancal curto (short bearing) com:
  - número de Sommerfeld modificado S\*,
  - excentricidade adimensional ε(S\*),
  - coeficientes auxiliares γ\_{ij} e β\_{ij},
  - rigidez k\_{ij} e amortecimento c\_{ij} dependentes da velocidade.

O script reproduz as principais etapas:

1. Cálculo do número de Sommerfeld (clássico e modificado);
2. Cálculo da excentricidade ε em função de S\* (via minimização numérica);
3. Cálculo de:
   - excentricidade física `exc`,
   - ângulo de atitude α,
   - fatores γ\_{ij} e β\_{ij},
   - fatores de rigidez k\_{ij} e amortecimento c\_{ij};
4. Integração temporal das equações do rotor:
   - Caso 1: torque suficientemente alto para cruzar rapidamente a 1ª ressonância;
   - Caso 2 e 3: torques menores, cruzando mais lentamente;
   - Casos com velocidade angular imposta (w constante) para várias rotações;
5. Geração de gráficos:
   - ε vs. velocidade angular;
   - S\* e α vs. ε;
   - γ\_{ij} vs. ε;
   - β\_{ij} vs. ε;
   - k\_{ij} e c\_{ij} vs. velocidade;
   - respostas temporais (Uy, Uz, Ym, Zm, rotação);
   - órbitas para rotações específicas.

## Dependências

- Python 3.10+
- `numpy`
- `scipy`
- `matplotlib`

Instalação:

```bash
pip install numpy scipy matplotlib
