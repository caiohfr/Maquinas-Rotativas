
## 3. Código Python – `trabalho5_fem_rotor_mancais.py`

import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import eig
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar


# ---------------------------------------------------------------------
# 1. Parâmetros e matrizes FEM (copiados da estrutura do MATLAB)
# ---------------------------------------------------------------------

def build_fem_matrices():
    """
    Monta Ke, Me, Ge globais, incluindo disco (Md, Gd),
    igual ao que é feito no Trabalho 5 em MATLAB.
    """
    # Dados do mancal e eixo
    Dm = 0.03     # diâmetro do mancal (m)
    R = Dm / 2.0
    Lm = 0.02     # comprimento do mancal (m)
    folga_r = 0.09e-3
    n_visc = 0.07  # viscosidade (Pa.s)

    De = 0.01     # diâmetro do eixo (m)
    Le = 0.8      # comprimento do eixo (m)
    Ste = np.pi * De**2 / 4.0  # área da seção

    # Dados do disco
    Dd = 0.09
    r = Dd / 2.0
    Ld = 0.015

    # Material
    E = 2e11
    ro = 7850.0
    v = 0.3
    g = 9.80665

    Ie = np.pi * De**4 / 64.0

    # Massa eixo e disco
    me = (np.pi * De**2 / 4.0) * Le * ro
    md = (np.pi * Dd**2 / 4.0) * Ld * ro

    # Desbalanceamento, momentos de inércia
    e = 1e-4 / md
    Id = 0.5 * md * (Dd / 2.0)**2 + md * e**2
    Ip = md * (Dd / 2.0)**2 / 2.0
    F0 = md * g / 2.0

    # Discretização FEM
    nelm = 5
    Lelm = Le / nelm

    # Matrizes locais (copiadas do MATLAB – 8x8)
    ke_local = (E * Ie / Lelm**3) * np.array([
        [12,      6*Lelm,   0,        0,       -12,     6*Lelm,   0,        0],
        [6*Lelm,  4*Lelm**2,0,        0,       -6*Lelm, 2*Lelm**2,0,        0],
        [0,       0,        12,      -6*Lelm,  0,       0,       -12,     -6*Lelm],
        [0,       0,       -6*Lelm,   4*Lelm**2,0,      0,        6*Lelm,  2*Lelm**2],
        [-12,    -6*Lelm,   0,        0,       12,     -6*Lelm,   0,        0],
        [6*Lelm,  2*Lelm**2,0,        0,      -6*Lelm,  4*Lelm**2,0,        0],
        [0,       0,       -12,       6*Lelm,  0,       0,        12,      6*Lelm],
        [0,       0,        6*Lelm,   2*Lelm**2,0,      0,        6*Lelm,  4*Lelm**2],
    ])

    me_local = (ro * Ste * Lelm / 420.0) * np.array([
        [156,   22*Lelm,   0,         0,       54,    -13*Lelm,   0,        0],
        [22*Lelm,4*Lelm**2,0,         0,      13*Lelm, -3*Lelm**2,0,        0],
        [0,      0,       156,      -22*Lelm, 0,       0,       54,      13*Lelm],
        [0,      0,      -22*Lelm,   4*Lelm**2,0,      0,      -13*Lelm, -3*Lelm**2],
        [54,    13*Lelm,   0,         0,      156,   -22*Lelm,   0,        0],
        [-13*Lelm,-3*Lelm**2,0,       0,     -22*Lelm,4*Lelm**2, 0,        0],
        [0,      0,       54,      -13*Lelm, 0,       0,      156,      22*Lelm],
        [0,      0,       13*Lelm,  -3*Lelm**2,0,     0,       22*Lelm,  4*Lelm**2],
    ])

    Ge_local = (ro * Ste * De**2 / (240.0 * Lelm)) * np.array([
        [0, 0,   36,  -3*Lelm, 0, 0,  -36, -3*Lelm],
        [0, 0,   3*Lelm, -4*Lelm**2, 0, 0, -3*Lelm,  Lelm**2],
        [-36, -3*Lelm, 0, 0,   36,  -3*Lelm, 0, 0],
        [3*Lelm, 4*Lelm**2, 0, 0, -3*Lelm, -Lelm**2, 0, 0],
        [0, 0,  -36,  3*Lelm, 0, 0,   36,  3*Lelm],
        [0, 0,   3*Lelm,  Lelm**2, 0, 0, -3*Lelm, -4*Lelm**2],
        [36, 3*Lelm, 0, 0,  -36, 3*Lelm, 0, 0],
        [3*Lelm, -Lelm**2, 0, 0, -3*Lelm, 4*Lelm**2, 0, 0],
    ])

    ndofs = (nelm + 1) * 4  # 6 nós * 4 GDL = 24
    Ke = np.zeros((ndofs, ndofs))
    Me = np.zeros((ndofs, ndofs))
    Ge = np.zeros((ndofs, ndofs))

    # Montagem global: para i = 1:4:((nelm+1)*4-7)
    for i in range(0, (nelm + 1) * 4 - 7, 4):
        s = np.arange(i, i + 8)
        Ke[np.ix_(s, s)] += ke_local
        Me[np.ix_(s, s)] += me_local
        Ge[np.ix_(s, s)] += Ge_local

    # Disco no nó "central": s1 = 9:12 (MATLAB 1-based) → [8,9,10,11]
    s1 = np.arange(8, 12)
    Gd = np.zeros((ndofs, ndofs))
    Gd_local = np.array([
        [0, 0, 0,   0],
        [0, 0, 0,  -Ip],
        [0, 0, 0,   0],
        [0, Ip, 0,  0],
    ])
    Gd[np.ix_(s1, s1)] = Gd_local
    G = Ge + Gd

    Md = np.zeros((ndofs, ndofs))
    Md_local = np.array([
        [md, 0,  0,  0],
        [0,  Id, 0,  0],
        [0,  0,  md, 0],
        [0,  0,  0,  Id],
    ])
    Md[np.ix_(s1, s1)] = Md_local

    # Amortecimento proporcional
    c_factor = 1e-4   # c = 1e-4 no código, Ce = c*Ke
    Ce = c_factor * Ke

    M = Me + Md

    params = dict(
        Dm=Dm, R=R, Lm=Lm, folga_r=folga_r, n_visc=n_visc,
        De=De, Le=Le, Ste=Ste, Dd=Dd, r=r, Ld=Ld,
        E=E, ro=ro, v=v, g=g,
        me=me, md=md, e=e, Id=Id, Ip=Ip, F0=F0,
        nelm=nelm, Lelm=Lelm
    )

    return Ke, M, G, Ce, params


# ---------------------------------------------------------------------
# 2. Mancal hidrodinâmico – epsilon(S), gama, beta, k_ij, c_ij
# ---------------------------------------------------------------------

def epsilon_residual(eps, S):
    """
    Residual da equação epsilon(e,S) do mancal curto.

    No MATLAB:
        eps(i,1) = fminsearch(@(e)epsilon(e,S(i-1,1)), eps(i-1));

    >>> TODO:
        Substituir esta função pela forma correta de epsilon(e,S)
        conforme a equação de Sommerfeld do material da disciplina.

    Por enquanto, usamos um placeholder monotônico apenas para manter
    o script executável, mas ele NÃO reproduz os valores reais.

    """
    eps_target = 1.0 - np.exp(-S)  # placeholder
    return eps - eps_target


def find_eccentricity(S, eps_initial=0.1):
    def obj(eps):
        return epsilon_residual(eps, S) ** 2

    res = minimize_scalar(obj, bounds=(1e-4, 0.99), method="bounded")
    return res.x


def compute_bearing_coeffs(omega, params):
    """
    Implementa o trecho:

        Fn = (n*Lm^3*R)/(2*folga_r^2)*omega;
        S0 = ((Lm/(2*R))^2*F0)./Fn;
        S = F0./Fn;
        eps(i) = fminsearch(@(e)epsilon(e,S(i-1)), eps(i-1));
        ...
        gama, beta, kii, cii

    omega: vetor (rad/s)
    """
    folga_r = params["folga_r"]
    Lm = params["Lm"]
    R = params["R"]
    n_visc = params["n_visc"]
    F0 = params["F0"]

    omega = np.asarray(omega).reshape(-1, 1)  # Nx1

    Fn = (n_visc * Lm**3 * R) / (2.0 * folga_r**2) * omega         # força de ref.
    S0 = ((Lm / (2.0 * R))**2 * F0) / Fn                           # Sommerfeld clássico
    S = F0 / Fn                                                    # Sommerfeld modificado

    eps_list = [0.0]
    for i in range(1, len(omega) + 1):
        S_i = S[i - 1, 0]
        eps_prev = eps_list[-1] if eps_list[-1] > 0 else 0.1
        eps_i = find_eccentricity(S_i, eps_prev)
        eps_list.append(eps_i)
    eps = np.array(eps_list[1:]).reshape(-1, 1)

    exc = eps * folga_r
    alfa = np.arctan((np.pi / 4.0) * np.sqrt(1.0 - eps**2) / eps)

    Ae = 4.0 / (np.pi**2 + (16.0 - np.pi**2) * eps**2)**1.5

    gama = np.zeros((len(eps), 4))
    gama[:, 0] = (2 * np.pi**2 + (16 - np.pi**2) * eps[:, 0]**2) * Ae[:, 0]
    gama[:, 1] = (np.pi / 4.0) * (
        np.pi**2 - 2 * np.pi**2 * eps[:, 0]**2 - (16 - np.pi**2) * eps[:, 0]**4
    ) / (eps[:, 0] * np.sqrt(1 - eps[:, 0]**2)) * Ae[:, 0]
    gama[:, 2] = -(np.pi / 4.0) * (
        np.pi**2 + (32 + np.pi**2) * eps[:, 0]**2 + (32 - 2 * np.pi**2) * eps[:, 0]**4
    ) / (eps[:, 0] * np.sqrt(1 - eps[:, 0]**2)) * Ae[:, 0]
    gama[:, 3] = (
        np.pi**2 + (32 + np.pi**2) * eps[:, 0]**2 + (32 - 2 * np.pi**2) * eps[:, 0]**4
    ) / np.sqrt(1 - eps[:, 0]**2) * Ae[:, 0]

    beta = np.zeros((len(eps), 4))
    beta[:, 0] = (np.pi / 2.0) * (np.sqrt(1 - eps[:, 0]**2) / eps[:, 0]) * (
        np.pi**2 + (2 * np.pi**2 - 16) * eps[:, 0]**2
    ) * Ae[:, 0]
    beta[:, 1] = -(2 * np.pi**2 + (4 * np.pi**2 - 32) * eps[:, 0]**2) * Ae[:, 0]
    beta[:, 2] = beta[:, 1]
    beta[:, 3] = (np.pi / 2.0) * (
        np.pi**2 + (48 - 2 * np.pi**2) * eps[:, 0]**2 + np.pi**2 * eps[:, 0]**4
    ) / (eps[:, 0] * np.sqrt(1 - eps[:, 0]**2)) * Ae[:, 0]

    kii = gama * (F0 / folga_r)
    cii = beta * (F0 / folga_r) / omega

    return eps[:, 0], exc[:, 0], alfa[:, 0], S[:, 0], kii, cii


# ---------------------------------------------------------------------
# 3. Montagem de K{i} e C{i} com mancais
# ---------------------------------------------------------------------

def build_KC_lists(Ke, Ce, kii, cii, omega):
    """
    Reproduz o laço:

        F=zeros(24,1);
        for i=1:5000
          K{i}=Ke; ... adicionar k11,k12,k21,k22 em nós 1 e 6
          C{i}=Ce; ... idem com c11,c12,c21,c22
          S = omega(i)*G + C{i};
          A{i} = [0 I; -M\K{i} -M\S];
          ...
        end

    Aqui só montamos as listas K_list e C_list para uso posterior.
    """
    ndofs = Ke.shape[0]
    K_list = []
    C_list = []

    for i in range(len(omega)):
        K_i = Ke.copy()
        C_i = Ce.copy()

        k11, k12, k21, k22 = kii[i, :]
        c11, c12, c21, c22 = cii[i, :]

        # Mancal na extremidade esquerda (nó 1 → dofs 0,1,2,3)
        K_i[0, 0] += k11
        K_i[0, 2] += k12
        K_i[2, 0] += k21
        K_i[2, 2] += k22

        C_i[0, 0] += c11
        C_i[0, 2] += c12
        C_i[2, 0] += c21
        C_i[2, 2] += c22

        # Mancal na extremidade direita (nó 6 → dofs 20,21,22,23)
        i21, i23 = 20, 22

        K_i[i21, i21] += k11
        K_i[i21, i23] += k12
        K_i[i23, i21] += k21
        K_i[i23, i23] += k22

        C_i[i21, i21] += c11
        C_i[i21, i23] += c12
        C_i[i23, i21] += c21
        C_i[i23, i23] += c22

        K_list.append(K_i)
        C_list.append(C_i)

    return K_list, C_list


# ---------------------------------------------------------------------
# 4. Diagrama de Campbell (aproximação do MATLAB)
# ---------------------------------------------------------------------

def compute_campbell(omega, M, G, K_list, C_list):
    """
    Para cada omega(i), monta A(Ω) = [[0 I], [-M^-1 K_i, -M^-1 S_i]]
    com S_i = omega(i)*G + C_i, e extrai autovalores.
    """
    nd = M.shape[0]
    Z = np.zeros((nd, nd))
    I = np.eye(nd)

    Wn = []  # lista de vetores de frequências naturais (rad/s)

    Minv = np.linalg.inv(M)

    for i in range(len(omega)):
        K_i = K_list[i]
        C_i = C_list[i]

        S_i = omega[i] * G + C_i

        A = np.block([
            [Z,              I],
            [-Minv @ K_i, -Minv @ S_i]
        ])

        eigvals, _ = eig(A)
        wn_i = np.abs(eigvals)  # magnitude de lambda
        wn_i = np.sort(wn_i)[:8]  # pega alguns modos principais
        Wn.append(wn_i)

    Wn = np.array(Wn)  # shape (len(omega), n_modes)
    return Wn


# ---------------------------------------------------------------------
# 5. Espaço de estados para resposta temporal – wcons
# ---------------------------------------------------------------------

def wcons_rhs(t, u, M, K, S, md, e, w):
    """
    Versão Python da função wcons(t,u):

        global md M k e w S
        F(9)  = md*e*w^2*cos(w*t);
        F(11) = md*e*w^2*sin(w*t);
        du = A*u + Q;

    Aqui:
        u = [q; qdot] com tamanho 48
    """
    nd = M.shape[0]
    q = u[:nd]
    qdot = u[nd:]

    F = np.zeros(nd)
    F[8] = md * e * w**2 * np.cos(w * t)   # DOF 9 (MATLAB)
    F[10] = md * e * w**2 * np.sin(w * t)  # DOF 11 (MATLAB)

    Minv = np.linalg.inv(M)

    # Sistema de 2ª ordem
    qddot = Minv @ (F - K @ q - S @ qdot)

    du = np.zeros_like(u)
    du[:nd] = qdot
    du[nd:] = qddot
    return du


# ---------------------------------------------------------------------
# 6. Função principal: junta tudo e gera alguns plots
# ---------------------------------------------------------------------

def main():
    # FEM básico + disco
    Ke, M, G, Ce, params = build_fem_matrices()

    # Malha de velocidades (rad/s)
    omega = 2.0 * np.pi * np.arange(0.1, 500.1, 0.1)

    # Mancais
    eps, exc, alfa, S_mod, kii, cii = compute_bearing_coeffs(omega, params)

    # Matrizes com mancais
    K_list, C_list = build_KC_lists(Ke, Ce, kii, cii, omega)

    # Campbell aproximado
    Wn = compute_campbell(omega, M, G, K_list, C_list)  # rad/s

    # Plot Campbell (só primeiros modos)
    plt.figure()
    for mode in range(min(4, Wn.shape[1])):
        plt.plot(omega / (2.0 * np.pi), Wn[:, mode] / (2.0 * np.pi))
    plt.xlabel("Velocidade de rotação Ω (Hz)")
    plt.ylabel("Frequência (Hz)")
    plt.title("Diagrama de Campbell (aproximado)")
    plt.grid(True)

    # Plot epsilon vs velocidade
    plt.figure()
    plt.plot(omega / (2.0 * np.pi), eps)
    plt.xlabel("Velocidade de rotação Ω (Hz)")
    plt.ylabel("Excentricidade adimensional ε")
    plt.title("Excentricidade vs velocidade")
    plt.grid(True)

    # Resposta temporal em algumas rotações (como no MATLAB)
    md = params["md"]
    e = params["e"]

    t_final = 5.0
    t_eval = np.linspace(0.0, t_final, 4000)

    for w_Hz in [15.9, 50.0, 100.0]:
        w = w_Hz * 2.0 * np.pi
        idx = int(w_Hz * 10) - 1  # omega = 2*pi*(0.1:0.1:500) → i = w_Hz/0.1
        if idx < 0 or idx >= len(omega):
            continue

        K_i = K_list[idx]
        C_i = C_list[idx]
        S_i = w * G + C_i

        u0 = np.zeros(2 * M.shape[0])
        sol = solve_ivp(
            lambda t, u: wcons_rhs(t, u, M, K_i, S_i, md, e, w),
            (0.0, t_final),
            u0,
            t_eval=t_eval,
            method="BDF",
            rtol=1e-6,
            atol=1e-9
        )

        t = sol.t
        q = sol.y[:M.shape[0], :]

        # Deslocamento vertical do disco (nó 3, DOF 9 em MATLAB → índice 8)
        y_disco = q[8, :]

        plt.figure()
        plt.plot(t, 1000.0 * y_disco)
        plt.xlabel("Tempo (s)")
        plt.ylabel("Deslocamento do disco (mm)")
        plt.title(f"Resposta temporal – w = {w_Hz:.1f} Hz")
        plt.grid(True)

    plt.show()


if __name__ == "__main__":
    main()
