import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import eig
from scipy.integrate import solve_ivp


# ============================================================
# 1) PARÂMETROS E DADOS – MESMOS DO MATLAB
# ============================================================

def build_parameters():
    # Mancal
    Dm = 0.03          # diâmetro do mancal (m)
    R = Dm / 2.0       # raio do mancal (m)
    Lm = 0.02          # comprimento do mancal (m)
    folga_r = 0.09e-3  # folga radial (m)
    n_visc = 0.07      # viscosidade absoluta (Pa.s)

    # Eixo / disco
    De = 0.01          # diâmetro do eixo (m)
    Le = 0.8           # comprimento do eixo (m)
    Ste = np.pi * De**2 / 4.0  # área seção transversal
    Dd = 0.09          # diâmetro do disco (m)
    r_disc = Dd / 2.0
    Ld = 0.015         # comprimento do disco (m)

    E = 2.0e11         # módulo de elasticidade (Pa)
    ro = 7850.0        # densidade (kg/m^3)
    v = 0.3            # coef. Poisson
    g = 9.80665        # gravidade (m/s²)

    Ie = np.pi * De**4 / 64.0  # momento de inércia de área do eixo

    c_struc = 1e-4            # amortecimento estrutural (mesmo escalar "c")

    me = Ste * Le * ro        # massa do eixo
    md = (np.pi * Dd**2 / 4.0) * Ld * ro   # massa do disco

    e = 1e-4 / md             # desbalanceamento (m)

    Id = 0.5 * md * (Dd / 2.0)**2 + md * e**2   # inércia polar (disco)
    Ip = md * (Dd / 2.0)**2 / 2.0               # inércia polar (outro uso)

    F0 = md * g / 2.0          # carga no mancal

    return dict(
        Dm=Dm, R=R, Lm=Lm, folga_r=folga_r, n_visc=n_visc,
        De=De, Le=Le, Ste=Ste, Dd=Dd, r_disc=r_disc, Ld=Ld,
        E=E, ro=ro, v=v, g=g,
        Ie=Ie, c_struc=c_struc,
        me=me, md=md, e=e, Id=Id, Ip=Ip, F0=F0
    )


# ============================================================
# 2) FUNÇÃO DA EXCENTRICIDADE – SEPARADA
# ============================================================

def sommerfeld_modificado_epsilon(eps: float) -> float:
    """
    S*(ε) da tua função epsilon.m original:

        erro = abs(S - ((pi/2)*(e/((1-e^2)^2))*sqrt(1-(e^2)+((4*e/pi)^2))));

    Aqui isolamos o S*(ε) (sem o 'abs' e sem S).
    """
    if eps <= 0.0 or eps >= 0.999999:
        return np.inf

    num = eps
    den = (1.0 - eps**2)**2
    inside = 1.0 - eps**2 + (4.0 * eps / np.pi)**2
    return (np.pi / 2.0) * (num / den) * np.sqrt(inside)


def resolver_epsilon(S_vec, tol=1e-8, max_iter=60):
    """
    Resolve S*(ε) = S, 0 < ε < 1, por bisseção.
    Usa EXATAMENTE a forma da tua função de epsilon,
    mas invertendo ao invés de usar fminsearch.
    """
    S_vec = np.asarray(S_vec)
    eps = np.zeros_like(S_vec)

    for i, S in enumerate(S_vec):
        if S <= 0:
            eps[i] = 0.0
            continue

        a, b = 1e-6, 0.999999
        for _ in range(max_iter):
            m = 0.5 * (a + b)
            Sm = sommerfeld_modificado_epsilon(m)
            if np.abs(Sm - S) < tol:
                break
            if Sm > S:
                b = m
            else:
                a = m
        eps[i] = 0.5 * (a + b)

    return eps


# ============================================================
# 3) MATRIZES DE ELEMENTO FINITO – MESMAS DO MATLAB
# ============================================================

def montar_matrizes_rotor(params, nelm=5):
    """
    Monta Ke, Me, Ge (24x24) exatamente como no MATLAB.
    DOFs por nó: [uy, thetaz, uz, thetay].
    """
    E = params["E"]
    Ie = params["Ie"]
    Ste = params["Ste"]
    ro = params["ro"]
    De = params["De"]
    Le = params["Le"]

    Lelm = Le / nelm

    # ke_local (8x8) – copiado da tua ke_local
    L = Lelm
    ke_local = (E * Ie / L**3) * np.array([
        [12,      6*L,     0,       0,     -12,     6*L,     0,       0],
        [6*L,  4*L**2,     0,       0,    -6*L,  2*L**2,     0,       0],
        [0,        0,    12,   -6*L,       0,       0,   -12,   -6*L],
        [0,        0,  -6*L,  4*L**2,      0,       0,   6*L,  2*L**2],
        [-12,   -6*L,     0,       0,     12,    -6*L,     0,       0],
        [6*L,  2*L**2,    0,       0,   -6*L,   4*L**2,    0,       0],
        [0,        0,   -12,   6*L,       0,       0,    12,    6*L],
        [0,        0,  -6*L,  2*L**2,     0,       0,   6*L,  4*L**2],
    ])

    # me_local (8x8)
    me_local = (ro * Ste * L / 420.0) * np.array([
        [156,       22*L,       0,        0,      54,    -13*L,      0,        0],
        [22*L,   4*L**2,       0,        0,    13*L,   -3*L**2,      0,        0],
        [0,          0,     156,    -22*L,      0,        0,     54,    13*L],
        [0,          0,   -22*L,   4*L**2,      0,        0,   -13*L,  -3*L**2],
        [54,       13*L,       0,        0,     156,   -22*L,      0,        0],
        [-13*L, -3*L**2,       0,        0,  -22*L,   4*L**2,      0,        0],
        [0,          0,      54,   -13*L,      0,        0,    156,    22*L],
        [0,          0,   13*L,  -3*L**2,      0,        0,    22*L,  4*L**2],
    ])

    # Ge_local (8x8)
    Ge_local = (ro * Ste * De**2 / (240.0 * L)) * np.array([
        [0,      0,     36,   -3*L,   0,     0,    -36,   -3*L],
        [0,      0,   3*L, -4*L**2,   0,     0,   -3*L,  L**2],
        [-36, -3*L,     0,      0,   36,   -3*L,    0,      0],
        [3*L, 4*L**2,   0,      0,  -3*L, -L**2,    0,      0],
        [0,      0,   -36,   3*L,    0,     0,    36,   3*L],
        [0,      0,   3*L,  L**2,    0,     0,   -3*L, -4*L**2],
        [36,   3*L,    0,      0,  -36,   3*L,    0,      0],
        [3*L, -L**2,   0,      0,  -3*L, 4*L**2,   0,      0],
    ])

    ndof = (nelm + 1) * 4  # 6 nós * 4 GDL = 24
    Ke = np.zeros((ndof, ndof))
    Me = np.zeros((ndof, ndof))
    Ge = np.zeros((ndof, ndof))

    # montagem global – loop i=1:4:((nelm+1)*4-7)
    for i in range(0, ndof - 7, 4):
        s = np.arange(i, i + 8)
        Ke[np.ix_(s, s)] += ke_local
        Me[np.ix_(s, s)] += me_local
        Ge[np.ix_(s, s)] += Ge_local

    # Disco no nó 3 → s1 = 9:12 (índices MATLAB) → 8..11 (0-based)
    md = params["md"]
    Id = params["Id"]
    Ip = params["Ip"]
    s1 = np.arange(8, 12)  # DOFs 9..12 em 0-based

    Gd = np.zeros_like(Ge)
    Gd_local = np.array([
        [0.0,  0.0, 0.0,  0.0],
        [0.0,  0.0, 0.0, -Ip],
        [0.0,  0.0, 0.0,  0.0],
        [0.0,  Ip,  0.0,  0.0],
    ])
    Gd[np.ix_(s1, s1)] += Gd_local
    G = Ge + Gd

    Md = np.zeros_like(Me)
    Md_local = np.array([
        [md, 0.0, 0.0, 0.0],
        [0.0, Id, 0.0, 0.0],
        [0.0, 0.0, md, 0.0],
        [0.0, 0.0, 0.0, Id],
    ])
    Md[np.ix_(s1, s1)] += Md_local

    c_struc = params["c_struc"]
    Ce = c_struc * Ke
    M = Me + Md

    return Ke, Me, Ge, G, Ce, M


# ============================================================
# 4) COEFICIENTES DO MANCAL (ε, γ, β, k, c)
# ============================================================

def coeficientes_mancal(params, omega):
    """
    omega: vetor [rad/s], mesmo que 'omega' do MATLAB (2*pi*0.1:0.1:500).
    """
    R = params["R"]
    Lm = params["Lm"]
    folga_r = params["folga_r"]
    n_visc = params["n_visc"]
    F0 = params["F0"]

    omega = np.asarray(omega)
    # Força de referência Fn
    Fn = (n_visc * Lm**3 * R) / (2.0 * folga_r**2) * omega *1/( 2.0 * np.pi )

    # Fator de Sommerfeld clássico e modificado
    S0 = ((Lm / (2.0 * R))**2 * F0) / Fn      # não usado depois, mas mantido
    S = F0 / Fn

    # Excentricidade ε(S)
    eps = resolver_epsilon(S)

    # Ângulo alfa, Ae, gama, beta – mesmíssimas fórmulas
    alfa = np.arctan((np.pi / 4.0) * (np.sqrt(1.0 - eps**2) / eps))
    Ae = 4.0 / (np.pi**2 + (16.0 - np.pi**2) * eps**2)**1.5

    gama11 = (2.0 * np.pi**2 + (16.0 - np.pi**2) * eps**2) * Ae
    gama12 = ((np.pi / 4.0) *
              (np.pi**2 - 2.0 * np.pi**2 * eps**2 - (16.0 - np.pi**2) * eps**4) /
              (eps * np.sqrt(1.0 - eps**2)) * Ae)
    gama21 = -((np.pi / 4.0) *
               (np.pi**2 + (32.0 + np.pi**2) * eps**2 +
                (32.0 - 2.0 * np.pi**2) * eps**4) /
               (eps * np.sqrt(1.0 - eps**2)) * Ae)
    gama22 = ((np.pi**2 + (32.0 + np.pi**2) * eps**2 +
              (32.0 - 2.0 * np.pi**2) * eps**4) /
              np.sqrt(1.0 - eps**2) * Ae)

    gama = np.column_stack((gama11, gama12, gama21, gama22))

    beta11 = ((np.pi / 2.0) * (np.sqrt(1.0 - eps**2) / eps) *
              (np.pi**2 + (2.0 * np.pi**2 - 16.0) * eps**2) * Ae)
    beta12 = -(2.0 * np.pi**2 + (4.0 * np.pi**2 - 32.0) * eps**2) * Ae
    beta21 = beta12.copy()
    beta22 = ((np.pi / 2.0) *
              (np.pi**2 + (48.0 - 2.0 * np.pi**2) * eps**2 + np.pi**2 * eps**4) /
              (eps * np.sqrt(1.0 - eps**2)) * Ae)

    beta = np.column_stack((beta11, beta12, beta21, beta22))

    # fatores de rigidez e amortecimento do mancal
    kii = gama * (F0 / folga_r)
    cii = beta * (F0 / folga_r) / omega[:, None]

    return eps, S, S0, alfa, gama, beta, kii, cii


# ============================================================
# 5) CAMPBELL + FRF (Hw, Hmancal1, Hdisco, Hmancal2)
# ============================================================

def mag2db(x):
    return 20.0 * np.log10(np.abs(x))


def construir_matrizes_por_velocidade(Ke, Ce, M, G, omega, kii, cii):
    """
    Constroi vetores de K[i], C[i] (i=0..len(omega)-1) incluindo mancais.
    """
    n_omega = omega.size
    ndof = Ke.shape[0]

    K_list = np.zeros((n_omega, ndof, ndof), dtype=float)
    C_list = np.zeros((n_omega, ndof, ndof), dtype=float)

    for i in range(n_omega):
        K_i = Ke.copy()
        C_i = Ce.copy()

        k11, k12, k21, k22 = kii[i, :]
        c11, c12, c21, c22 = cii[i, :]

        # Mancal 1 – DOFs 1,3 -> idx 0,2
        K_i[0, 0] += k11
        K_i[0, 2] += k12
        K_i[2, 0] += k21
        K_i[2, 2] += k22

        C_i[0, 0] += c11
        C_i[0, 2] += c12
        C_i[2, 0] += c21
        C_i[2, 2] += c22

        # Mancal 2 – DOFs 21,23 -> idx 20,22
        K_i[20, 20] += k11
        K_i[20, 22] += k12
        K_i[22, 20] += k21
        K_i[22, 22] += k22

        C_i[20, 20] += c11
        C_i[20, 22] += c12
        C_i[22, 20] += c21
        C_i[22, 22] += c22

        K_list[i] = K_i
        C_list[i] = C_i

    return K_list, C_list

def calcular_campbell_e_frfs(params, M, G, K_list, C_list, omega):
    """
    Replica o loop i=1:5000 do MATLAB:
      - A{i}, autovalores → Wn, qsi, Wd
      - Zw, Hw, Xw
      - Hmancal1, Hdisco, Hmancal2

    IMPORTANTE: o sistema de estado tem 48 autovalores,
    mas a gente só guarda os 24 primeiros (igual Wn(:,1:24) no .m).
    """
    ndof = M.shape[0]          # 24
    n_omega = omega.size

    # Só 24 modos guardados (como no MATLAB: 1:24)
    Wn = np.zeros((n_omega, ndof))
    qsi = np.zeros((n_omega, ndof))
    Wd = np.zeros((n_omega, ndof))

    Hw = np.zeros((n_omega, ndof, ndof), dtype=complex)
    Xw = np.zeros((ndof, n_omega), dtype=complex)

    md = params["md"]
    e = params["e"]

    for i in range(n_omega):
        K_i = K_list[i]
        C_i = C_list[i]
        w_i = omega[i]

        S_i = w_i * G + C_i

        # Matriz de estado A{i} = [0 I; -M\K -M\S]
        Z = np.zeros_like(M)
        I = np.eye(ndof)
        Minv = np.linalg.inv(M)
        A_i = np.block([
            [Z,              I],
            [-Minv @ K_i, -Minv @ S_i],
        ])

        eigvals, _ = eig(A_i)
        eigvals = np.asarray(eigvals).ravel()   # 48 autovalores

        # Frequência natural "bruta"
        wn = np.sqrt(np.real(eigvals)**2 + np.imag(eigvals)**2)  # 48

        # Fator de amortecimento para cada autovalor
        qsi_i = -np.real(eigvals) / np.maximum(wn, 1e-12)        # 48

        # Ordena por frequência e PEGA SÓ OS 24 PRIMEIROS
        idx_sort = np.argsort(wn)
        wn_sorted = wn[idx_sort]
        qsi_sorted = qsi_i[idx_sort]

        wn_sel = wn_sorted[:ndof]         # 24
        qsi_sel = qsi_sorted[:ndof]       # 24

        Wn[i, :] = wn_sel
        qsi[i, :] = qsi_sel

        # Frequência amortecida Wd
        Wd_i = wn_sel * np.sqrt(np.maximum(1.0 - qsi_sel**2, 0.0))
        Wd[i, :] = Wd_i

        # ---------- FRF / Hw ----------
        Zw = -M * w_i**2 + 1j * w_i * S_i + K_i
        Zw_inv = np.linalg.inv(Zw)
        Hw[i, :, :] = Zw_inv

        F = np.zeros(ndof, dtype=complex)
        # F(9) e F(11) no MATLAB → índices 8 e 10 (0-based)
        F[8] = md * e * w_i
        F[10] = md * e * w_i
        Xw[:, i] = Zw_inv @ F

    # ---------- Hmancal1, Hmancal2, Hdisco ----------
    # Hmancal1: i=1:4, j=1:4 → 0..3
    Hmancal1 = []
    for i in range(0, 4):
        for j in range(0, 4):
            Hmancal1.append(Hw[:, i, j])
    Hmancal1 = np.array(Hmancal1)  # (16, n_omega)

    # Hmancal2: i=21:24, j=21:24 → 20..23
    Hmancal2 = []
    for i in range(20, 24):
        for j in range(20, 24):
            Hmancal2.append(Hw[:, i, j])
    Hmancal2 = np.array(Hmancal2)

    # Hdisco: i=9:12, j=1:4 → 8..11, 0..3
    Hdisco = []
    for i in range(8, 12):
        for j in range(0, 4):
            Hdisco.append(Hw[:, i, j])
    Hdisco = np.array(Hdisco)

    return Wn, qsi, Wd, Hw, Xw, Hmancal1, Hdisco, Hmancal2

# ============================================================
# 6) SIMULAÇÃO EM ROTAÇÃO CONSTANTE – wcons
# ============================================================

def simular_wconst(params, M, G, K_list, C_list, omega, w, idx, t_final=10.0):
    """
    Equivalente ao:

        global md M k e w S
        S = w*G + C{i};
        k = K{i};
        [t,x] = ode15s(@wcons, [0 t], zeros(48,1));

    idx: índice do 'i' (MATLAB) - 1, então 159 → 158 etc.
    """
    ndof = M.shape[0]
    md = params["md"]
    e = params["e"]

    K_i = K_list[idx]
    C_i = C_list[idx]
    S = w * G + C_i

    Minv = np.linalg.inv(M)

    def rhs(t, u):
        x = u[:ndof]
        xdot = u[ndof:]

        F = np.zeros(ndof)
        # F(9) e F(11) → 8 e 10
        F[8] = md * e * w**2 * np.cos(w * t)
        F[10] = md * e * w**2 * np.sin(w * t)

        xddot = Minv @ (F - K_i @ x - S @ xdot)
        return np.concatenate((xdot, xddot))

    y0 = np.zeros(2 * ndof)
    t_eval = np.linspace(0.0, t_final, 10001)

    sol = solve_ivp(rhs, (0.0, t_final), y0, t_eval=t_eval,
                    method="BDF", rtol=1e-6, atol=1e-9)

    # devolve t e estados (igual ao x(t) do MATLAB – só a metade "x")
    x = sol.y[:ndof, :].T  # shape (Nt, 24)
    return sol.t, x


# ============================================================
# 7) MAIN – REPLICA QUASE TUDO DO .M
# ============================================================

def main():
    params = build_parameters()

    # Matrizes FE
    Ke, Me, Ge, G, Ce, M = montar_matrizes_rotor(params)

    # omega (2*pi*0.1:0.1:500)
    Omega1_Hz = np.arange(0.1, 500.0 + 0.1, 0.1)
    omega = 2.0 * np.pi * Omega1_Hz  # rad/s

    # Mancais
    eps, S, S0, alfa, gama, beta, kii, cii = coeficientes_mancal(params, omega)

    # Monta K{i}, C{i}
    K_list, C_list = construir_matrizes_por_velocidade(Ke, Ce, M, G, omega, kii, cii)

    # Campbell + FRFs
    Wn, qsi, Wd, Hw, Xw, Hmancal1, Hdisco, Hmancal2 = calcular_campbell_e_frfs(
        params, M, G, K_list, C_list, omega
    )

    # ----------------- PLOTS ESTÁTICOS (ε, γ, β, k, c) -----------------
    # Excentricidade x Ω
    plt.figure()
    plt.plot(omega / (2*np.pi), eps)
    plt.xlim(0, 50)
    plt.ylim(0, 1)
    plt.xlabel("Velocidade angular Ω (Hz)")
    plt.ylabel("Excentricidade ε")
    plt.title("Excentricidade em função da velocidade angular")
    plt.grid(True)

    # S* x ε e alfa x ε (plotyy simplificado)
    fig, ax1 = plt.subplots()
    ax1.plot(eps, S, "b")
    ax1.set_xlabel("Excentricidade ε")
    ax1.set_ylabel("S*", color="b")
    ax1.set_xlim(0, 1)
    ax1.grid(True)

    ax2 = ax1.twinx()
    ax2.plot(eps, alfa, "r")
    ax2.set_ylabel("Ângulo α (rad)", color="r")
    fig.suptitle("S* e α em função da excentricidade ε")

    # gama_ik
    plt.figure()
    plt.semilogy(eps, gama[:, 0],
                 eps, np.abs(gama[:, 1]),
                 eps, -gama[:, 2],
                 eps, gama[:, 3])
    plt.axis([0, 1, 1e-4, 1e2])
    plt.xlabel("Excentricidade ε")
    plt.ylabel("γ_ik")
    plt.title("Fatores γ_ik em função da Excentricidade ε")
    plt.legend([r"γ_{11}", r"γ_{12}", r"-γ_{21}", r"γ_{22}"])
    plt.grid(True)

    # beta_ik
    plt.figure()
    plt.semilogy(eps, beta[:, 0],
                 eps, -beta[:, 1],
                 eps, -beta[:, 2],
                 eps, beta[:, 3])
    plt.axis([0, 1, 1e0, 1e3])
    plt.xlabel("Excentricidade ε")
    plt.ylabel("β_ik")
    plt.title("Fatores β_ik em função da Excentricidade ε")
    plt.legend([r"β_{11}", r"-β_{12}", r"-β_{21}", r"β_{22}"])
    plt.grid(True)

    # kii x Ω
    plt.figure()
    plt.plot(omega/(2*np.pi), kii[:, 0],
             omega/(2*np.pi), kii[:, 1],
             omega/(2*np.pi), kii[:, 2],
             omega/(2*np.pi), kii[:, 3])
    plt.axis([0, 50, -1e6, 1e6])
    plt.xlabel("Velocidade angular Ω (Hz)")
    plt.ylabel("k_ik")
    plt.title("Fatores k_ik em função da Velocidade angular Ω")
    plt.legend([r"k_{11}", r"k_{12}", r"k_{21}", r"k_{22}"])
    plt.grid(True)

    # cii x Ω
    plt.figure()
    plt.plot(omega/(2*np.pi), cii[:, 0],
             omega/(2*np.pi), cii[:, 1],
             omega/(2*np.pi), cii[:, 2],
             omega/(2*np.pi), cii[:, 3])
    plt.axis([0, 50, -5e5, 1.5e6])
    plt.xlabel("Velocidade angular Ω (Hz)")
    plt.ylabel("c_ik")
    plt.title("Fatores c_ik em função da Velocidade angular Ω")
    plt.legend([r"c_{11}", r"c_{12}", r"c_{21}", r"c_{22}"])
    plt.grid(True)

    # ----------------- DIAGRAMA DE CAMPBELL -----------------
    plt.figure()
    plt.plot(omega/(2*np.pi), omega/(2*np.pi), "--", linewidth=2)
    plt.plot(omega/(2*np.pi), omega/np.pi, "--", linewidth=2)
    plt.plot(omega/(2*np.pi), Wn[:, :24] / (2*np.pi))
    plt.axis([0, 500, 0, 1000])
    plt.title("Diagrama de Campbell")
    plt.xlabel("Velocidade angular Ω (Hz)")
    plt.ylabel("Frequência ω (Hz)")
    plt.legend([r"ω = Ω", r"ω = 2Ω"])
    plt.grid(True)

    # ----------------- FRFs (mancais e disco) -----------------
    so = [0, 2, 8, 10]  # DOFs 1,3,9,11 (MATLAB) -> 0,2,8,10

    # Nó 1 – Mancal 1
    plt.figure()
    plt.subplot(3, 2, 1)
    plt.plot(omega/(2*np.pi), mag2db(Hmancal1[so, :].T))
    plt.xlabel("Frequência (Hz)")
    plt.ylabel("H amplitude (dB)")
    plt.title("Nó 1 - Mancal 1 FRF amplitude")

    plt.subplot(3, 2, 2)
    plt.plot(omega/(2*np.pi), 180/np.pi * np.angle(Hmancal1[so, :].T))
    plt.xlabel("Frequência (Hz)")
    plt.ylabel("H fase (graus)")
    plt.title("Nó 1 - Mancal 1 FRF fase")

    # Nó 6 – Mancal 2
    plt.subplot(3, 2, 3)
    plt.plot(omega/(2*np.pi), mag2db(Hmancal2[so, :].T))
    plt.xlabel("Frequência (Hz)")
    plt.ylabel("H amplitude (dB)")
    plt.title("Nó 6 - Mancal 2 FRF amplitude")

    plt.subplot(3, 2, 4)
    plt.plot(omega/(2*np.pi), 180/np.pi * np.angle(Hmancal2[so, :].T))
    plt.xlabel("Frequência (Hz)")
    plt.ylabel("H fase (graus)")
    plt.title("Nó 6 - Mancal 2 FRF fase")

    # Nó 3 – Disco
    plt.subplot(3, 2, 5)
    plt.plot(omega/(2*np.pi), mag2db(Hdisco[so, :].T))
    plt.xlabel("Frequência (Hz)")
    plt.ylabel("H amplitude (dB)")
    plt.title("Nó 3 - Disco FRF amplitude")

    plt.subplot(3, 2, 6)
    plt.plot(omega/(2*np.pi), 180/np.pi * np.angle(Hdisco[so, :].T))
    plt.xlabel("Frequência (Hz)")
    plt.ylabel("H fase (graus)")
    plt.title("Nó 3 - Disco FRF fase")

    # ----------------- RESPOSTAS NO TEMPO E ÓRBITAS -----------------
    # Mesmos casos do MATLAB: w e i "na mão"
    t_final = 10.0

    # Caso 2 primeiro (w=15.9*2π, i=159)
    w2 = 15.9 * 2.0 * np.pi
    idx2 = 159 - 1
    t2, x2 = simular_wconst(params, M, G, K_list, C_list, omega, w2, idx2, t_final)

    # Caso 1 (w=10, i=100)
    w1 = 10.0
    idx1 = 100 - 1
    t1, x1 = simular_wconst(params, M, G, K_list, C_list, omega, w1, idx1, t_final)

    # w=50*2π, i=500
    w3 = 50.0 * 2.0 * np.pi
    idx3 = 500 - 1
    t3, x3 = simular_wconst(params, M, G, K_list, C_list, omega, w3, idx3, t_final)

    # w=100*2π, i=1000
    w4 = 100.0 * 2.0 * np.pi
    idx4 = 1000 - 1
    t4, x4 = simular_wconst(params, M, G, K_list, C_list, omega, w4, idx4, t_final)

    # w=109.3*2π, i=1093
    w5 = 109.3 * 2.0 * np.pi
    idx5 = 1093 - 1
    t5, x5 = simular_wconst(params, M, G, K_list, C_list, omega, w5, idx5, t_final)

    # w=180*2π, i=1800
    w6 = 180.0 * 2.0 * np.pi
    idx6 = 1800 - 1
    t6, x6 = simular_wconst(params, M, G, K_list, C_list, omega, w6, idx6, t_final)

    # Exemplo de órbita (tipo plots com rplot e theta)
    # Vou só mostrar para o caso w5 (109.3 Hz) como no teu "3"
    num = -1  # último passo de tempo
    rplot = np.sqrt(x5[num, 0:21:4]**2 + x5[num, 2:23:4]**2)  # nós 1..6

    theta = np.linspace(0.0, 2*np.pi, 100)
    y = np.zeros((6, theta.size))
    z = np.zeros_like(y)
    for k in range(6):
        y[k, :] = rplot[k] * np.cos(theta)
        z[k, :] = rplot[k] * np.sin(theta)

    # 3D plot da forma (igual teu plot3)
    xs = np.linspace(0.0, 0.8, 6)
    plt.figure()
    for k in range(6):
        plt.plot(xs[k] * np.ones_like(theta), y[k, :], z[k, :], "-k", linewidth=2)
    plt.plot(xs, np.zeros(6), np.zeros(6), "-b", linewidth=3)
    plt.xlabel("Distância desde o mancal 1 x (m)")
    plt.ylabel("Deslocamento y (m)")
    plt.zlabel = "Deslocamento z (m)"  # matplotlib não tem zlabel diretamente em 2D; se quiser, transforma em 3D Axes3D
    plt.grid(True)

    # Plots de respostas no tempo (só um resumo, dá pra detalhar como no .m)
    # Ex.: w=10 (t1,x1) – nós 1, 3, 6 em y e z
    plt.figure()
    plt.subplot(3, 2, 1)
    plt.plot(t1, 1000 * x1[:, 0])
    plt.xlabel("Tempo (s)")
    plt.ylabel("Nó 1, y (mm)")
    plt.title("Nó 1 - Mancal 1 w=10 rad/s")

    plt.subplot(3, 2, 2)
    plt.plot(t1, 1000 * x1[:, 2])
    plt.xlabel("Tempo (s)")
    plt.ylabel("Nó 1, z (mm)")

    plt.subplot(3, 2, 3)
    plt.plot(t1, 1000 * x1[:, 8])
    plt.xlabel("Tempo (s)")
    plt.ylabel("Nó 3, y (mm)")

    plt.subplot(3, 2, 4)
    plt.plot(t1, 1000 * x1[:, 10])
    plt.xlabel("Tempo (s)")
    plt.ylabel("Nó 3, z (mm)")

    plt.subplot(3, 2, 5)
    plt.plot(t1, 1000 * x1[:, 20])
    plt.xlabel("Tempo (s)")
    plt.ylabel("Nó 6, y (mm)")

    plt.subplot(3, 2, 6)
    plt.plot(t1, 1000 * x1[:, 22])
    plt.xlabel("Tempo (s)")
    plt.ylabel("Nó 6, z (mm)")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
