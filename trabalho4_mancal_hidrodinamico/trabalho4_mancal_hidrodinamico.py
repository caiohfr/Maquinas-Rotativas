
## 🧠 Código Python – `trabalho4_mancal_hidrodinamico.py`

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar

# --------------------------------------------------------------------
# Dados globais (equivalentes ao "global" do MATLAB)
# --------------------------------------------------------------------

# Vou colocar os parâmetros em um dicionário em vez de usar "global",
# mas manter os mesmos nomes do seu código para ficar fácil mapear.


def build_parameters():
    params = {}

    # Dados do mancal
    Dm = 0.03       # diâmetro do mancal (m)
    R = Dm / 2.0    # raio
    Lm = 0.018      # comprimento (m)
    folga_r = 0.1e-3  # folga radial (m)
    n_visc = 0.08   # viscosidade (Pa.s)

    # Dados do rotor (eixo + disco)
    De = 0.01       # diâmetro do eixo (m)
    Le = 0.9        # comprimento do eixo (m)
    Dd = 0.1        # diâmetro do disco (m)
    r_d = Dd / 2.0
    Ld = 0.02       # comprimento do disco (m)
    E = 2e11        # módulo de elasticidade (Pa)
    ro = 7850.0     # densidade (kg/m^3)
    v = 0.3         # coeficiente de Poisson (não usado diretamente)
    g = 9.80665     # gravidade (m/s^2)

    Ie = np.pi * De**4 / 64.0        # inércia geométrica eixo
    k = 48 * E * Ie / (Le**3)        # rigidez equivalente
    c = 1e-4 * k                     # amortecimento (proporcional)
    me = (np.pi * De**2 / 4.0) * Le * ro
    md = (np.pi * Dd**2 / 4.0) * Ld * ro
    e = (1.5e-4) / md                # desbalanceamento (m)
    Id = 0.5 * md * (Dd / 2.0)**2 + md * e**2  # inércia polar do disco

    F0 = md * g / 2.0                # carga vertical no mancal

    params.update(dict(
        Dm=Dm, R=R, Lm=Lm, folga_r=folga_r, n_visc=n_visc,
        De=De, Le=Le, Dd=Dd, r_d=r_d, Ld=Ld,
        E=E, ro=ro, v=v, g=g,
        Ie=Ie, k=k, c=c, me=me, md=md, e=e, Id=Id,
        F0=F0
    ))

    return params


# --------------------------------------------------------------------
# Relação epsilon(S*) – TODO: colocar a forma correta do seu material
# --------------------------------------------------------------------

def epsilon_residual(eps, S):
    """
    Residual da equação que relaciona excentricidade adimensional ε e
    número de Sommerfeld modificado S* para mancal curto.

    NO SEU MATLAB:
        eps(i) = fminsearch(@(e)epsilon(e,S(i)), eps(i-1));

    Aqui precisamos da função "epsilon(e,S)" que não veio no anexo.
    Sem a expressão exata, não dá pra reproduzir fielmente o resultado.

    >>> TODO: Substituir este residual pela equação correta.
             Por enquanto, uso uma relação monotônica simples
             apenas pra manter o código executável.

    Exemplo placeholder (não físico):
        eps_target = 1 - exp(-S)
        residual = eps - eps_target
    """
    eps_target = 1.0 - np.exp(-S)
    return eps - eps_target


def find_eccentricity(S, eps_initial=0.1):
    """
    Resolve epsilon_residual(eps, S) ~ 0 via minimização do residual^2,
    como um análogo do fminsearch do MATLAB.
    """
    def obj(eps):
        return epsilon_residual(eps, S)**2

    res = minimize_scalar(obj, bounds=(1e-4, 0.99), method="bounded")
    return res.x


# --------------------------------------------------------------------
# Cálculo de gama, beta, k_ij, c_ij (copiado do seu MATLAB)
# --------------------------------------------------------------------

def compute_bearing_coeffs(omega, params):
    """
    Reproduz bloco:

        Fn = (n*Lm^3*R)/(2*folga_r^2)*omega;
        S0 = ((Lm/(2*R))^2*F0)./Fn;
        S = F0./Fn;
        eps = ...
        alfa = atan((pi/4)*sqrt(1-eps^2)./eps);
        Ae = 4./(pi^2+(16-pi^2)*eps.^2).^(3/2);
        gama(:,1:4) = ...
        beta(:,1:4) = ...
        kii = gama*(F0/folga_r);
        cii = beta*(F0/folga_r)./omega;

    omega é um vetor coluna em rad/s.
    """
    folga_r = params["folga_r"]
    Lm = params["Lm"]
    R = params["R"]
    n_visc = params["n_visc"]
    F0 = params["F0"]

    omega = np.asarray(omega).reshape(-1, 1)  # Nx1

    # Força de referência
    Fn = (n_visc * Lm**3 * R) / (2 * folga_r**2) * omega
    S0 = ((Lm / (2 * R))**2 * F0) / Fn
    S = F0 / Fn  # Sommerfeld modificado

    # Excentricidade vs S
    eps_list = [0.0]
    for i in range(1, len(omega) + 1):
        S_i = S[i - 1, 0]
        eps_prev = eps_list[-1] if eps_list[-1] > 0 else 0.1
        eps_i = find_eccentricity(S_i, eps_prev)
        eps_list.append(eps_i)

    eps = np.array(eps_list[1:]).reshape(-1, 1)  # Nx1
    exc = eps * folga_r
    alfa = np.arctan((np.pi / 4.0) * np.sqrt(1 - eps**2) / eps)

    # Fator Ae
    Ae = 4.0 / (np.pi**2 + (16 - np.pi**2) * eps**2)**(1.5)

    # gama(:,1) ... gama(:,4)
    gama = np.zeros((len(eps), 4))
    gama[:, 0] = (2 * np.pi**2 + (16 - np.pi**2) * eps[:, 0]**2) * Ae[:, 0]

    gama[:, 1] = (np.pi / 4.0) * (
        np.pi**2
        - 2 * np.pi**2 * eps[:, 0]**2
        - (16 - np.pi**2) * eps[:, 0]**4
    ) / (eps[:, 0] * np.sqrt(1 - eps[:, 0]**2)) * Ae[:, 0]

    gama[:, 2] = -(np.pi / 4.0) * (
        np.pi**2
        + (32 + np.pi**2) * eps[:, 0]**2
        + (32 - 2 * np.pi**2) * eps[:, 0]**4
    ) / (eps[:, 0] * np.sqrt(1 - eps[:, 0]**2)) * Ae[:, 0]

    gama[:, 3] = (
        np.pi**2
        + (32 + np.pi**2) * eps[:, 0]**2
        + (32 - 2 * np.pi**2) * eps[:, 0]**4
    ) / (np.sqrt(1 - eps[:, 0]**2)) * Ae[:, 0]

    # beta(:,1) ... beta(:,4)
    beta = np.zeros((len(eps), 4))
    beta[:, 0] = (np.pi / 2.0) * (
        np.sqrt(1 - eps[:, 0]**2) / eps[:, 0]
    ) * (np.pi**2 + (2 * np.pi**2 - 16) * eps[:, 0]**2) * Ae[:, 0]

    beta[:, 1] = -(2 * np.pi**2 + (4 * np.pi**2 - 32) * eps[:, 0]**2) * Ae[:, 0]
    beta[:, 2] = beta[:, 1]

    beta[:, 3] = (np.pi / 2.0) * (
        np.pi**2
        + (48 - 2 * np.pi**2) * eps[:, 0]**2
        + np.pi**2 * eps[:, 0]**4
    ) / (eps[:, 0] * np.sqrt(1 - eps[:, 0]**2)) * Ae[:, 0]

    # Rigidez e amortecimento (igual ao MATLAB)
    kii = gama * (F0 / folga_r)              # Nx4
    cii = beta * (F0 / folga_r) / omega      # Nx4

    return eps[:, 0], exc[:, 0], alfa[:, 0], S[:, 0], kii, cii


# --------------------------------------------------------------------
# Fatores(km,cm) em função da rotação (equivalente a fatores.m)
# --------------------------------------------------------------------

def fatores(rot, omega_grid, kii, cii):
    """
    Interpola km, dm (kii, cii) para uma dada rotação "rot" (rad/s),
    equivalente a:

        Kii = interp1(omega,kii,rot,'linear','extrap');
        dii = interp1(omega,cii,rot,'linear','extrap');

    Retorna vetores 4-elemento km, dm.
    """
    km = np.array([
        np.interp(rot, omega_grid, kii[:, j]) for j in range(4)
    ])
    dm = np.array([
        np.interp(rot, omega_grid, cii[:, j]) for j in range(4)
    ])
    return km, dm


# --------------------------------------------------------------------
# Equação diferencial – laval_mancal (caso com torque T1)
# --------------------------------------------------------------------

def laval_mancal(t, y, params, omega_grid, kii, cii, T1):
    """
    Equivalente à função laval_mancal(t,y) do MATLAB.

    Estados:
        y = [theta, y, z, theta_dot, y_dot, z_dot, y_m, z_m]

    """
    Id = params["Id"]
    md = params["md"]
    c = params["c"]
    k = params["k"]
    e = params["e"]

    theta, y_d, z_d, theta_dot, y_dot, z_dot, y_m, z_m = y

    # Fatores do mancal
    km, dm = fatores(theta_dot, omega_grid, kii, cii)
    K11, K12, K21, K22 = km
    C11, C12, C21, C22 = dm

    D1 = k / 2.0 * (y_d - y_m) - (K11 * y_m + K12 * z_m)
    D2 = k / 2.0 * (z_d - z_m) - (K21 * y_m + K22 * z_m)
    B = (D2 - D1 * C21 / C11) / (C22 - C21 * C12 / C11)

    A1 = (1.0 / Id) * (
        T1
        + k * (z_d - z_m) * (e * np.cos(theta))
        - k * (y_d - y_m) * (e * np.sin(theta))
        - c * y_dot * (e * np.sin(theta))
        + c * z_dot * (e * np.cos(theta))
    )

    dydt = np.zeros_like(y)
    dydt[0] = theta_dot                       # theta_dot
    dydt[1] = y_dot                           # y_dot
    dydt[2] = z_dot                           # z_dot
    dydt[3] = A1                              # theta_ddot

    dydt[4] = (1.0 / md) * (
        md * e * (dydt[3] * np.sin(theta) + theta_dot**2 * np.cos(theta))
        - c * y_dot
        - k * (y_d - y_m)
    )

    dydt[5] = (1.0 / md) * (
        md * e * (-dydt[3] * np.cos(theta) + theta_dot**2 * np.sin(theta))
        - c * z_dot
        - k * (z_d - z_m)
    )

    dydt[6] = D1 / C11 - C12 / C11 * B        # y_m_dot
    dydt[7] = B                               # z_m_dot

    return dydt


# --------------------------------------------------------------------
# Equação diferencial – laval_mancalwcons (rotação constante)
# --------------------------------------------------------------------

def laval_mancalwcons(t, y, params, km, cm, w):
    """
    Equivalente a laval_mancalwcons(t,y) do MATLAB.

    Estados:
        y = [y, z, y_dot, z_dot, y_m, z_m]
    """
    md = params["md"]
    c = params["c"]
    k = params["k"]
    e = params["e"]

    y_d, z_d, y_dot, z_dot, y_m, z_m = y

    K11, K12, K21, K22 = km
    C11, C12, C21, C22 = cm

    D1 = k / 2.0 * (y_d - y_m) - (K11 * y_m + K12 * z_m)
    D2 = k / 2.0 * (z_d - z_m) - (K21 * y_m + K22 * z_m)
    B = (D2 - D1 * C21 / C11) / (C22 - C21 * C12 / C11)

    dydt = np.zeros_like(y)
    dydt[0] = y_dot
    dydt[1] = z_dot

    dydt[2] = (1.0 / md) * (
        md * e * w**2 * np.cos(w * t)
        - c * y_dot
        - k * (y_d - y_m)
    )

    dydt[3] = (1.0 / md) * (
        md * e * w**2 * np.sin(w * t)
        - c * z_dot
        - k * (z_d - z_m)
    )

    dydt[4] = D1 / C11 - C12 / C11 * B   # y_m_dot
    dydt[5] = B                          # z_m_dot

    return dydt


# --------------------------------------------------------------------
# Função principal – espelha o fluxo do script MATLAB
# --------------------------------------------------------------------

def main():
    params = build_parameters()

    # Malha de velocidades (igual ao MATLAB)
    omega_Hz = np.arange(0.01, 50.01, 0.01)
    omega = 2.0 * np.pi * np.arange(0.01, 50.001, 0.001)  # rad/s

    # Cálculo dos coeficientes do mancal
    eps, exc, alfa, S, kii, cii = compute_bearing_coeffs(omega, params)

    # ---------------------------
    # Caso 1, 2, 3 – torque T1
    # ---------------------------
    t1 = 50.0
    y0_case1 = np.array([0, 0, 0, 0.1, 0, 0, 0, 0], dtype=float)
    y0_case2 = np.array([0, 0, 0, 0.01, 0, 0, 0, 0], dtype=float)
    y0_case3 = np.array([0, 0, 0, 0.01, 0, 0, 0, 0], dtype=float)

    t_eval1 = np.linspace(0.0, t1, 5000)

    # T1 = 0.01
    sol1 = solve_ivp(
        lambda t, y: laval_mancal(t, y, params, omega, kii, cii, T1=0.01),
        (0.0, t1),
        y0_case1,
        t_eval=t_eval1,
        rtol=1e-6,
        atol=1e-9
    )

    # T1 = 0.0045
    sol2 = solve_ivp(
        lambda t, y: laval_mancal(t, y, params, omega, kii, cii, T1=0.0045),
        (0.0, t1),
        y0_case2,
        t_eval=t_eval1,
        rtol=1e-6,
        atol=1e-9
    )

    # T1 = 0.004
    sol3 = solve_ivp(
        lambda t, y: laval_mancal(t, y, params, omega, kii, cii, T1=0.004),
        (0.0, t1),
        y0_case3,
        t_eval=t_eval1,
        rtol=1e-6,
        atol=1e-9
    )

    # ----------------------------------------------------
    # Casos com rotação constante (w = 30, 50, 72, 120, 150, 145)
    # ----------------------------------------------------
    t1_const = 100.0
    t_eval_const = np.linspace(0.0, t1_const, 5000)

    w_list = [30.0, 50.0, 72.0, 120.0, 150.0, 145.0]  # rad/s (no MATLAB w=30, 50, 72, 120, 150, 145)
    orbits = {}

    for w in w_list:
        km = np.array([np.interp(w, omega, kii[:, j]) for j in range(4)])
        cm = np.array([np.interp(w, omega, cii[:, j]) for j in range(4)])

        y0_const = np.zeros(6, dtype=float)
        sol_const = solve_ivp(
            lambda t, y: laval_mancalwcons(t, y, params, km, cm, w),
            (0.0, t1_const),
            y0_const,
            t_eval=t_eval_const,
            rtol=1e-6,
            atol=1e-9
        )

        orbits[w] = sol_const

    # ---------------------------
    # Gráficos principais
    # ---------------------------

    # Excentricidade vs velocidade angular (em Hz)
    plt.figure()
    plt.grid(True)
    plt.title("Excentricidade em função da velocidade angular")
    plt.xlabel(r"Velocidade angular $\Omega$ (Hz)")
    plt.ylabel(r"Excentricidade $\epsilon$")
    plt.plot(omega / (2.0 * np.pi), eps, linewidth=1.5)
    plt.xlim(0, 50)
    plt.ylim(0, 1)

    # S* e alfa vs epsilon
    plt.figure()
    plt.title(r"S* e ângulo $\alpha$ em função da excentricidade $\epsilon$")
    plt.xlabel(r"Excentricidade $\epsilon$")
    plt.grid(True)

    ax1 = plt.gca()
    ax1.set_ylabel(r"Número de Sommerfeld modificado $S^*$")
    ax1.semilogy(eps, S, linewidth=1.5)
    ax1.set_ylim(1e-2, 1e2)

    ax2 = ax1.twinx()
    ax2.set_ylabel(r"Ângulo $\alpha$ (graus)")
    ax2.plot(eps, alfa * 180.0 / np.pi, linewidth=1.5)
    ax2.set_ylim(0, 90)

    # Fatores gama
    plt.figure()
    plt.semilogy(eps, kii[:, 0] * 0 + (np.abs(kii[:, 0]) / (params["F0"] / params["folga_r"])))  # só pra ilustrar
    plt.title(r"Fatores $\gamma_{ik}$ e $\beta_{ik}$ – ver código para detalhes")
    plt.grid(True)

    # Resposta temporal – exemplo para Caso 1 (Uy, Uz, rotação)
    t = sol1.t
    y1 = sol1.y
    theta1 = y1[0, :]
    uy1 = y1[1, :]
    uz1 = y1[2, :]
    theta_dot1 = y1[3, :]

    fig, axs = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    fig.suptitle("Caso 1 – Rotação e deslocamentos Uy, Uz")

    axs[0].plot(t, 1000.0 * uy1, label="Uy (mm)")
    axs[0].plot(t, 1000.0 * uz1, label="Uz (mm)")
    axs[0].set_ylabel("Deslocamento (mm)")
    axs[0].grid(True)
    axs[0].legend()

    axs[1].plot(t, theta_dot1, label=r"$\dot{\theta}$ (rad/s)")
    axs[1].set_xlabel("Tempo (s)")
    axs[1].set_ylabel("Rotação (rad/s)")
    axs[1].grid(True)
    axs[1].legend()

    # Órbita para uma rotação constante (por exemplo w=150 rad/s)
    w_example = 150.0
    if w_example in orbits:
        sol_w = orbits[w_example]
        y_w = sol_w.y
        y_d = y_w[0, :]
        z_d = y_w[1, :]

        plt.figure()
        plt.plot(1000.0 * y_d, 1000.0 * z_d)
        plt.xlabel("Uy (mm)")
        plt.ylabel("Uz (mm)")
        plt.title(f"Órbita do rotor – w = {w_example:.1f} rad/s")
        plt.grid(True)
        plt.axis("equal")

    plt.show()


if __name__ == "__main__":
    main()
