import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
from scipy.integrate import solve_ivp


# -------------------------------------------------------------
# Parâmetros do mancal e rotor  (mesmos do MATLAB)
# -------------------------------------------------------------
def build_parameters():
    Dm = 0.03       # diâmetro do mancal (m)
    R = Dm / 2.0    # raio do mancal (m)
    Lm = 0.018      # comprimento do mancal (m)
    folga_r = 0.1e-3  # folga radial (m)
    n_visc = 0.08   # viscosidade absoluta (Pa.s)

    De = 0.01       # diâmetro do eixo (m)
    Le = 0.9        # comprimento do eixo (m)
    Dd = 0.10       # diâmetro do disco (m)
    Ld = 0.02       # comprimento do disco (m)

    E = 2.0e11      # módulo de elasticidade (Pa)
    ro = 7850.0     # densidade (kg/m³)
    v = 0.3         # Poisson (não entra aqui)
    g = 9.80665     # gravidade (m/s²)

    Ie = np.pi * De**4 / 64.0         # momento de inércia de área do eixo
    k = 48.0 * E * Ie / (Le**3)       # rigidez do eixo
    c = 1e-4 * k                      # amortecimento viscoso do eixo

    me = (np.pi * De**2 / 4.0) * Le * ro
    md = (np.pi * Dd**2 / 4.0) * Ld * ro
    e = 1.5e-4 / md                   # desbalanceamento (m)

    Id = 0.5 * md * (Dd / 2.0) ** 2 + md * e**2  # inércia polar do disco

    F0 = md * g / 2.0                 # carga no mancal (axial -> reação)

    return {
        "Dm": Dm, "R": R, "Lm": Lm, "folga_r": folga_r, "n_visc": n_visc,
        "De": De, "Le": Le, "Dd": Dd, "Ld": Ld,
        "E": E, "ro": ro, "v": v, "g": g,
        "Ie": Ie, "k": k, "c": c, "me": me, "md": md, "e": e, "Id": Id,
        "F0": F0
    }


# -------------------------------------------------------------
# Sommerfeld modificado S*(ε) – exatamente eq. (15) / epsilon.m
# S*(ε) = (π/2) * ε / (1-ε²)² * sqrt(1-ε² + (4ε/π)²)
# -------------------------------------------------------------
def S_star_from_eps(eps):
    inside = 1.0 - eps**2 + (4.0 * eps / np.pi) ** 2
    return (np.pi / 2.0) * (eps / (1.0 - eps**2) ** 2) * np.sqrt(inside)


def eps_from_S_star(S_target, tol=1e-8):
    """
    Resolve S*(ε) = S_target para 0<ε<1 via busca em intervalo.
    S*(ε) é crescente em (0,1), então dá pra usar busca binária
    embutida no minimize_scalar com bounds.
    """
    if S_target <= 0.0:
        return 0.0

    def obj(e):
        return abs(S_star_from_eps(e) - S_target)

    res = minimize_scalar(obj, bounds=(1e-6, 0.999999), method="bounded")
    return float(res.x)


# -------------------------------------------------------------
# Cálculo de ε(Ω), S*, γik, βik, kik, cik
# -------------------------------------------------------------
def compute_bearing_coeffs(params):
    R = params["R"]
    Lm = params["Lm"]
    folga_r = params["folga_r"]
    n_visc = params["n_visc"]
    F0 = params["F0"]

    # mesmo range de frequência da parte de torque no MATLAB (em Hz)
    omega_Hz = np.arange(0.01, 50.0 + 0.01, 0.01)  # 0.01:0.01:50
    omega = 2.0 * np.pi * omega_Hz                 # rad/s

    # para os coeficientes o MATLAB usa um passo mais fino (0.001),
    # mas para Python isso aqui já é bem próximo e bem mais leve.
    # Se quiser idêntico, troque para 0.001 aqui e em omega_Hz.

    # Fη (eq. 15): Fη = η L³ R / (2 δ²) * Ω
    F_eta = (n_visc * Lm**3 * R) / (2.0 * folga_r**2) * omega_Hz

    # Sommerfelds
    S0 = ((Lm / (2.0 * R)) ** 2 * F0) / F_eta   # clássico (não usado depois)
    S_star = F0 / F_eta                         # modificado (eq. 15)

    # ε(S*)
    eps = np.zeros_like(omega)
    for i in range(omega.size):
        eps[i] = eps_from_S_star(S_star[i])

    # Excentricidade dimensional (se quiser comparar com exc do MATLAB)
    exc = eps * folga_r

    # α(ε) – eq. (16)
    alfa = np.arctan((np.pi / 4.0) * (np.sqrt(1.0 - eps**2) / eps))

    # A(ε) – mesma forma do MATLAB
    Ae = 4.0 / (np.pi**2 + (16.0 - np.pi**2) * eps**2) ** 1.5

    # γik
    gama11 = (2.0 * np.pi**2 + (16.0 - np.pi**2) * eps**2) * Ae
    gama12 = (
        (np.pi / 4.0)
        * (np.pi**2 - 2.0 * np.pi**2 * eps**2 - (16.0 - np.pi**2) * eps**4)
        / (eps * np.sqrt(1.0 - eps**2))
        * Ae
    )
    gama21 = -(
        (np.pi / 4.0)
        * (np.pi**2 + (32.0 + np.pi**2) * eps**2
           + (32.0 - 2.0 * np.pi**2) * eps**4)
        / (eps * np.sqrt(1.0 - eps**2))
        * Ae
    )
    gama22 = (
        (np.pi**2 + (32.0 + np.pi**2) * eps**2
         + (32.0 - 2.0 * np.pi**2) * eps**4)
        / np.sqrt(1.0 - eps**2)
        * Ae
    )
    gama = np.column_stack((gama11, gama12, gama21, gama22))

    # βik
    beta11 = (
        (np.pi / 2.0)
        * (np.sqrt(1.0 - eps**2) / eps)
        * (np.pi**2 + (2.0 * np.pi**2 - 16.0) * eps**2)
        * Ae
    )
    beta12 = -(2.0 * np.pi**2 + (4.0 * np.pi**2 - 32.0) * eps**2) * Ae
    beta21 = beta12.copy()
    beta22 = (
        (np.pi / 2.0)
        * (np.pi**2 + (48.0 - 2.0 * np.pi**2) * eps**2 + np.pi**2 * eps**4)
        / (eps * np.sqrt(1.0 - eps**2))
        * Ae
    )
    beta = np.column_stack((beta11, beta12, beta21, beta22))

    # Rigidezes e amortecimentos (eq. 17)
    kii = gama * (F0 / folga_r)                    # N/m
    cii = beta * (F0 / folga_r) / omega[:, None]   # N·s/m

    return omega_Hz, omega, eps, S_star, alfa, gama, beta, kii, cii


# -------------------------------------------------------------
# Interpolação dos fatores de mancal (equivalente a fatores.m)
# -------------------------------------------------------------
def bearing_factors(rot, omega, kii, cii):
    """
    rot: velocidade angular instantânea (rad/s)
    omega: vetor (rad/s) usado para tabular kii/cii
    retorna K11..K22, C11..C22
    """
    K = np.array([np.interp(rot, omega, kii[:, j]) for j in range(4)])
    C = np.array([np.interp(rot, omega, cii[:, j]) for j in range(4)])
    return K, C


# -------------------------------------------------------------
# ODE – laval_mancal (caso com torque T1, igual ao MATLAB)
# -------------------------------------------------------------
def laval_mancal_rhs(t, y, params, omega, kii, cii, T1):
    Id = params["Id"]
    md = params["md"]
    c = params["c"]
    k = params["k"]
    e = params["e"]

    # estados: y = [theta, uy, uz, theta_dot, uy_dot, uz_dot, ym, zm]
    theta, uy, uz, theta_dot, uy_dot, uz_dot, ym, zm = y

    # coeficientes do mancal para a rotação atual
    km, dm = bearing_factors(theta_dot, omega, kii, cii)
    K11, K12, K21, K22 = km
    C11, C12, C21, C22 = dm

    # D1, D2 e B (mesma álgebra do MATLAB)
    D1 = k / 2.0 * (uy - ym) - (K11 * ym + K12 * zm)
    D2 = k / 2.0 * (uz - zm) - (K21 * ym + K22 * zm)
    B = (D2 - D1 * C21 / C11) / (C22 - C21 * C12 / C11)

    # equação angular
    A1 = (1.0 / Id) * (
        T1
        + k * (uz - zm) * (e * np.cos(theta))
        - k * (uy - ym) * (e * np.sin(theta))
        - c * uy_dot * (e * np.sin(theta))
        + c * uz_dot * (e * np.cos(theta))
    )

    dydt = np.zeros_like(y)
    dydt[0] = theta_dot
    dydt[1] = uy_dot
    dydt[2] = uz_dot
    dydt[3] = A1
    dydt[4] = (
        1.0 / md * (
            md * e * (dydt[3] * np.sin(theta) + theta_dot**2 * np.cos(theta))
            - c * uy_dot
            - k * (uy - ym)
        )
    )
    dydt[5] = (
        1.0 / md * (
            md * e * (-dydt[3] * np.cos(theta) + theta_dot**2 * np.sin(theta))
            - c * uz_dot
            - k * (uz - zm)
        )
    )
    dydt[6] = D1 / C11 - C12 / C11 * B
    dydt[7] = B
    return dydt


# -------------------------------------------------------------
# ODE – laval_mancalwcons (Ω constante)
# -------------------------------------------------------------
def laval_mancal_const_rhs(t, y, params, w, km, cm):
    md = params["md"]
    c = params["c"]
    k = params["k"]
    e = params["e"]

    # estados: y = [uy, uz, uy_dot, uz_dot, ym, zm]
    uy, uz, uy_dot, uz_dot, ym, zm = y

    K11, K12, K21, K22 = km
    C11, C12, C21, C22 = cm

    D1 = k / 2.0 * (uy - ym) - (K11 * ym + K12 * zm)
    D2 = k / 2.0 * (uz - zm) - (K21 * ym + K22 * zm)
    B = (D2 - D1 * C21 / C11) / (C22 - C21 * C12 / C11)

    dydt = np.zeros_like(y)
    dydt[0] = uy_dot
    dydt[1] = uz_dot
    dydt[2] = (1.0 / md) * (
        md * e * w**2 * np.cos(w * t) - c * uy_dot - k * (uy - ym)
    )
    dydt[3] = (1.0 / md) * (
        md * e * w**2 * np.sin(w * t) - c * uz_dot - k * (uz - zm)
    )
    dydt[4] = D1 / C11 - C12 / C11 * B
    dydt[5] = B
    return dydt


# -------------------------------------------------------------
# Plots estáticos (ε, S*, α, γik, βik, kik, cik)
# -------------------------------------------------------------
def plot_static_results(omega_Hz, omega, eps, S_star, alfa, gama, beta, kii, cii):
    # ε x Ω
    plt.figure()
    plt.plot(omega_Hz, eps, linewidth=1.5)
    plt.grid(True)
    plt.xlim(0.0, 50.0)
    plt.ylim(0.0, 1.0)
    plt.xlabel("Velocidade angular Ω (Hz)")
    plt.ylabel("Excentricidade ε")
    plt.title("Excentricidade em função da velocidade angular")

    # S* e α x ε
    fig, ax1 = plt.subplots()
    ax1.semilogy(eps, S_star, "b", linewidth=1.5, label="S*")
    ax1.set_xlabel("Excentricidade ε")
    ax1.set_ylabel("S*", color="b")
    ax1.grid(True)
    ax1.set_xlim(0.0, 1.0)
    ax1.set_ylim(1e-2, 1e2)

    ax2 = ax1.twinx()
    ax2.plot(eps, np.degrees(alfa), "r", linewidth=1.5, label="α")
    ax2.set_ylabel("Ângulo α (graus)", color="r")
    ax2.set_ylim(0.0, 90.0)

    fig.suptitle("S* e ângulo α em função da excentricidade ε")

    # γik x ε
    plt.figure()
    labels_g = ["γ11", "γ12", "-γ21", "γ22"]
    series_g = [gama[:, 0], gama[:, 1], -gama[:, 2], gama[:, 3]]
    for s, lab in zip(series_g, labels_g):
        plt.semilogy(eps, s, linewidth=1.5, label=lab)
    plt.grid(True)
    plt.xlim(0.0, 1.0)
    plt.ylim(10**-4, 10**2)
    plt.xlabel("Excentricidade ε")
    plt.ylabel("γik")
    plt.title("Fatores γik em função da excentricidade ε")
    plt.legend()

    # βik x ε
    plt.figure()
    labels_b = ["β11", "-β12", "-β21", "β22"]
    series_b = [beta[:, 0], -beta[:, 1], -beta[:, 2], beta[:, 3]]
    for s, lab in zip(series_b, labels_b):
        plt.semilogy(eps, s, linewidth=1.5, label=lab)
    plt.grid(True)
    plt.xlim(0.0, 1.0)
    plt.ylim(10**-0.5, 10**2)
    plt.xlabel("Excentricidade ε")
    plt.ylabel("βik")
    plt.title("Fatores βik em função da excentricidade ε")
    plt.legend()

    # k_ik x Ω (usei Ω em Hz no eixo, fica mais “físico” que o MATLAB)
    plt.figure()
    for j, lab in enumerate(["k11", "k12", "k21", "k22"]):
        plt.plot(omega_Hz, kii[:, j], linewidth=1.5, label=lab)
    plt.grid(True)
    plt.xlim(0.0, 50.0)
    plt.ylim(-1e6, 1e6)
    plt.xlabel("Velocidade angular Ω (Hz)")
    plt.ylabel("kik (N/m)")
    plt.title("Rigidezes kik em função da velocidade angular Ω")
    plt.legend()

    # c_ik x Ω
    plt.figure()
    for j, lab in enumerate(["c11", "c12", "c21", "c22"]):
        plt.plot(omega_Hz, cii[:, j], linewidth=1.5, label=lab)
    plt.grid(True)
    plt.xlim(0.0, 50.0)
    plt.ylim(-5e5, 1.5e6)
    plt.xlabel("Velocidade angular Ω (Hz)")
    plt.ylabel("cik (N·s/m)")
    plt.title("Amortecimentos cik em função da velocidade angular Ω")
    plt.legend()


# -------------------------------------------------------------
# Simulações de torque (3 casos) e órbitas a Ω constante
# -------------------------------------------------------------
def simulate_torque_cases(params, omega, kii, cii):
    # mesmos T1 e t1 do MATLAB
    t_end = 50.0
    T_list = [0.01, 0.0045, 0.004]
    theta0_list = [0.1, 0.01, 0.01]

    res = []
    for T1, th0 in zip(T_list, theta0_list):
        y0 = np.array([0.0, 0.0, 0.0, th0, 0.0, 0.0, 0.0, 0.0])
        sol = solve_ivp(
            lambda t, y: laval_mancal_rhs(t, y, params, omega, kii, cii, T1),
            (0.0, t_end),
            y0,
            method="RK45",
            max_step=0.05,
            rtol=1e-6,
            atol=1e-9,
        )
        res.append((T1, sol.t, sol.y))

    # plot simples: velocidade angular x tempo
    plt.figure()
    for T1, t, y in res:
        plt.plot(t, y[3, :], label=f"T1 = {T1:.4f} N·m")
    plt.xlabel("Tempo (s)")
    plt.ylabel("Velocidade angular θ̇ (rad/s)")
    plt.title("Passagem pela 1ª velocidade crítica – diferentes torques")
    plt.grid(True)
    plt.legend()

    # deslocamentos em y para cada caso
    plt.figure()
    for T1, t, y in res:
        plt.plot(t, 1000.0 * y[1, :], label=f"T1 = {T1:.4f} N·m")
    plt.xlabel("Tempo (s)")
    plt.ylabel("Deslocamento uy (mm)")
    plt.title("Deslocamento radial uy para diferentes torques")
    plt.grid(True)
    plt.legend()


def simulate_constant_speed_orbits(params, omega, kii, cii):
    # mesmos w do MATLAB: 30, 50, 72, 120, 150 rad/s
    w_list = [30.0, 50.0, 72.0, 120.0, 150.0]
    t_end = 100.0

    for w in w_list:
        # pega km, cm fixos para essa rotação
        km = np.array([np.interp(w, omega, kii[:, j]) for j in range(4)])
        cm = np.array([np.interp(w, omega, cii[:, j]) for j in range(4)])

        u0 = np.zeros(6)  # [uy, uz, uy_dot, uz_dot, ym, zm]
        sol = solve_ivp(
            lambda t, u: laval_mancal_const_rhs(t, u, params, w, km, cm),
            (0.0, t_end),
            u0,
            method="RK45",
            max_step=0.05,
            rtol=1e-6,
            atol=1e-9,
        )

        t = sol.t
        uy = sol.y[0, :]
        uz = sol.y[1, :]

        # ==== FILTRO DE TRANSIENTE ====
        T = 2.0 * np.pi / w          # período da rotação
        Nper = 5                     # quantos períodos finais queremos
        t_min = t_end - Nper * T     # começo da janela em regime

        mask = t >= t_min            # pega só o final
        uy_ss = uy[mask]
        uz_ss = uz[mask]
        # ==============================

        plt.figure()
        plt.plot(1000.0 * uy_ss, 1000.0 * uz_ss, linewidth=1.5)
        plt.xlabel("uy (mm)")
        plt.ylabel("uz (mm)")
        plt.title(f"Órbita no plano y–z (Ω = {w:.0f} rad/s)")
        plt.axis("equal")
        plt.grid(True)

# -------------------------------------------------------------
# main
# -------------------------------------------------------------
def main():
    params = build_parameters()
    omega_Hz, omega, eps, S_star, alfa, gama, beta, kii, cii = compute_bearing_coeffs(params)

    # parte estática
    plot_static_results(omega_Hz, omega, eps, S_star, alfa, gama, beta, kii, cii)

    # resposta temporal com torque (3 casos)
    simulate_torque_cases(params, omega, kii, cii)

    # órbitas para velocidades constantes
    simulate_constant_speed_orbits(params, omega, kii, cii)

    plt.show()


if __name__ == "__main__":
    main()
