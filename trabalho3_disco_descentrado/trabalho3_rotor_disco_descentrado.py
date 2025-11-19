import numpy as np
from numpy.linalg import eig, inv
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# Dados geométricos e de material (do enunciado / MATLAB)
# ----------------------------------------------------------------------

De = 0.012   # diâmetro do eixo (m)
Le = 0.85    # comprimento do eixo (m)
Dd = 0.10    # diâmetro do disco (m)
Ld = 0.01    # comprimento do disco (m)
E = 2.0e11   # módulo de elasticidade (Pa)
rho = 7850.0 # densidade (kg/m^3)

# Segmentação do eixo
a = (2.0 / 3.0) * Le
b = (1.0 / 3.0) * Le

# Inércia de área do eixo
Ie = np.pi * De**4 / 64.0

# Massa e inércias do disco
md = (np.pi * Dd**2 / 4.0) * Ld * rho
r_disk = Dd / 2.0
Ip = 0.5 * md * r_disk**2
Id = md * r_disk**2 / 4.0 + md * Ld**2 / 12.0

# Rigidezes equivalentes (mesmas expressões do MATLAB)
K = 3.0 * E * Ie * (a**3 + b**3) / (a**3 * b**3)
k_cpl = 3.0 * E * Ie * Le * (a - b) / (a**2 * b**2)
Kr = 3.0 * E * Ie * Le / (a * b)

# Matrizes de massa e rigidez (4 GDL: uy, uz, phiy, phiz)
M = np.diag([md, md, Id, Id])

Kmat = np.array([
    [ K,      0.0,    0.0,    -k_cpl],
    [ 0.0,    K,      k_cpl,   0.0  ],
    [ 0.0,    k_cpl,  Kr,      0.0  ],
    [-k_cpl,  0.0,    0.0,     Kr   ],
])


# ----------------------------------------------------------------------
# Matrizes G(Ω) e A(Ω)
# ----------------------------------------------------------------------

def G_matrix(Omega_rad):
    """
    Matriz giroscópica G(Ω).
    No MATLAB: G(3,4) = Ip*Ω, G(4,3) = -Ip*Ω.
    """
    G = np.zeros((4, 4))
    G[2, 3] = Ip * Omega_rad
    G[3, 2] = -Ip * Omega_rad
    return G


def state_matrix(M, K, C, G):
    """
    Matriz de espaço de estados:

        z_dot = A z

    com:

        A = [[0, I],
             [-M^-1 K, -M^-1 (C+G)]]
    """
    n = M.shape[0]
    Z = np.zeros((n, n))
    I = np.eye(n)
    Minv = inv(M)
    S = C + G

    A = np.block([
        [Z,           I],
        [-Minv @ K,  -Minv @ S]
    ])
    return A


# ----------------------------------------------------------------------
# Campbell sem amortecimento
# ----------------------------------------------------------------------

def campbell_diagram(wmax_Hz=1000.0, dph=100):
    """
    Diagrama de Campbell sem amortecimento (C = 0).
    wmax_Hz: velocidade máxima (Hz)
    dph: fator de discretização tal que Omega_Hz = i/dph
    """
    max_idx = int(wmax_Hz * dph)
    n_modes = 4
    omega1 = np.zeros((n_modes, max_idx + 1))

    C0 = np.zeros((4, 4))

    for i in range(max_idx + 1):
        Omega_Hz = i / dph
        Omega_rad = 2.0 * np.pi * Omega_Hz

        G = G_matrix(Omega_rad)
        A = state_matrix(M, Kmat, C0, G)
        eigvals, _ = eig(A)

        # Mantém apenas autovalores com parte imaginária >= 0
        eigvals = eigvals[np.imag(eigvals) >= 0.0]
        omega_n = np.sort(np.abs(eigvals))   # |lambda|
        omega1[:, i] = omega_n

    return omega1


# ----------------------------------------------------------------------
# Campbell com amortecimento proporcional C = alpha K, e zeta(Ω)
# ----------------------------------------------------------------------

def campbell_with_damping(wmax_Hz=1000.0, dph=100, alpha=2.0e-4):
    max_idx = int(wmax_Hz * dph)
    n_modes = 4
    omega2 = np.zeros((n_modes, max_idx + 1))
    zeta = np.zeros((n_modes, max_idx + 1))

    C = alpha * Kmat

    for i in range(max_idx + 1):
        Omega_Hz = i / dph
        Omega_rad = 2.0 * np.pi * Omega_Hz

        G = G_matrix(Omega_rad)
        A = state_matrix(M, Kmat, C, G)
        eigvals, _ = eig(A)

        eigvals = eigvals[np.imag(eigvals) >= 0.0]

        # frequências amortecidas (parte imaginária)
        wn = np.imag(eigvals)
        omega2[:, i] = np.sort(wn)

        # fatores de amortecimento (mesma expressão do MATLAB)
        for j in range(n_modes):
            if np.real(eigvals[j]) == 0.0:
                zeta[j, i] = 0.0
            else:
                zeta[j, i] = ((np.imag(eigvals[j]) / np.real(eigvals[j]))**2 + 1.0) ** (-0.5)

    return omega2, zeta


# ----------------------------------------------------------------------
# Cálculo das FRFs H(Ω) para excitação síncrona (w = Ω)
# ----------------------------------------------------------------------

def compute_frf_vs_speed(wmax_Hz=1000.0, dph=100, C=None):
    """
    Calcula FRFs H(Ω) tal como no MATLAB:

        Z(ω) = -M ω^2 + j ω (C + G) + K

    com ω = Ω (excitação síncrona) e G = G(Ω).

    Retorna dicionário com H11, H22, H33, H44, H12, H13, H14, H23, H24, H34.
    """
    if C is None:
        C = np.zeros((4, 4))

    max_idx = int(wmax_Hz * dph)

    H11 = np.zeros(max_idx + 1, dtype=complex)
    H22 = np.zeros(max_idx + 1, dtype=complex)
    H33 = np.zeros(max_idx + 1, dtype=complex)
    H44 = np.zeros(max_idx + 1, dtype=complex)
    H12 = np.zeros(max_idx + 1, dtype=complex)
    H13 = np.zeros(max_idx + 1, dtype=complex)
    H14 = np.zeros(max_idx + 1, dtype=complex)
    H23 = np.zeros(max_idx + 1, dtype=complex)
    H24 = np.zeros(max_idx + 1, dtype=complex)
    H34 = np.zeros(max_idx + 1, dtype=complex)

    for idx in range(max_idx + 1):
        Omega_Hz = idx / dph
        Omega_rad = 2.0 * np.pi * Omega_Hz

        G = G_matrix(Omega_rad)
        S = C + G
        w = Omega_rad  # excitação síncrona

        Z = -M * w**2 + 1j * w * S + Kmat
        H = inv(Z)

        H11[idx] = H[0, 0]
        H22[idx] = H[1, 1]
        H33[idx] = H[2, 2]
        H44[idx] = H[3, 3]
        H12[idx] = H[0, 1]
        H13[idx] = H[0, 2]
        H14[idx] = H[0, 3]
        H23[idx] = H[1, 2]
        H24[idx] = H[1, 3]
        H34[idx] = H[2, 3]

    frfs = {
        "H11": H11, "H22": H22, "H33": H33, "H44": H44,
        "H12": H12, "H13": H13, "H14": H14,
        "H23": H23, "H24": H24, "H34": H34,
    }
    return frfs


# ----------------------------------------------------------------------
# Funções de plot – Campbell, zeta e FRFs
# ----------------------------------------------------------------------

def plot_campbell_pair(omega_arr, wmax_Hz=1000.0, dph=100, title_prefix="Campbell"):
    """
    Cria figura com 2 subplots:
    - esquerda: faixa completa (0–1000 Hz)
    - direita: zoom (0–40 Hz)
    com mesmas linhas ω = Ω e ω = 2Ω.
    """
    max_idx = int(wmax_Hz * dph)
    Omega_idx = np.arange(0, max_idx + 1)
    Omega_Hz = Omega_idx / dph

    freqs_Hz = omega_arr / (2.0 * np.pi)

    plt.figure()
    # Subplot 1 – full
    ax1 = plt.subplot(1, 2, 1)
    for i in range(omega_arr.shape[0]):
        ax1.plot(Omega_Hz, freqs_Hz[i, :], linewidth=1.5)
    ax1.plot(Omega_Hz, Omega_Hz, "--", linewidth=2.0, label="ω = Ω")
    ax1.plot(Omega_Hz, 2 * Omega_Hz, "--", linewidth=2.0, label="ω = 2Ω")
    ax1.set_title(f"{title_prefix}")
    ax1.set_xlabel("Velocidade angular Ω (Hz)")
    ax1.set_ylabel("Frequência ω (Hz)")
    ax1.grid(True)
    ax1.set_xlim(0.0, 500.0)
    ax1.set_ylim(0.0, 1000.0)

    # Subplot 2 – zoom
    ax2 = plt.subplot(1, 2, 2)
    for i in range(omega_arr.shape[0]):
        ax2.plot(Omega_Hz, freqs_Hz[i, :], linewidth=1.5)
    ax2.plot(Omega_Hz, Omega_Hz, "--", linewidth=2.0, label="ω = Ω")
    ax2.plot(Omega_Hz, 2 * Omega_Hz, "--", linewidth=2.0, label="ω = 2Ω")
    ax2.set_title(f"{title_prefix} (zoom)")
    ax2.set_xlabel("Velocidade angular Ω (Hz)")
    ax2.set_ylabel("Frequência ω (Hz)")
    ax2.grid(True)
    ax2.set_xlim(0.0, 500.0)
    ax2.set_ylim(0.0, 40.0)
    ax2.legend()


def plot_zeta(zeta_arr, wmax_Hz=1000.0, dph=100, title="Fator de amortecimento vs Ω"):
    max_idx = int(wmax_Hz * dph)
    Omega_idx = np.arange(0, max_idx + 1)
    Omega_Hz = Omega_idx / dph

    plt.figure()
    for i in range(zeta_arr.shape[0]):
        plt.plot(Omega_Hz, zeta_arr[i, :], linewidth=1.5)
    plt.xlabel("Velocidade angular Ω (Hz)")
    plt.ylabel("Fator de amortecimento ζ")
    plt.title(title)
    plt.grid(True)
    plt.xlim(0.0, 500.0)
    plt.ylim(0.0, 0.3)
    plt.legend([r"$\zeta_1$", r"$\zeta_2$", r"$\zeta_3$", r"$\zeta_4$"])


def plot_frf(frfs, wmax_Hz=1000.0, dph=100, title_prefix="FRF", amp_ylim=(-350, 50)):
    """
    Plota FRF: amplitude em dB (subplot esquerdo) e fase em graus (direito),
    no mesmo estilo das figuras 4 e 5 do MATLAB.
    """
    max_idx = int(wmax_Hz * dph)
    Omega_idx = np.arange(0, max_idx + 1)
    Omega_Hz = Omega_idx / dph

    # Empilha todas as FRFs para plotar
    keys_order = ["H11", "H12", "H13", "H14", "H22", "H23", "H24", "H33", "H34", "H44"]
    H_list = [frfs[k] for k in keys_order]

    # Amplitude (dB)
    plt.figure()
    ax1 = plt.subplot(1, 2, 1)
    for H in H_list:
        mag = np.abs(H)
        # mag2db: 20*log10(mag)
        with np.errstate(divide="ignore"):
            mag_db = 20.0 * np.log10(mag)
        ax1.plot(Omega_Hz, mag_db, linewidth=1.5)
    ax1.set_title(f"{title_prefix} – Amplitude")
    ax1.set_xlabel("Velocidade angular Ω (Hz)")
    ax1.set_ylabel("Amplitude (dB)")
    ax1.grid(True)
    ax1.set_xlim(0.0, 500.0)
    ax1.set_ylim(amp_ylim[0], amp_ylim[1])

    # Fase (graus)
    ax2 = plt.subplot(1, 2, 2)
    for H in H_list:
        phase_deg = np.angle(H, deg=True)
        ax2.plot(Omega_Hz, phase_deg, linewidth=1.5)
    ax2.set_title(f"{title_prefix} – Fase")
    ax2.set_xlabel("Velocidade angular Ω (Hz)")
    ax2.set_ylabel("Fase (graus)")
    ax2.grid(True)
    ax2.set_xlim(0.0, 500.0)
    ax2.set_ylim(-200.0, 200.0)

    ax1.legend(keys_order, fontsize=8)


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    # Mesmos parâmetros conceituais do MATLAB:
    # wmax_Hz = 1000, dph = 100
    # Aqui mantemos wmax_Hz = 1000 mas os plots focam em 0–500 Hz (axis).
    wmax_Hz = 1000.0
    dph = 100
    alpha = 2.0e-4

    # Campbell sem amortecimento
    omega1 = campbell_diagram(wmax_Hz=wmax_Hz, dph=dph)
    plot_campbell_pair(omega1, wmax_Hz=wmax_Hz, dph=dph,
                       title_prefix="Diagrama de Campbell - sem amortecimento")

    # Campbell com amortecimento + zeta
    omega2, zeta = campbell_with_damping(wmax_Hz=wmax_Hz, dph=dph, alpha=alpha)
    plot_campbell_pair(omega2, wmax_Hz=wmax_Hz, dph=dph,
                       title_prefix="Diagrama de Campbell - com amortecimento")
    plot_zeta(zeta, wmax_Hz=wmax_Hz, dph=dph,
              title="Fator de amortecimento conforme a rotação Ω")

    # FRFs – caso não amortecido (C = 0, S = G)
    C0 = np.zeros((4, 4))
    frfs_undamped = compute_frf_vs_speed(wmax_Hz=wmax_Hz, dph=dph, C=C0)
    plot_frf(frfs_undamped, wmax_Hz=wmax_Hz, dph=dph,
             title_prefix="FRF (sem amortecimento)", amp_ylim=(-350, 50))

    # FRFs – caso amortecido (C = alpha*K)
    C_damped = alpha * Kmat
    frfs_damped = compute_frf_vs_speed(wmax_Hz=wmax_Hz, dph=dph, C=C_damped)
    plot_frf(frfs_damped, wmax_Hz=wmax_Hz, dph=dph,
             title_prefix="FRF (com amortecimento)", amp_ylim=(-300, 0))

    plt.show()


if __name__ == "__main__":
    main()
