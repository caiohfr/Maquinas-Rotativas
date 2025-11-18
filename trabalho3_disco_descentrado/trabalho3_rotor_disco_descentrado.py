
## `trabalho3_rotor_disco_descentrado/trabalho3_rotor_disco_descentrado.py`


import numpy as np
from numpy.linalg import eig
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# Geometric and material data (from the original assignment)
# ----------------------------------------------------------------------

De = 0.012   # shaft diameter (m)
Le = 0.85    # shaft length (m)
Dd = 0.10    # disk diameter (m)
Ld = 0.01    # disk length (m)
E = 2.0e11   # Young's modulus (Pa)
rho = 7850.0 # density (kg/m^3)
nu = 0.3     # Poisson ratio (not directly used here)

# Shaft segmentation (same structure as original code: lengths a and b)
a = (2.0 / 3.0) * Le
b = (1.0 / 3.0) * Le

Ie = np.pi * De**4 / 64.0

# Disk mass and inertias
md = (np.pi * Dd**2 / 4.0) * Ld * rho
r_disk = Dd / 2.0
Ip = 0.5 * md * r_disk**2
Id = md * r_disk**2 / 4.0 + md * Ld**2 / 12.0

# From the original formulation (see your MATLAB code):
K = 3.0 * E * Ie * (a**3 + b**3) / (a**3 * b**3)
k_cpl = 3.0 * E * Ie * Le * (a - b) / (a**2 * b**2)
Kr = 3.0 * E * Ie * Le / (a * b)

# Mass and stiffness matrices (4 DOF: uy, uz, phiy, phiz)
M = np.diag([md, md, Id, Id])

Kmat = np.array([
    [ K,     0.0,    0.0,   -k_cpl],
    [ 0.0,   K,      k_cpl,  0.0  ],
    [ 0.0,   k_cpl,  Kr,     0.0  ],
    [-k_cpl, 0.0,    0.0,    Kr   ],
])


def G_matrix(Omega):
    """
    Gyroscopic matrix G(Ω).
    In the original MATLAB, G(3,4) = Ip*Ω and G(4,3) = -Ip*Ω.
    """
    G = np.zeros((4, 4))
    G[2, 3] = Ip * Omega
    G[3, 2] = -Ip * Omega
    return G


def state_matrix(M, K, C, G):
    """
    Build state-space matrix A(Ω):

        z_dot = A z

    with:

        A = [[0, I],
             [-M^-1 K, -M^-1 (C+G)]]
    """
    n = M.shape[0]
    Z = np.zeros((n, n))
    I = np.eye(n)
    Minv = np.linalg.inv(M)
    S = C + G

    A = np.block([
        [Z,          I],
        [-Minv @ K, -Minv @ S]
    ])
    return A


def campbell_diagram(wmax=1000.0, dph=100):
    """
    Computes undamped Campbell diagram, similar to the original MATLAB code.

    Omega1 = 0:dph:wmax*dph  (Hz)
    w = Omega1/dph           (Hz)
    """
    max_Omega_idx = int(wmax * dph)
    n_modes = 4
    omega1 = np.zeros((n_modes, max_Omega_idx + 1))

    C0 = np.zeros((4, 4))

    for i in range(max_Omega_idx + 1):
        Omega_Hz = i / dph
        Omega_rad = 2.0 * np.pi * Omega_Hz

        G = G_matrix(Omega_rad)
        A = state_matrix(M, Kmat, C0, G)
        eigvals, _ = eig(A)

        # keep eigenvalues with Im(lambda) >= 0
        eigvals = eigvals[np.imag(eigvals) >= 0.0]
        omega_n = np.sort(np.abs(eigvals))  # |lambda|
        omega1[:, i] = omega_n

    return omega1


def campbell_with_damping(wmax=1000.0, dph=100, alpha=2.0e-4):
    """
    Computes Campbell diagram with proportional damping C = alpha*K,
    and extracts modal damping ratios ζ(Ω).
    """
    max_Omega_idx = int(wmax * dph)
    n_modes = 4
    omega2 = np.zeros((n_modes, max_Omega_idx + 1))
    zeta = np.zeros((n_modes, max_Omega_idx + 1))

    C = alpha * Kmat

    for i in range(max_Omega_idx + 1):
        Omega_Hz = i / dph
        Omega_rad = 2.0 * np.pi * Omega_Hz

        G = G_matrix(Omega_rad)
        A = state_matrix(M, Kmat, C, G)
        eigvals, _ = eig(A)

        eigvals = eigvals[np.imag(eigvals) >= 0.0]

        # Damped frequencies and damping factors from eigenvalues lambda
        wn = np.imag(eigvals)
        omega2[:, i] = np.sort(wn)

        for j in range(n_modes):
            if np.real(eigvals[j]) == 0.0:
                zeta[j, i] = 0.0
            else:
                zeta[j, i] = (
                    (np.imag(eigvals[j]) / np.real(eigvals[j])) ** 2 + 1.0
                ) ** (-0.5)

    return omega2, zeta


def compute_frf(Omega_Hz, freq_Hz):
    """
    Compute FRF matrix H(jω) = Z^-1, with:

        Z(ω) = -M ω^2 + j ω (C + G) + K

    for a given spin speed Ω (Hz) and an array of excitation frequencies freq_Hz.
    """
    C = np.zeros((4, 4))
    Omega_rad = 2.0 * np.pi * Omega_Hz
    G = G_matrix(Omega_rad)

    freq_Hz = np.asarray(freq_Hz)
    H_all = np.zeros((len(freq_Hz), 4, 4), dtype=complex)

    for i, f in enumerate(freq_Hz):
        w = 2.0 * np.pi * f
        Z = -M * w**2 + 1j * w * (C + G) + Kmat
        H_all[i, :, :] = np.linalg.inv(Z)

    return H_all


def plot_campbell(omega_arr, wmax=1000.0, dph=100, title="Campbell diagram", zoom=False):
    Omega_idx = np.arange(0, int(wmax * dph) + 1)
    Omega_Hz = Omega_idx / dph

    plt.figure()
    for i in range(omega_arr.shape[0]):
        plt.plot(Omega_Hz, omega_arr[i, :] / (2.0 * np.pi), linewidth=1.5)
    plt.plot(Omega_Hz, Omega_Hz, "--", linewidth=2.0, label="ω = Ω")
    plt.plot(Omega_Hz, 2 * Omega_Hz, "--", linewidth=2.0, label="ω = 2Ω")
    plt.xlabel("Spin speed Ω (Hz)")
    plt.ylabel("Frequency ω (Hz)")
    plt.title(title)
    plt.grid(True)
    plt.legend()
    if zoom:
        plt.ylim(0.0, 40.0)


def plot_zeta(zeta_arr, wmax=1000.0, dph=100, title="Modal damping vs speed"):
    Omega_idx = np.arange(0, int(wmax * dph) + 1)
    Omega_Hz = Omega_idx / dph

    plt.figure()
    for i in range(zeta_arr.shape[0]):
        plt.plot(Omega_Hz, zeta_arr[i, :], linewidth=1.5)
    plt.xlabel("Spin speed Ω (Hz)")
    plt.ylabel("Damping ratio ζ")
    plt.title(title)
    plt.grid(True)
    plt.ylim(0.0, 0.3)


def main():
    # Parameters for Campbell diagram (can be tuned to match your original script)
    wmax = 100.0  # maximum frequency index
    dph = 10      # resolution (Ω steps in Hz)

    # Undamped Campbell
    omega1 = campbell_diagram(wmax=wmax, dph=dph)
    plot_campbell(omega1, wmax=wmax, dph=dph, title="Campbell diagram (undamped)", zoom=True)

    # Damped Campbell
    omega2, zeta = campbell_with_damping(wmax=wmax, dph=dph, alpha=2.0e-4)
    plot_campbell(omega2, wmax=wmax, dph=dph, title="Campbell diagram (with damping)", zoom=True)
    plot_zeta(zeta, wmax=wmax, dph=dph)

    plt.show()


if __name__ == "__main__":
    main()
