

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# ============================================================
# Matrizes M, C, K
# ============================================================

def build_mck():
    """
    Build mass, damping and stiffness matrices for the 3-DOF system.

    Values taken from the original Octave/MATLAB assignment:
        M = diag(2, 5, 1)
        C, K as in original code.
    """
    M = np.array([
        [2.0, 0.0, 0.0],
        [0.0, 5.0, 0.0],
        [0.0, 0.0, 1.0],
    ])

    C = np.array([
        [300.0, -200.0,   0.0],
        [-200.0, 350.0, -150.0],
        [0.0,  -150.0,  150.0],
    ])

    K = np.array([
        [12000.0, -8000.0,     0.0],
        [-8000.0, 14000.0, -6000.0],
        [0.0,     -6000.0,  6000.0],
    ])

    return M, C, K


# ============================================================
# Equações de movimento no tempo (ODE)
# ============================================================

def vib_rhs(t, y, M, C, K, w_exc=19.5, F0=50.0):
    """
    Right-hand side of the 2nd-order system in state-space form.

    M x¨ + C x˙ + K x = F(t)
    F(t) = [0, F0*sin(w_exc*t), 0]^T

    State vector:
        y = [x; x_dot], where x and x_dot are 3-vectors.
    """
    n = M.shape[0]
    x = y[:n]
    xdot = y[n:]

    # Harmonic force applied at DOF 2 – real-valued version
    f_t = F0 * np.sin(w_exc * t)
    F = np.array([0.0, f_t, 0.0])

    # More numerically robust than explicit inverse
    rhs = F - C @ xdot - K @ x
    xddot = np.linalg.solve(M, rhs)

    dy = np.concatenate((xdot, xddot))
    return dy


def simulate_time_response(
    M, C, K,
    w_exc=19.5,
    F0=50.0,
    t_final=2.5,
    n_points=2000,
    x0=None,
    xdot0=None,
):
    """
    Integrate the system in time using solve_ivp.
    """
    n = M.shape[0]

    if x0 is None:
        x0 = np.zeros(n)
    if xdot0 is None:
        xdot0 = np.zeros(n)

    y0 = np.concatenate((x0, xdot0))
    t_eval = np.linspace(0.0, t_final, n_points)

    sol = solve_ivp(
        vib_rhs,
        (0.0, t_final),
        y0,
        t_eval=t_eval,
        args=(M, C, K, w_exc, F0),
        method="RK45",
        rtol=1e-6,
        atol=1e-9,
    )

    if not sol.success:
        raise RuntimeError(f"Time integration failed: {sol.message}")

    t = sol.t
    x = sol.y[:n, :]
    xdot = sol.y[n:, :]

    return t, x, xdot


# ============================================================
# FRF e resposta em frequência
# ============================================================

def compute_frf(M, C, K, w_min=0.2, w_max=200.0, dw=0.2):
    """
    Compute the Frequency Response Function:

        Z(w) = -M w^2 + j C w + K
        H(w) = Z(w)^(-1)

    Returns
    -------
    w : (N,) array
        Frequency vector [rad/s].
    H : (3, 3, N) complex array
        FRF matrix for each frequency.
    """
    w = np.arange(w_min, w_max + dw / 2.0, dw)
    n = M.shape[0]
    H = np.zeros((n, n, w.size), dtype=complex)

    for i, wi in enumerate(w):
        Z = - (wi**2) * M + 1j * wi * C + K
        H[:, :, i] = np.linalg.inv(Z)

    return w, H


def compute_steady_state_response(H, F):
    """
    Given H(w) and a force vector F, compute X(w) = H(w) * F
    for each frequency.

    Parameters
    ----------
    H : (3, 3, N) complex array
    F : (3,) real or complex array

    Returns
    -------
    Xs : (3, N) complex array
         Steady-state complex amplitudes for each DOF and frequency.
    """
    n, _, n_freq = H.shape
    Xs = np.zeros((n, n_freq), dtype=complex)
    for i in range(n_freq):
        Xs[:, i] = H[:, :, i] @ F
    return Xs


# ============================================================
# Plots
# ============================================================

def plot_displacements(t, x):
    """
    Plot displacements x1, x2, x3 vs time.
    """
    plt.figure(figsize=(8, 4))
    for i in range(x.shape[0]):
        plt.plot(t, x[i, :], label=f"x{i+1}")
    plt.title("Deslocamento de cada GDL pelo tempo")
    plt.xlabel("Tempo (s)")
    plt.ylabel("Deslocamento (m)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()


def plot_velocities(t, xdot):
    """
    Plot velocities x1_dot, x2_dot, x3_dot vs time.
    """
    plt.figure(figsize=(8, 4))
    for i in range(xdot.shape[0]):
        plt.plot(t, xdot[i, :], label=f"x{i+1}_dot")
    plt.title("Velocidade de cada GDL pelo tempo")
    plt.xlabel("Tempo (s)")
    plt.ylabel("Velocidade (m/s)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()


def plot_frf(w, H):
    """
    Plot H_ij(w) magnitudes, similar to the original Octave script.
    """
    H11 = np.abs(H[0, 0, :])
    H22 = np.abs(H[1, 1, :])
    H33 = np.abs(H[2, 2, :])
    H12 = np.abs(H[0, 1, :])
    H23 = np.abs(H[1, 2, :])
    H13 = np.abs(H[0, 2, :])

    plt.figure(figsize=(8, 4))
    plt.plot(w, H11, label="H11")
    plt.plot(w, H22, label="H22")
    plt.plot(w, H33, label="H33")
    plt.plot(w, H12, label="H12")
    plt.plot(w, H23, label="H23")
    plt.plot(w, H13, label="H13")
    plt.title("Resposta de H(ω)")
    plt.xlabel("ω (rad/s)")
    plt.ylabel("|H(ω)|")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()


def plot_frequency_response_X(w, Xs):
    """
    Plot amplitude and phase of X(w) for each DOF, similar to Octave.

    Xs : (3, N) complex array
    """
    Xs1 = np.abs(Xs[0, :])
    Xs2 = np.abs(Xs[1, :])
    Xs3 = np.abs(Xs[2, :])

    # Amplitude em mm
    plt.figure(figsize=(8, 4))
    plt.plot(w, 1000.0 * Xs1, label="X1")
    plt.plot(w, 1000.0 * Xs2, label="X2")
    plt.plot(w, 1000.0 * Xs3, label="X3")
    plt.title("Resposta de X(ω) em amplitude")
    plt.xlabel("ω (rad/s)")
    plt.ylabel("Amplitude (mm)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    # Fase em graus
    Xsf1 = np.angle(Xs[0, :], deg=True)
    Xsf2 = np.angle(Xs[1, :], deg=True)
    Xsf3 = np.angle(Xs[2, :], deg=True)

    plt.figure(figsize=(8, 4))
    plt.plot(w, Xsf1, label="X1")
    plt.plot(w, Xsf2, label="X2")
    plt.plot(w, Xsf3, label="X3")
    plt.title("Resposta de X(ω) em fase")
    plt.xlabel("ω (rad/s)")
    plt.ylabel("Fase (graus)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()


# ============================================================
# Modal analysis (opcional, mas próximo do script original)
# ============================================================

def modal_analysis(M, C, K):
    """
    Build state-space matrix A and compute eigenvalues/eigenvectors,
    similar to the Octave code.

    A = [0   I;
         -M^-1 K   -M^-1 C]
    """
    n = M.shape[0]
    Z = np.zeros((n, n))
    I = np.eye(n)

    MK = np.linalg.solve(M, K)   # M\K in Octave
    MC = np.linalg.solve(M, C)   # M\C in Octave

    A = np.block([
        [Z,         I],
        [-MK,    -MC],
    ])

    eigvals, eigvecs = np.linalg.eig(A)

    # Damping ratio and natural frequency like in Octave
    # epis(i) = ((imag/real)^2 + 1)^(-0.5)
    epis = ( (np.imag(eigvals) / np.real(eigvals))**2 + 1.0 )**(-0.5)
    w_n = -np.real(eigvals) / epis

    # Modes with damping ratio < 1 (subamortecidos)
    mask = epis < 1.0
    modos = np.where(mask)[0]
    mode_shapes = eigvecs[:, modos]

    return eigvals, w_n, epis, modos, mode_shapes


# ============================================================
# Main
# ============================================================

def main():
    # Build matrices
    M, C, K = build_mck()

    # 1) Modal analysis (opcional, mas bom pra relatório)
    eigvals, w_n, epis, modos, mode_shapes = modal_analysis(M, C, K)
    # Você pode imprimir ou salvar isso se quiser:
    # print("Autovalores:", eigvals)
    # print("Damping ratios:", epis)
    # print("Natural frequencies:", w_n)

    # 2) FRF e resposta em frequência
    w, H = compute_frf(M, C, K)
    plot_frf(w, H)

    # Força harmônica aplicada no GDL 2 (mesmo vetor do Octave)
    F_vec = np.array([0.0, 50.0, 0.0])
    Xs = compute_steady_state_response(H, F_vec)
    plot_frequency_response_X(w, Xs)

    # 3) Resposta no tempo (mesmas CIs do script Octave)
    x0 = np.array([0.01, 0.005, 0.005])
    xdot0 = np.array([0.0, 1.0, 0.0])

    t, x, xdot = simulate_time_response(
        M, C, K,
        w_exc=19.5,
        F0=50.0,
        t_final=2.5,
        n_points=2000,
        x0=x0,
        xdot0=xdot0,
    )

    plot_displacements(t, x)
    plot_velocities(t, xdot)

    plt.show()


if __name__ == "__main__":
    main()
