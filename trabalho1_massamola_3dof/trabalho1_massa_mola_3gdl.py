

## 🧾 Código Python do Trabalho 1  
# trabalho1_massa_mola_3gdl.py`


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


def build_mck():
    """
    Build mass, damping and stiffness matrices for the 3-DOF system.

    These values are taken from the original Octave/MATLAB assignment:
        M = diag(2, 5, 1)
        C and K as in the report/code.
    """
    M = np.array([
        [2.0,   0.0,   0.0],
        [0.0,   5.0,   0.0],
        [0.0,   0.0,   1.0],
    ])

    C = np.array([
        [300.0,  -200.0,    0.0],
        [-200.0,  350.0, -150.0],
        [0.0,    -150.0,  150.0],
    ])

    K = np.array([
        [12000.0,  -8000.0,      0.0],
        [-8000.0,  14000.0,  -6000.0],
        [0.0,      -6000.0,   6000.0],
    ])

    return M, C, K


def vib_rhs(t, y, M, C, K, w_exc=19.5, F0=50.0):
    """
    Right-hand side of the state-space model equivalent to the Octave function vib(t,y).

    Original idea (in Octave):
        M x¨ + C x˙ + K x = F(t)
    with:
        F(t) = [0; 50*exp(j*w*t); 0]  (we use the real part here)

    State vector:
        y = [x; x_dot], where x and x_dot are 3-vectors.
    """
    n = M.shape[0]
    x = y[:n]
    xdot = y[n:]

    # Harmonic force applied at DOF 2 – real-valued version
    f_t = F0 * np.sin(w_exc * t)
    F = np.array([0.0, f_t, 0.0])

    Minv = np.linalg.inv(M)
    xddot = Minv @ (F - C @ xdot - K @ x)

    dy = np.concatenate((xdot, xddot))
    return dy


def simulate_time_response(
    M, C, K,
    w_exc=19.5,
    F0=50.0,
    t_final=5.0,
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

    t = sol.t
    x = sol.y[:n, :]
    xdot = sol.y[n:, :]

    return t, x, xdot


def plot_displacements(t, x):
    """
    Plot displacements x1, x2, x3 vs time.
    """
    plt.figure()
    for i in range(x.shape[0]):
        plt.plot(t, x[i, :], label=f"x{i+1}")
    plt.title("Displacement of each DOF vs. time")
    plt.xlabel("Time (s)")
    plt.ylabel("Displacement (m)")
    plt.grid(True)
    plt.legend()


def plot_velocities(t, xdot):
    """
    Plot velocities x1_dot, x2_dot, x3_dot vs time.
    """
    plt.figure()
    for i in range(xdot.shape[0]):
        plt.plot(t, xdot[i, :], label=f"x{i+1}_dot")
    plt.title("Velocity of each DOF vs. time")
    plt.xlabel("Time (s)")
    plt.ylabel("Velocity (m/s)")
    plt.grid(True)
    plt.legend()


def main():
    # Build matrices
    M, C, K = build_mck()

    # Time response
    t, x, xdot = simulate_time_response(M, C, K, w_exc=19.5, F0=50.0, t_final=5.0)

    # Plots (equivalente aos gráficos do Octave)
    plot_displacements(t, x)
    plot_velocities(t, xdot)

    plt.show()


if __name__ == "__main__":
    main()
