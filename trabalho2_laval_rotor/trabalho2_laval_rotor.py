
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# ----------------------------------------------------------------------
# Parameters (from original OCTAVE code / assignment)
# ----------------------------------------------------------------------

# Geometric and material data
De = 0.012   # shaft diameter (m)
Le = 0.85    # shaft length (m)
Dd = 0.10    # disk diameter (m)
Ld = 0.015   # disk length (m)
E = 2.0e11   # Young's modulus (Pa)
rho = 7850.0 # density (kg/m^3)
g = 9.81     # gravity (m/s^2)

# Section properties and masses
Ie = np.pi * De**4 / 64.0                 # shaft second moment of area
me = (np.pi * De**2 / 4.0) * Le * rho     # shaft mass (not used directly)
md = (np.pi * Dd**2 / 4.0) * Ld * rho     # disk mass

# Eccentricity such that md * e = 1.5e-3 (as in the original code)
e = 1.5e-3 / md

# Equivalent linear stiffness and damping at midspan
k = 48.0 * E * Ie / (Le**3)
c = 1.0e-3 * k

# Polar inertia about rotation axis, including eccentricity contribution
Id = 0.5 * md * (Dd / 2.0) ** 2 + md * e**2


# ----------------------------------------------------------------------
# ODE definitions – torque-driven rotor (full 6-DOF state)
# ----------------------------------------------------------------------

def laval_rhs(t, y, T):
    """
    State:
        y = [theta, uy, uz, theta_dot, uy_dot, uz_dot]
    """
    theta, uy, uz, theta_dot, uy_dot, uz_dot = y

    # Angular acceleration (Id * theta_ddot = ...)
    theta_ddot = (1.0 / Id) * (
        T
        + k * uz * e * np.cos(theta)
        - k * uy * e * np.sin(theta)
        - c * uy_dot * e * np.sin(theta)
        + c * uz_dot * e * np.cos(theta)
    )

    # Translational accelerations
    uy_ddot = (1.0 / md) * (
        md * e * (theta_ddot * np.sin(theta) + theta_dot**2 * np.cos(theta))
        - c * uy_dot
        - k * uy
    )

    uz_ddot = (1.0 / md) * (
        md * e * (-theta_ddot * np.cos(theta) + theta_dot**2 * np.sin(theta))
        - c * uz_dot
        - k * uz
        - md * g
    )

    dydt = np.zeros_like(y)
    dydt[0] = theta_dot    # d(theta)/dt
    dydt[1] = uy_dot       # d(uy)/dt
    dydt[2] = uz_dot       # d(uz)/dt
    dydt[3] = theta_ddot   # d(theta_dot)/dt
    dydt[4] = uy_ddot      # d(uy_dot)/dt
    dydt[5] = uz_ddot      # d(uz_dot)/dt

    return dydt


# ----------------------------------------------------------------------
# ODE with imposed constant angular speed (reduced 4-DOF state)
# ----------------------------------------------------------------------

def vibwconst_rhs(t, u, w):
    """
    u = [uy, uz, uy_dot, uz_dot]
    theta(t) = w t, theta_dot = w
    """
    uy, uz, uy_dot, uz_dot = u

    uy_ddot = (1.0 / md) * (
        md * e * w**2 * np.cos(w * t)
        - c * uy_dot
        - k * uy
    )

    uz_ddot = (1.0 / md) * (
        md * e * w**2 * np.sin(w * t)
        - md * g
        - c * uz_dot
        - k * uz
    )

    dudt = np.zeros_like(u)
    dudt[0] = uy_dot
    dudt[1] = uz_dot
    dudt[2] = uy_ddot
    dudt[3] = uz_ddot

    return dudt


# ----------------------------------------------------------------------
# Simulation helpers
# ----------------------------------------------------------------------

def simulate_case(T, t_end, n_points=2000):
    """ode45(@laval, [0 t_end], [0;0;0;0;0;0])"""
    y0 = np.zeros(6)
    t_eval = np.linspace(0.0, t_end, n_points)
    sol = solve_ivp(
        lambda t, y: laval_rhs(t, y, T),
        (0.0, t_end),
        y0,
        t_eval=t_eval,
        method="RK45",
        rtol=1e-6,
        atol=1e-9,
    )
    return sol.t, sol.y


def simulate_case_1():
    return simulate_case(T=2.0, t_end=1.0)    # T1, t1


def simulate_case_2():
    return simulate_case(T=0.3, t_end=4.0)    # T2, t2


def simulate_case_3():
    return simulate_case(T=0.25, t_end=4.0)   # T3, t3


def simulate_orbit_constant_speed(w, t_final=2.0, n_points=4000, transient_fraction=0.6):
    """ode45(@vibwconst_rotor, [0 tfinal], [0;0;0;0]) + corte 60% inicial."""
    u0 = np.zeros(4)
    t_eval = np.linspace(0.0, t_final, n_points)
    sol = solve_ivp(
        lambda t, u: vibwconst_rhs(t, u, w),
        (0.0, t_final),
        u0,
        t_eval=t_eval,
        method="RK45",
        rtol=1e-6,
        atol=1e-9,
    )

    t = sol.t
    uy = sol.y[0, :]
    uz = sol.y[1, :]

    # Ponto G (centro de massa do disco)
    ug_y = uy + e * np.cos(w * t)
    ug_z = uz + e * np.sin(w * t)

    # Mantém só os últimos 40% (como i = 6*L/10 : L)
    i0 = int(transient_fraction * len(t))
    return t[i0:], uy[i0:], uz[i0:], ug_y[i0:], ug_z[i0:]


# ----------------------------------------------------------------------
# Plot helpers
# ----------------------------------------------------------------------

def plot_case(t, y, title_str, disp_limit_mm):
    theta = y[0, :]
    uy = y[1, :]
    uz = y[2, :]
    theta_dot = y[3, :]

    fig, axs = plt.subplots(1, 2, figsize=(12, 4))
    fig.suptitle(title_str)

    # Y
    ax1 = axs[0]
    ax1.set_xlabel("Tempo (s)")
    ax1.set_ylabel("Deslocamento do centro em Y (mm)")
    ax1.plot(t, 1000.0 * uy, label="Uy (mm)")
    ax1.grid(True)
    ax1.set_xlim(t[0], t[-1])
    ax1.set_ylim(-disp_limit_mm, disp_limit_mm)

    ax1b = ax1.twinx()
    ax1b.set_ylabel("Rotação (rad/s)")
    ax1b.plot(t, theta_dot, color="r", alpha=0.7, label="θ̇ (rad/s)")

    # Z
    ax2 = axs[1]
    ax2.set_xlabel("Tempo (s)")
    ax2.set_ylabel("Deslocamento do centro em Z (mm)")
    ax2.plot(t, 1000.0 * uz, label="Uz (mm)")
    ax2.grid(True)
    ax2.set_xlim(t[0], t[-1])
    ax2.set_ylim(-disp_limit_mm, disp_limit_mm)

    ax2b = ax2.twinx()
    ax2b.set_ylabel("Rotação (rad/s)")
    ax2b.plot(t, theta_dot, color="r", alpha=0.7, label="θ̇ (rad/s)")

    plt.tight_layout()


def plot_orbit(uy, uz, ug_y, ug_z, w=None):
    plt.figure(figsize=(5, 5))
    plt.plot(1000.0 * uy, 1000.0 * uz, label="W (center)")
    plt.plot(1000.0 * ug_y, 1000.0 * ug_z, label="G (disk center)")
    plt.xlabel("uy (mm)")
    plt.ylabel("uz (mm)")
    title = "Órbita no plano Y–Z"
    if w is not None:
        title += f" (w = {w:.1f} rad/s)"
    plt.title(title)
    plt.axis("equal")
    plt.grid(True)
    plt.legend()
    plt.xlim(-2.0, 2.0)
    plt.ylim(-2.0, 2.0)


def compare_amplitudes(uy_orbit, w):
    # Teórico (Kramer) em mm
    r_theoretical_mm = 1000.0 * (md * e * w**2) / np.sqrt(
        (k - md * w**2)**2 + (c * w)**2
    )
    # Numérico (máx |uy|)
    r_numeric_mm = np.max(1000.0 * np.abs(uy_orbit))
    return r_theoretical_mm, r_numeric_mm


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    # Caso 1
    t1, y1 = simulate_case_1()
    plot_case(t1, y1, "Caso 1 – Alto torque (cruzamento rápido da 1ª crítica)", disp_limit_mm=6.0)

    # Caso 2
    t2, y2 = simulate_case_2()
    plot_case(t2, y2, "Caso 2 – Torque moderado", disp_limit_mm=10.0)

    # Caso 3
    t3, y3 = simulate_case_3()
    plot_case(t3, y3, "Caso 3 – Torque reduzido", disp_limit_mm=12.0)

    # Órbita com velocidade constante (w = 1000 rad/s, tfinal = 2 s)
    w_const = 1000.0
    t_orb, uy, uz, ug_y, ug_z = simulate_orbit_constant_speed(w_const, t_final=2.0)
    plot_orbit(uy, uz, ug_y, ug_z, w=w_const)

    # Comparação de amplitudes (Kramer x numérico)
    r_teo, r_num = compare_amplitudes(uy, w_const)
    print(f"Amplitude teórica (Kramer) em Y: {r_teo:.3f} mm")
    print(f"Amplitude numérica (ode) em Y:   {r_num:.3f} mm")

    plt.show()


if __name__ == "__main__":
    main()
