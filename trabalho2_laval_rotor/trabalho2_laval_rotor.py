
## `trabalho2_laval_rotor/trabalho2_laval_rotor.py`


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
nu = 0.3     # Poisson ratio (not directly used here)
g = 9.81     # gravity (m/s^2)

# Section properties
Ie = np.pi * De**4 / 64.0          # shaft second moment of area
me = (np.pi * De**2 / 4.0) * Le * rho
md = (np.pi * Dd**2 / 4.0) * Ld * rho

# Eccentricity such that md * e = 1.5e-3 (as in the original code)
e = 1.5e-3 / md

# Equivalent linear stiffness and damping at midspan
k = 48.0 * E * Ie / (Le**3)
c = 1.0e-3 * k

# Polar inertia about rotation axis, including eccentricity contribution
Id = 0.5 * md * (Dd / 2.0) ** 2 + md * e**2


# ----------------------------------------------------------------------
# ODE definitions
# ----------------------------------------------------------------------

def laval_rhs(t, y, T):
    """
    Right-hand side for the Laval rotor with torque-driven angular acceleration.

    State vector:
        y = [theta, uy, uz, theta_dot, uy_dot, uz_dot]^T

    Equations:
        Id * theta_ddot = T + k*uz*e*cos(theta) - k*uy*e*sin(theta)
                          - c*uy_dot*e*sin(theta) + c*uz_dot*e*cos(theta)

        md * uy_ddot = md*e*(theta_ddot*sin(theta) + theta_dot^2*cos(theta))
                       - c*uy_dot - k*uy

        md * uz_ddot = md*e*(-theta_ddot*cos(theta) + theta_dot^2*sin(theta))
                       - c*uz_dot - k*uz - md*g
    """
    theta, uy, uz, theta_dot, uy_dot, uz_dot = y

    # Angular acceleration
    A = (1.0 / Id) * (
        T
        + k * uz * e * np.cos(theta)
        - k * uy * e * np.sin(theta)
        - c * uy_dot * e * np.sin(theta)
        + c * uz_dot * e * np.cos(theta)
    )

    theta_ddot = A

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
    dydt[0] = theta_dot
    dydt[1] = uy
    dydt[2] = uz
    dydt[3] = theta_ddot
    dydt[4] = uy_ddot
    dydt[5] = uz_ddot

    return dydt


def vibwconst_rhs(t, u, w):
    """
    Right-hand side for the Laval rotor with imposed constant angular speed w.

    Here we use a reduced 4-DOF system:

        u = [uy, uz, uy_dot, uz_dot]^T

    assuming theta(t) = w t and theta_dot = w (no angular dynamics).
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
    """
    Simulate the torque-driven case for a given torque T and final time t_end.
    Corresponds to the original 'ode45(@laval, [0 t_end], y0)'.
    """
    y0 = np.zeros(6)  # [theta, uy, uz, theta_dot, uy_dot, uz_dot]
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
    # T1 = 2 N·m, t1 = 1 s
    T1 = 2.0
    t1 = 1.0
    return simulate_case(T1, t1)


def simulate_case_2():
    # T2 = 0.3 N·m, t2 = 4 s
    T2 = 0.3
    t2 = 4.0
    return simulate_case(T2, t2)


def simulate_case_3():
    # T3 = 0.25 N·m, t3 = 4 s
    T3 = 0.25
    t3 = 4.0
    return simulate_case(T3, t3)


def simulate_orbit_constant_speed(w, t_final=2.0, n_points=4000):
    """
    Simulate orbital motion at constant angular speed w (rad/s),
    using vibwconst_rhs (equivalent to vibwconst_rotor in the OCTAVE code).
    """
    u0 = np.zeros(4)  # [uy, uz, uy_dot, uz_dot]
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

    # Coordinates of point G (mass center of disk)
    ug_y = uy + e * np.cos(w * t)
    ug_z = uz + e * np.sin(w * t)

    return t, uy, uz, ug_y, ug_z


# ----------------------------------------------------------------------
# Plot helpers
# ----------------------------------------------------------------------

def plot_case(t, y, title_str):
    """
    Plot displacement and rotation vs. time, similar to MATLAB plotyy.
    y has shape (6, N): rows = [theta, uy, uz, theta_dot, uy_dot, uz_dot].
    """
    theta = y[0, :]
    uy = y[1, :]
    uz = y[2, :]

    fig, axs = plt.subplots(1, 2, figsize=(12, 4))
    fig.suptitle(title_str)

    # uy x time + theta x time
    ax1 = axs[0]
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("Uy (mm)")
    ax1.plot(t, 1000.0 * uy, label="Uy (mm)")
    ax1.grid(True)

    ax1b = ax1.twinx()
    ax1b.set_ylabel("Angular velocity (rad/s)")
    ax1b.plot(t, np.gradient(theta, t), color="r", alpha=0.7, label="θ̇ (rad/s)")

    # uz x time + theta x time
    ax2 = axs[1]
    ax2.set_xlabel("Time (s)")
    ax2.set_ylabel("Uz (mm)")
    ax2.plot(t, 1000.0 * uz, label="Uz (mm)")
    ax2.grid(True)

    ax2b = ax2.twinx()
    ax2b.set_ylabel("Angular velocity (rad/s)")
    ax2b.plot(t, np.gradient(theta, t), color="r", alpha=0.7, label="θ̇ (rad/s)")

    plt.tight_layout()


def plot_orbit(uy, uz, ug_y=None, ug_z=None, w=None):
    """
    Plot orbit of the rotor center and optionally of point G.
    """
    plt.figure()
    plt.plot(1000.0 * uy, 1000.0 * uz, label="Center")
    if ug_y is not None and ug_z is not None:
        plt.plot(1000.0 * ug_y, 1000.0 * ug_z, label="Point G")
    plt.xlabel("Uy (mm)")
    plt.ylabel("Uz (mm)")
    title = "Orbit in the Y–Z plane"
    if w is not None:
        title += f" (w = {w:.1f} rad/s)"
    plt.title(title)
    plt.axis("equal")
    plt.grid(True)
    plt.legend()


# ----------------------------------------------------------------------
# Main script
# ----------------------------------------------------------------------

def main():
    # Case 1
    t1, y1 = simulate_case_1()
    plot_case(t1, y1, "Case 1 – High torque (fast crossing of 1st critical speed)")

    # Case 2
    t2, y2 = simulate_case_2()
    plot_case(t2, y2, "Case 2 – Moderate torque")

    # Case 3
    t3, y3 = simulate_case_3()
    plot_case(t3, y3, "Case 3 – Lower torque")

    # Constant-speed orbit
    w_const = 1000.0  # rad/s (example; adjust as in the original assignment)
    t_orb, uy, uz, ug_y, ug_z = simulate_orbit_constant_speed(w_const, t_final=0.5)
    plot_orbit(uy, uz, ug_y, ug_z, w=w_const)

    plt.show()


if __name__ == "__main__":
    main()
