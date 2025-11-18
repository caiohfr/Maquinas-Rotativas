# Trabalho 1 – 3-DOF Mass–Spring–Damper System

This project implements the numerical simulation of a 3-DOF mass–spring–damper system,
based on a graduate Rotating Machinery coursework assignment.

The equations of motion are written in matrix form:

\[
M \ddot{x} + C \dot{x} + K x = F(t)
\]

where:

- \(x = [x_1, x_2, x_3]^T\) are the displacements,
- \(M\) is the mass matrix,
- \(C\) is the damping matrix,
- \(K\) is the stiffness matrix,
- \(F(t)\) is a harmonic force applied to the second degree of freedom.

In the original assignment, the system is excited by a harmonic force:

\[
F(t) = [0,\; F_0 \sin(\omega t),\; 0]^T
\]

with \(F_0 = 50\) and \(\omega = 19.5\ \text{rad/s}\).

The script:

- Defines matrices \(M\), \(C\) and \(K\);
- Builds the state-space model;
- Integrates the system response in time using a numerical ODE solver;
- Plots displacements and velocities of the three DOFs.

## Requirements

- Python 3.10+
- `numpy`
- `scipy`
- `matplotlib`

Install:

```bash
pip install numpy scipy matplotlib
