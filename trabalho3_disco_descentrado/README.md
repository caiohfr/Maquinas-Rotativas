# Trabalho 3 – Frequency-Domain Analysis of a Rotor with Eccentric Disk

This project implements in Python the frequency-domain analysis of a
single-disk rotor with eccentricity, based on the MATLAB script of
"Trabalho 3 – Análise no domínio da frequência".

The model is a 4-DOF rotor:

- `uy, uz` – lateral displacements at the disk location;
- `phiy, phiz` – rotations about z and y axes respectively.

The equations of motion, in matrix form, are:

\[
M \ddot{q} + (C + G(\Omega)) \dot{q} + K q = F(t)
\]

where:

- \(q = [u_y, u_z, \phi_y, \phi_z]^T\)
- \(M\) is the mass/inertia matrix (disk translational and rotational inertia);
- \(K\) is the lateral/rotational stiffness matrix of the shaft;
- \(C\) is the viscous damping matrix (often proportional to \(K\));
- \(G(\Omega)\) is the gyroscopic matrix, proportional to the polar inertia and rotor speed.

The script reproduces:

1. **Campbell diagram (undamped):**
   - Natural frequencies vs. spin speed Ω
   - Forward and backward modes

2. **Campbell diagram (with proportional damping):**
   - Damped natural frequencies vs. Ω
   - Modal damping ratios vs. Ω

3. **Frequency response functions (FRFs):**
   - For selected input–output DOF pairs
   - As a function of excitation frequency at fixed spin speeds.

## Implemented features

The Python script includes:

- Construction of matrices `M`, `K`, `C` from the geometric and material data:
  - Shaft segment lengths `a`, `b`;
  - Shaft diameter, Young’s modulus and second moment of area `Ie`;
  - Disk translational mass `md` and inertias `Id` (transverse), `Ip` (polar).
- Definition of a gyroscopic matrix `G(Omega)` proportional to `Ip * Omega`;
- State-space assembly:

\[
\dot{z} = A(\Omega) z
\]

with:

\[
A(\Omega) =
\begin{bmatrix}
0 & I \\
- M^{-1}K & -M^{-1}(C + G(\Omega))
\end{bmatrix}
\]

- Functions to compute:
  - `campbell_diagram()` – undamped frequencies vs. Ω;
  - `campbell_with_damping()` – damped frequencies and ζ vs. Ω;
  - `compute_frf()` – FRF matrix \(H(j\omega)\) for a given spin speed.

## Requirements

- Python 3.10+
- `numpy`
- `matplotlib`

Optionally, `scipy` if you extend the script for time-domain simulation.

Install:

```bash
pip install numpy matplotlib
