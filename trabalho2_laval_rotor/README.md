# Trabalho 2 – Laval Rotor with Different Acceleration Profiles

This project implements in Python the numerical study of a Laval rotor
subjected to different torque/acceleration profiles, based on the original
OCTAVE script.

The model considers:

- A rigid disk with mass `md` mounted at the mid-span of a flexible shaft;
- Shaft flexibility represented by an equivalent lateral stiffness `k`;
- Viscous damping `c` at the disk location;
- Unbalance `e` (mass eccentricity);
- Gravity acting in the vertical direction.

The main equations (in the disk-fixed coordinates) are:

\[
I_d \ddot{\theta} =
T + k\,u_z\,e\cos\theta - k\,u_y\,e\sin\theta
- c\,\dot{u}_y e\sin\theta + c\,\dot{u}_z e\cos\theta
\]

\[
m_d \ddot{u}_y =
m_d e\left(\ddot{\theta}\sin\theta + \dot{\theta}^2\cos\theta\right)
- c\,\dot{u}_y - k\,u_y
\]

\[
m_d \ddot{u}_z =
m_d e\left(-\ddot{\theta}\cos\theta + \dot{\theta}^2\sin\theta\right)
- c\,\dot{u}_z - k\,u_z - m_d g
\]

where:

- \(\theta\) is the angular position of the disk,
- \(u_y, u_z\) are lateral displacements at the disk center,
- \(T\) is the applied torque.

The original assignment simulates three cases:

1. **Case 1:** torque `T1 = 2` N·m, time span `t1 = 1 s`  
2. **Case 2:** torque `T2 = 0.3` N·m, time span `t2 = 4 s`  
3. **Case 3:** torque `T3 = 0.25` N·m, time span `t3 = 4 s`

and also analyzes:

- The steady-state response for imposed constant angular speed \(\omega\);
- The orbit of the rotor center and of a generic point \(W\) relative to the mass center \(G\);
- Comparison of the amplitude with Kramer’s analytical expression.

## Implemented features

The Python script:

- Defines the geometric and material properties of the shaft and disk;
- Computes:
  - Shaft mass, disk mass and eccentricity;
  - Shaft flexural stiffness `k` and viscous damping `c`;
  - Disk polar inertia `Id`;
- Implements the ODE right-hand sides:
  - `laval_rhs(t, y, T)` – torque-driven angular acceleration;
  - `vibwconst_rhs(t, u, w)` – motion with imposed constant speed `w`;
- Integrates the equations in time using `scipy.integrate.solve_ivp`;
- Provides helper functions to:
  - Simulate each case (`simulate_case_1/2/3`);
  - Simulate constant-speed response (`simulate_orbit_constant_speed`);
  - Plot time histories and orbits.

## Requirements

- Python 3.10+
- `numpy`
- `scipy`
- `matplotlib`

Install:

```bash
pip install numpy scipy matplotlib
