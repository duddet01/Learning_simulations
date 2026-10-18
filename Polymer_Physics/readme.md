# Polymer Physics: Molecular Dynamics and Langevin Simulations

![Language](https://img.shields.io/badge/language-Fortran%2090-734f96)
![Parallel](https://img.shields.io/badge/parallel-OpenMP-blue)
![Integrator](https://img.shields.io/badge/integrator-BAOAB-green)
![Field](https://img.shields.io/badge/field-polymer%20physics-orange)

A Fortran-based computational study of bead–spring polymer models using molecular dynamics (MD) and Langevin dynamics. The repository explores polymer conformations, stochastic dynamics, relaxation times, and scaling behavior for **linear** and **ring** polymers.

The project combines numerical integration, statistical mechanics, and computational polymer physics to investigate how chain length and polymer architecture influence equilibrium structure and dynamical properties.

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Physical Model](#physical-model)
- [Simulation Methods](#simulation-methods)
- [Simulation Programs](#simulation-programs)
- [Output Files](#output-files)
- [Compilation and Execution](#compilation-and-execution)
- [Analysis and Scaling Laws](#analysis-and-scaling-laws)
- [Dependencies](#dependencies)
- [Future Directions](#future-directions)
- [Motivation](#motivation)

## Overview

Polymers are many-body systems whose macroscopic properties emerge from the interactions and collective dynamics of their constituent monomers. Molecular simulations make it possible to investigate these systems beyond simplified analytical models.

This repository contains three complementary simulation programs:

| Program | Method | Main objective |
|---|---|---|
| [`MD_Simple_polymer.f90`](MD_Simple_polymer.f90) | Molecular dynamics | Polymer trajectories, energy variation, and momentum conservation |
| [`MD_Langevin_polymer_w_parallel.f90`](MD_Langevin_polymer_w_parallel.f90) | Langevin dynamics (BAOAB) | Linear-polymer conformations and relaxation times |
| [`MD_Langevin_ring_polymer_w_parallel.f90`](MD_Langevin_ring_polymer_w_parallel.f90) | Langevin dynamics (BAOAB) | Ring-polymer structure, internal relaxation, and monomer motion |

The Langevin programs use **OpenMP** to run independent simulation jobs in parallel, so multiple polymer sizes or realizations can be studied efficiently.

## Repository Structure

```text
Polymer_Physics/
├── MD_Simple_polymer.f90                      # MD of a linear bead–spring polymer
├── MD_Langevin_polymer_w_parallel.f90         # Langevin dynamics of linear polymers (BAOAB)
├── MD_Langevin_ring_polymer_w_parallel.f90    # Langevin dynamics of ring polymers (BAOAB)
└── README.md
```

## Physical Model

### 1. Bead–Spring Polymer

A polymer is represented as a chain of $N$ interacting beads connected by harmonic springs. Non-bonded beads interact through a purely repulsive excluded-volume potential. The total potential energy is

$$
U = U_{\mathrm{WCA}} + U_{\mathrm{spring}}
$$

where $U_{\mathrm{WCA}}$ describes excluded-volume interactions and $U_{\mathrm{spring}}$ describes the bonded interactions between neighboring monomers.

### 2. Weeks–Chandler–Andersen (WCA) Potential

The WCA potential is a purely repulsive, truncated and shifted Lennard–Jones potential:

$$
U_{\mathrm{WCA}}(r) =
\begin{cases}
4\varepsilon\left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right] + \varepsilon, & r < 2^{1/6}\sigma \\
0, & r \ge 2^{1/6}\sigma
\end{cases}
$$

where

- $\varepsilon$ is the interaction energy scale,
- $\sigma$ is the characteristic monomer diameter,
- $r$ is the distance between interacting monomers.

The WCA interaction prevents monomers from overlapping and models short-range excluded-volume effects without an attractive Lennard–Jones tail.

### 3. Harmonic Bond Potential

Neighboring monomers are connected through harmonic springs:

$$
U_{\mathrm{spring}} = \frac{\kappa}{2} \sum_{\langle i,j \rangle} \left(r_{ij} - r_0\right)^2
$$

where

- $\kappa$ is the spring constant,
- $r_{ij}$ is the distance between bonded monomers,
- $r_0$ is the equilibrium bond length.

For a **linear** polymer, neighboring beads are connected sequentially. For a **ring** polymer, an additional bond connects the last monomer to the first, closing the chain.

## Simulation Methods

### 1. Molecular Dynamics

The first program evolves bead positions and velocities according to Newton's equations of motion:

$$
m\,\frac{d^2\mathbf{r}_i}{dt^2} = \mathbf{F}_i, \qquad \mathbf{F}_i = -\nabla_i U
$$

It records polymer trajectories and tracks kinetic energy, potential-energy contributions, total energy, and total momentum. Without an external thermostat, these diagnostics help assess numerical stability and energy/momentum conservation.

### 2. Langevin Dynamics

The second and third programs add friction and thermal noise:

$$
m\,\frac{d\mathbf{v}_i}{dt} = \mathbf{F}_i - \gamma\,\mathbf{v}_i + \boldsymbol{\eta}_i(t)
$$

where $\gamma$ is the friction coefficient and $\boldsymbol{\eta}_i(t)$ is a stochastic force representing thermal fluctuations. For Gaussian white noise, the fluctuation–dissipation relation is

$$
\left\langle \eta_{i\alpha}(t)\,\eta_{j\beta}(t') \right\rangle = 2\gamma k_B T\,\delta_{ij}\,\delta_{\alpha\beta}\,\delta(t - t')
$$

Here $k_B$ is the Boltzmann constant and $T$ the temperature. The Langevin thermostat lets the polymer exchange energy with a thermal environment.

### 3. BAOAB Integration

Langevin dynamics is integrated with the **BAOAB** splitting scheme:

| Step | Name | Action |
|:---:|---|---|
| **B** | Half kick | Update velocities using the forces |
| **A** | Half drift | Advance positions by half a timestep |
| **O** | Ornstein–Uhlenbeck | Apply friction and stochastic thermal noise |
| **A** | Half drift | Complete the position update |
| **B** | Half kick | Update velocities using the new forces |

BAOAB is well suited to Langevin systems because it offers favorable configurational sampling properties at suitable timesteps.

### 4. Radius of Gyration

The radius of gyration measures the spatial extent of a polymer relative to its center of mass:

$$
R_g^2 = \frac{1}{N}\sum_{i=1}^{N} \left\lvert \mathbf{r}_i - \mathbf{r}_{\mathrm{CM}} \right\rvert^2
$$

where $N$ is the number of monomers, $\mathbf{r}_i$ is the position of monomer $i$, and $\mathbf{r}_{\mathrm{CM}}$ is the center-of-mass position. It is a key observable for characterizing conformations and how polymer size changes with chain length.

### 5. Rouse Relaxation Dynamics

The Rouse model describes the internal relaxation of a free-draining polymer, neglecting hydrodynamic interactions. The relaxation time of the slowest internal mode measures how quickly the polymer loses memory of its initial conformation.

For ideal Rouse dynamics, the longest relaxation time scales as

$$
\tau_R \sim N^2
$$

For free-draining dynamics with excluded-volume scaling, a commonly used estimate is

$$
\tau_R \sim N^{1+2\nu}
$$

where $\nu$ is the polymer size-scaling exponent. The Langevin programs estimate relaxation times from the autocorrelation of the slowest internal Rouse mode, allowing dynamical scaling to be investigated numerically.

## Simulation Programs

### Program 1: Basic Molecular Dynamics

**File:** [`MD_Simple_polymer.f90`](MD_Simple_polymer.f90)

Molecular dynamics of a linear bead–spring polymer.

**Features**

- Three-dimensional polymer trajectories
- WCA excluded-volume interactions
- Harmonic bonded interactions
- Kinetic and potential-energy tracking
- Total-energy variation over time
- Total-momentum tracking

**Objective:** Examine polymer motion under conservative interactions and assess numerical stability through energy and momentum diagnostics.

### Program 2: Langevin Dynamics of Linear Polymers

**File:** [`MD_Langevin_polymer_w_parallel.f90`](MD_Langevin_polymer_w_parallel.f90)

Extends the bead–spring model to Langevin dynamics and investigates how conformation and relaxation change with chain length.

**Features**

- BAOAB integration
- Friction and stochastic thermal forces
- OpenMP-parallel independent simulation jobs
- Radius-of-gyration measurements and time-dependent traces
- Slowest internal Rouse-mode autocorrelation
- Relaxation-time estimation
- Scaling analysis across polymer lengths

Used to investigate the conformational scaling relation $R_g \sim N^{\nu}$ and the dependence of relaxation time on chain length.

**Objective:** Study the relationship between polymer size, excluded-volume interactions, and internal relaxation dynamics.

### Program 3: Langevin Dynamics of Ring Polymers

**File:** [`MD_Langevin_ring_polymer_w_parallel.f90`](MD_Langevin_ring_polymer_w_parallel.f90)

Studies cyclic bead–spring polymers using Langevin dynamics and BAOAB integration. Unlike linear chains, rings have closed connectivity, which changes their internal modes and can influence conformational and dynamical properties.

**Features**

- Closed-chain ring connectivity
- WCA excluded-volume interactions
- Harmonic bonds, including the closing bond
- BAOAB Langevin integration
- OpenMP-parallel independent simulation jobs
- Radius-of-gyration traces
- Internal ring-mode autocorrelation and relaxation-time estimation
- Mean-squared displacement analysis
- Selected monomer and center-of-mass trajectories

**Objective:** Investigate the conformational statistics and relaxation behavior of ring polymers and establish a basis for comparison with linear chains.

## Output Files

All programs write plain-text numerical data for subsequent analysis and visualization.

### Basic Molecular Dynamics

| Output | Description |
|---|---|
| `polymer_trajectory.dat` | Polymer coordinates as a function of simulation time |
| `energy_momentum.dat` | Energy components, total energy, and momentum diagnostics |

### Linear-Polymer Langevin Dynamics

| Output | Description |
|---|---|
| `scaling_results3.dat` | Summary of polymer-size scaling and relaxation measurements |
| `all_runs3.dat` | Results from individual simulation runs and associated seeds |
| `Rg_trace_N*_r*.dat` | Time-dependent radius of gyration for different chain lengths and runs |
| `Rouse_AC_N*_r*.dat` | Rouse-mode autocorrelation data used for relaxation analysis |

### Ring-Polymer Langevin Dynamics

| Output | Description |
|---|---|
| `ring_scaling_results.dat` | Summary of radius-of-gyration and relaxation measurements |
| `ring_Rg_trace_N*.dat` | Radius of gyration as a function of simulation time |
| `ring_Rouse_AC_N*.dat` | Internal ring-mode autocorrelation data |
| `ring_MSD_N*.dat` | Mean-squared displacement data |
| `ring_track_traj_N*.dat` | Trajectories of selected beads and the center of mass |

> [!NOTE]
> Asterisks (`*`) stand for values inserted into filenames by the programs, such as polymer size or run index. Refer to the output headers and the corresponding Fortran `write` statements for exact column definitions.

## Compilation and Execution

A Fortran compiler such as GNU Fortran (`gfortran`) is required. OpenMP support is needed for the parallel Langevin programs.

**1. Basic MD program**

```bash
gfortran -O3 MD_Simple_polymer.f90 -o md_polymer
./md_polymer
```

**2. Linear-polymer Langevin program**

```bash
gfortran -O3 -fopenmp MD_Langevin_polymer_w_parallel.f90 -o md_langevin
OMP_NUM_THREADS=4 ./md_langevin
```

**3. Ring-polymer Langevin program**

```bash
gfortran -O3 -fopenmp MD_Langevin_ring_polymer_w_parallel.f90 -o md_ring
OMP_NUM_THREADS=4 ./md_ring
```

Replace `4` with the desired number of OpenMP threads. The commands assume the source files are in the current working directory; polymer sizes, timesteps, simulation lengths, and output configuration are set by the programs' internal parameters.

> [!WARNING]
> Output files with the same names are overwritten when a program is run again. Back up important results before starting a new run.

## Analysis and Scaling Laws

Outputs can be analyzed with Python, NumPy, pandas, and Matplotlib.

### 1. Polymer Size Scaling

For a polymer whose size follows a power law, $R_g \sim N^{\nu}$, taking logarithms gives

$$
\log R_g = \nu \log N + \text{const.}
$$

so $\nu$ is the slope of a log–log plot of $R_g$ against $N$. For a three-dimensional self-avoiding polymer in a good solvent, the asymptotic theoretical exponent is $\nu \approx 0.588$.

### 2. Relaxation-Time Scaling

The relaxation time of the slowest internal mode generally increases with chain length:

| Regime | Scaling |
|---|---|
| Ideal Rouse | $\tau_R \sim N^2$ |
| Free-draining with excluded volume | $\tau_R \sim N^{1+2\nu}$ |

A log–log plot of relaxation time against chain length yields the dynamical scaling exponent. Finite-size effects, incomplete equilibration, statistical uncertainty, and numerical parameters can all influence the measured exponent.

### 3. Mean-Squared Displacement

For a tracked monomer,

$$
\mathrm{MSD}(t) = \left\langle \left\lvert \mathbf{r}_i(t) - \mathbf{r}_i(0) \right\rvert^2 \right\rangle
$$

and for the center of mass,

$$
\mathrm{MSD}_{\mathrm{CM}}(t) = \left\langle \left\lvert \mathbf{r}_{\mathrm{CM}}(t) - \mathbf{r}_{\mathrm{CM}}(0) \right\rvert^2 \right\rangle
$$

The monomer MSD includes both internal conformational motion and overall translation, whereas the center-of-mass MSD isolates translational motion.

### 4. Numerical Diagnostics

The basic MD program provides energy and momentum diagnostics. The Langevin programs provide complementary observables: radius-of-gyration traces, mode autocorrelations, relaxation times, and MSD measurements. Together they allow numerical stability, equilibration, conformation, and dynamics to be examined.

### Example: fitting $\nu$

A minimal sketch for extracting the size exponent. Adjust the column indices to match the headers in your `scaling_results3.dat`.

```python
import numpy as np

data = np.loadtxt("scaling_results3.dat", comments="#")
N, Rg = data[:, 0], data[:, 1]          # check column order against the file header

nu, intercept = np.polyfit(np.log(N), np.log(Rg), 1)
print(f"nu ≈ {nu:.3f}  (self-avoiding walk: ≈ 0.588)")
```

## Dependencies

**Required**

- A Fortran compiler, such as GNU Fortran
- OpenMP support (for the parallel Langevin programs)

**Optional**

- Python 3
- NumPy — numerical analysis
- pandas — processing tabular output
- Matplotlib — plotting trajectories, energy variation, scaling laws, and relaxation curves

The Fortran programs do not require any Python packages to run.

## Future Directions

- [ ] Incorporate hydrodynamic interactions beyond free-draining Langevin dynamics
- [ ] Systematically compare linear and ring polymers
- [ ] Investigate finite-size effects in conformational and dynamical scaling
- [ ] Improve statistical uncertainty estimates through independent realizations
- [ ] Analyze additional Rouse modes and time-dependent correlation functions
- [ ] Develop automated trajectory visualization and post-processing
- [ ] Profile and optimize force calculations and neighbor-list methods

## Motivation

This repository explores the connection between microscopic stochastic dynamics and emergent polymer properties. By combining molecular dynamics, Langevin integration, and scaling analysis, it provides a computational framework for investigating soft matter, statistical mechanics, and the physics of interacting many-particle systems.

---

**Topics:** `polymer-physics` · `molecular-dynamics` · `langevin-dynamics` · `baoab` · `wca-potential` · `bead-spring-model` · `rouse-dynamics` · `radius-of-gyration` · `relaxation-times` · `fortran` · `openmp`
