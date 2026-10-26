# Heart blood-flow simulation: numerical study guide

This folder is a **learning implementation inspired by** Charles S. Peskin's paper:

> C. S. Peskin, “Numerical Analysis of Blood Flow in the Heart,” *Journal of Computational Physics*, 25 (1977), 220–252. DOI: [10.1016/0021-9991(77)90100-0](https://doi.org/10.1016/0021-9991(77)90100-0).

The paper models blood, valve leaflets, and muscular heart wall in one incompressible-flow framework. Its key numerical idea is the **immersed boundary (IB) method**: solve fluid equations on a fixed Cartesian grid, represent the moving tissue with Lagrangian points and links, and exchange forces/velocities between those two descriptions using a regularised delta kernel. It combines this with a Chorin-style finite-difference/projection method for incompressible Navier–Stokes, a pressure Poisson/Laplace solve, and a nonlinear solve for the boundary configuration.

## Contents

- `peskin_1977_numerical_work.ipynb` — a guided, executable notebook that derives the equations, implements small numerical experiments, and asks you to check the results.
- This README — scope, equation map, setup, and a suggested study order.

## Important scope and fidelity note

This is a **teaching notebook, not a faithful reimplementation of every calculation or heart geometry in the 1977 paper**. It derives and demonstrates the mathematical building blocks on small grids. The code is deliberately lightweight and intended to run on an ordinary Windows PC.

The notebook covers:
1. Incompressible Navier–Stokes and dimensional/physical meaning.
2. Eulerian (fluid-grid) versus Lagrangian (boundary-point) variables.
3. Regularised delta functions: spreading boundary forces to the grid and interpolating grid velocity back to the boundary.
4. Elastic links and their force law.
5. A discrete pressure projection and the pressure Poisson equation.
6. A periodic-grid Laplacian solved with Fourier methods.
7. A small coupling demonstration, limitations, and validation questions.

The original paper's exact geometry, complete heart muscle/valve model, source/sink treatment, full time-step ordering, and published flow visualisations are **not** reproduced as a single validated simulator here. The code examples isolate components so that each can be understood and tested. The paper's 1977 discretisation and notation differ in places from modern textbook projection methods; the notebook labels simplified demonstrations accordingly.

## Numerical model at a glance

For constant-density incompressible flow, the paper uses the Navier–Stokes equations with a force density representing immersed tissue:

```math
\rho\left(\frac{\partial \mathbf u}{\partial t}+\mathbf u\cdot\nabla\mathbf u\right)
=-\nabla p+\mu\nabla^2\mathbf u+\mathbf F,\qquad
\nabla\cdot\mathbf u=0.
```

The boundary is represented by moving points $\mathbf X(s,t)$. A regularised delta kernel transfers quantities between the boundary and the fixed grid:

```math
\mathbf F(\mathbf x,t)=\int \mathbf f(s,t)\,\delta_h(\mathbf x-\mathbf X(s,t))\,ds,
```
```math
\frac{\partial\mathbf X}{\partial t}
=\int \mathbf u(\mathbf x,t)\,\delta_h(\mathbf x-\mathbf X(s,t))\,d\mathbf x.
```

The force law $\mathbf f$ depends on the boundary configuration. The paper discretises the boundary as a link structure; for a simple spring-like link, tension acts along the line joining its endpoints. At each time step, the boundary configuration and fluid velocity are updated in a coupled sequence, and a pressure solve enforces discrete incompressibility.

## Run the notebook

Recommended: Python 3.10+.

```bash
python -m venv .venv
# Windows PowerShell:
.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install numpy matplotlib jupyter
jupyter notebook peskin_1977_numerical_work.ipynb
```

The examples are intentionally small; no CFD package, GPU, or external solver is required. Run cells from top to bottom. If you use VS Code, open the notebook and select the environment where the packages were installed.

## Suggested learning route

1. Read the equations in the notebook before executing the cells.
2. Predict each result in words or on paper.
3. Run the cell and compare your prediction with the output.
4. Change one parameter at a time (grid spacing, viscosity, link stiffness, or time step).
5. Record which quantities should be conserved and test them numerically.
6. Only after the building blocks are understood, combine them into a stable two-dimensional IB solver.

## Validation checklist

- Does the discrete delta kernel have the intended support and approximately unit sum?
- Does force spreading preserve the total force to discretisation accuracy?
- Does velocity interpolation return a constant velocity unchanged?
- Do equal-and-opposite spring forces sum to zero?
- Does the projection reduce the discrete divergence?
- Does the periodic Poisson solve satisfy its equation up to solver tolerance?
- Does decreasing grid spacing or time step change the result in a systematic way?

## Primary source and further reading

- Peskin (1977), publisher record and abstract: [ScienceDirect](https://www.sciencedirect.com/science/article/pii/0021999177901000).
- A publicly indexed copy of the paper: [PDF hosted by LJLL](https://www.ljll.fr/~frey/papers/Navier-Stokes/Peskin%20C.S.%2C%20Numerical%20analysis%20of%20blood%20flow%20in%20the%20heart.pdf).
- Peskin's later review, “The Immersed Boundary Method” (2002): [Cambridge University Press](https://www.cambridge.org/core/journals/acta-numerica/article/immersed-boundary-method/95ECDAC5D1824285563270D6DD70DA9A).

Use the original paper for its full derivation, assumptions, historical discretisation, and figures. This notebook is an independent educational explanation and does not reproduce the paper's copyrighted text or figures.
