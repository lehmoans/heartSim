# Pulsatile Membrane Chamber

A lightweight, code-first fluid–structure interaction project that builds toward the immersed-boundary ideas used in Peskin's heart-flow work. The goal is to understand and validate each numerical component on a modest computer before attempting a full heart geometry.

## Project goal

Model a two-dimensional fluid domain driven by a prescribed oscillating membrane. Start with a controlled, grid-aligned moving-wall problem to validate the fluid solver, then introduce Lagrangian membrane markers and immersed-boundary coupling.

This is a learning and numerical-methods project, not a physiological heart model. The first stages prescribe membrane motion; they do not predict cardiac muscle activation or solve a fully coupled elastic membrane.

## Learning outcomes

By completing the stages, you should be able to:

- explain incompressible Navier–Stokes and the role of pressure;
- derive a pressure-projection step and verify discrete divergence;
- calculate a prescribed periodic membrane position and velocity;
- understand the distinction between Eulerian fluid-grid values and Lagrangian membrane markers;
- implement a regularised delta kernel for velocity interpolation and force spreading;
- estimate advection and explicit-viscosity time-step restrictions;
- validate the solver with conservation checks, convergence tests, and simple limiting cases.

## Recommended implementation stages

### Stage 0 — Establish units and assumptions

Use a two-dimensional, incompressible Newtonian fluid with constant density and viscosity. Begin with a small structured grid (for example, 64 × 64). Use dimensionless or clearly stated physical parameters and keep the first model deliberately simple.

### Stage 1 — Prescribe the pulsation

For a membrane with mean position $y_0$, amplitude $A$, and frequency $f$, define

```math
y_m(t) = y_0 + A \sin(2\pi f t)
```

Its velocity is

```math
v_m(t) = 2\pi f A \cos(2\pi f t)
```

These are prescribed kinematics: they tell the solver where the membrane is and how fast it moves. They do not calculate the motion from forces.

### Stage 2 — Validate the fluid solver separately

Implement the incompressible momentum update and pressure projection before adding a moving internal membrane. For a predictor velocity $\mathbf{u}^*$, the projection solves

```math
\nabla^2 p = \frac{\rho}{\Delta t}\nabla\cdot\mathbf{u}^*
```

and corrects the velocity:

```math
\mathbf{u}^{n+1}=\mathbf{u}^* - \frac{\Delta t}{\rho}\nabla p
```

Use compatible discrete gradient, divergence, and Laplacian operators. Test that projection reduces the discrete divergence to the expected solver tolerance. The exact pressure boundary conditions depend on the chosen domain and wall treatment.

### Stage 3 — Introduce membrane motion

First use a simple grid-aligned moving-wall benchmark to verify that a prescribed wall velocity can drive the fluid. Be explicit about the geometry: a moving top wall is a useful solver test, but it is not yet an internal membrane.

Then represent the internal membrane by Lagrangian markers $\mathbf{X}(s,t)$ inside the fixed Eulerian fluid grid. The marker position is prescribed in this stage.

### Stage 4 — Add immersed-boundary coupling

Use a regularised delta kernel $\delta_h$ to transfer quantities between the moving boundary and fixed grid. Force spreading is

```math
\mathbf{F}(\mathbf{x},t)=\int \mathbf{f}(s,t)\,\delta_h(\mathbf{x}-\mathbf{X}(s,t))\,ds
```

Velocity interpolation is

```math
\frac{\partial\mathbf{X}}{\partial t}=\int \mathbf{u}(\mathbf{x},t)\,\delta_h(\mathbf{x}-\mathbf{X}(s,t))\,d\mathbf{x}
```

For prescribed membrane motion, the interpolated marker velocity must be made consistent with the target velocity. Possible numerical approaches include direct-forcing or penalty methods, but their stability and accuracy need to be analysed; force spreading and interpolation alone do not enforce a prescribed velocity automatically.

### Stage 5 — Add membrane mechanics only after the prescribed-motion case works

Replace the prescribed motion with an elastic boundary model. Begin with springs between adjacent markers, then consider bending resistance if required. Derive marker forces from an elastic energy and keep force units consistent with whether values represent force per marker or force per unit boundary length.

### Stage 6 — Validation and convergence

- Check that the prescribed position and velocity agree analytically.
- Check that the regularised kernel weights sum to the expected value and that force spreading preserves the integrated force, accounting for grid/marker spacing.
- Measure the divergence before and after projection.
- Monitor maximum velocity, pressure, boundary position, and numerical energy.
- Repeat at smaller time steps and finer grids; compare the same physical quantities.
- Check limiting cases such as zero membrane amplitude and zero forcing.
- Do not claim physiological accuracy without suitable geometry, parameters, and validation data.

## Time-step estimates

For an explicit advection scheme, a common estimate is

```math
\Delta t_{\mathrm{adv}} \le C_{\mathrm{adv}}\frac{h}{U_{\max}}
```

For forward-Euler diffusion on a two-dimensional square grid with spacing $h$,

```math
\Delta t_{\mathrm{visc}} \le \frac{h^2}{4\nu}
```

where $\nu=\mu/\rho$. These bounds are scheme-dependent and are not a complete stability guarantee. Immersed-boundary forcing, wall motion, and accuracy requirements may impose additional restrictions.

## Notebook

Open [pulsatile_membrane_numerical_work.ipynb](pulsatile_membrane_numerical_work.ipynb) and work through the derivations and small numerical checks in order. It is a teaching notebook: it establishes the kinematics, discrete operators, projection idea, and kernel checks, but it is not yet a complete production-ready Navier–Stokes/immersed-boundary solver.

## Suggested repository layout

- pulsatileMembraneChamber/README.md
- pulsatileMembraneChamber/pulsatile_membrane_numerical_work.ipynb

## Core references

- Peskin, C. S. (1977), *Numerical Analysis of Blood Flow in the Heart*, Journal of Computational Physics, 25(3), 220–252. DOI: [10.1016/0021-9991(77)90100-0](https://doi.org/10.1016/0021-9991(77)90100-0).
- The notebook derives the simplified projection and immersed-boundary concepts used for learning; consult the original paper for its specific formulation and numerical choices.
