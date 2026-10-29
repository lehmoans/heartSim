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


## Extended notebook: end-to-end flow benchmark and plots

The notebook now includes an educational oscillating-lid cavity prototype that advances a simplified 2D incompressible Navier–Stokes flow and plots:

- final speed-magnitude field with velocity vectors;
- horizontal and vertical velocity-component fields, $u(x,y)$ and $v(x,y)$;
- velocity streamlines;
- relative-pressure contour field with velocity vectors overlaid;
- velocity-vector field coloured by speed magnitude;
- prescribed wall speed over time;
- interior RMS divergence over time;
- maximum fluid speed over time.

This moving-lid cavity is an intermediate benchmark for testing the fluid update, pressure solve, wall conditions, and diagnostics. It is **not** the final internal pulsatile membrane simulation and is not a faithful reproduction of Peskin's complete heart-flow solver. The notebook labels the collocated finite-difference projection scheme as educational and lists its remaining numerical validation work.

## Remaining work for the internal membrane

Before treating the project as a complete internal-membrane simulation:

1. Upgrade or independently validate the fluid discretisation with compatible pressure/velocity operators.
2. Choose physically consistent chamber geometry and inlet/outlet or volume-preserving motion.
3. Add internal Lagrangian membrane markers and a documented direct-forcing/penalty coupling method.
4. Track marker slip, pressure-solver residual, discrete divergence, integrated force, and volume flux.
5. Generate spatial plots of both velocity components, speed magnitude, velocity vectors/streamlines, pressure contours with velocity vectors overlaid, membrane position, marker slip, divergence, and flux over at least one cycle.
6. Demonstrate time-step and grid-refinement convergence.
7. Add elastic membrane mechanics only after the prescribed-motion case passes the tests.

A simulation should not be considered complete merely because it produces smooth plots. Record the boundary conditions, parameters, tolerances, and convergence results.


## Standalone Python script

[`model.py`](model.py) consolidates the notebook's code cells and their comments into one Python script, in notebook execution order. It includes the analytical kinematics example, spectral projection test, immersed-boundary kernel and force-spreading checks, interpolation example, time-step estimates, oscillating-lid flow prototype, and plotting/diagnostic code. Install `numpy` and `matplotlib`, then run it with `python model.py`.

The script runs these educational examples sequentially and later sections reuse some variable names. It remains a learning prototype; the oscillating-lid scheme has not been independently validated and is not yet a full internal-membrane simulation.


## Stage 2 prototype: prescribed immersed-boundary coupling

The new [stage2_immersed_boundary_demo.py](stage2_immersed_boundary_demo.py) integrates several previously separate ideas in one small experiment:

1. a periodic 2D incompressible fluid update with advection and viscosity;
2. a moving set of Lagrangian markers with prescribed sinusoidal motion;
3. regularised-delta interpolation of fluid velocity to the markers;
4. penalty feedback to calculate marker forces;
5. force spreading back to the Eulerian grid;
6. pressure projection and diagnostics.

Run it with:

```bash
python stage2_immersed_boundary_demo.py
```

Dependencies remain NumPy and Matplotlib. The demonstration uses a periodic line translating tangentially, not a closed membrane that changes chamber volume. It is intended to demonstrate the coupling data path, not to represent pumping or validated heart flow. Its numerical behaviour has not yet been independently verified; complete the Stage 2 checklist before interpreting its plots. In particular, study pressure-solver convergence, discrete divergence, time-step sensitivity, penalty-gain sensitivity, and force-transfer consistency.

## Next implementation milestones

1. **Verify Stage 2:** test the kernel and transfer operators independently; record divergence and pressure residuals; vary the time step, pressure iteration count, grid resolution, and penalty gain.
2. **Move to a closed membrane:** choose a geometry and prescribed deformation that is compatible with incompressible volume conservation and the domain boundary conditions.
3. **Add chamber walls and ports:** define inlet/outlet boundary conditions and measure net volume flux over a cycle.
4. **Add membrane mechanics:** derive elastic marker forces from a documented energy or force law.
5. **Close the fluid–structure feedback loop:** calculate membrane forces from deformation and fluid traction, then check coupling stability and conservation.
6. **Validate before claiming pump performance:** compare against analytical limiting cases or published benchmarks and demonstrate grid/time-step convergence.
