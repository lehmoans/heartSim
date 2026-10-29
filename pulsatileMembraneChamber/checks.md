# Pulsatile Membrane Chamber — Verification Checklist

## Current scope warning

The current model.py is an educational numerical prototype containing prescribed membrane kinematics, immersed-boundary helper-function experiments, and an oscillating-lid cavity / pressure–velocity projection prototype. These components are **not yet evidence of a validated internal pulsatile-membrane simulation**.

Keep the tests below separated by component. Passing a helper-function test does not mean the complete coupled solver is correct.

## Status key

- [ ] Not started
- [x] Completed and supported by recorded evidence
- [!] Failed, blocked, or needs investigation

## A. Establish the baseline

- [ ] Read each section of model.py and identify its inputs, outputs, equations, units, and assumptions.
- [ ] Run python model.py from a fresh terminal using the documented dependencies.
- [ ] Record Python, NumPy, and Matplotlib versions and the exact run command.
- [ ] Confirm execution finishes without exceptions or unexpected NaN/infinity values.
- [ ] Identify which plots/results belong to the cavity prototype and which belong to standalone immersed-boundary helper tests.
- [ ] Record the current grid, time step, fluid properties, geometry, and membrane-motion parameters.

**Do not interpret successful execution or visually smooth plots as proof of accuracy.**

## B. Prescribed membrane kinematics

- [ ] Write down the prescribed displacement function and its amplitude, frequency, phase, and units.
- [ ] Check displacement at known phases (for example, zero crossing and maximum displacement).
- [ ] Check the computed membrane velocity against the analytical time derivative of displacement.
- [ ] Check periodicity over one complete cycle.
- [ ] Confirm amplitude and frequency changes affect the kinematics as predicted.
- [ ] Check that the chosen time step resolves the motion adequately.

## C. Regularized delta kernel

- [ ] Plot the kernel over its support and confirm its expected symmetry.
- [ ] Check that the kernel is zero outside its documented support.
- [ ] Numerically test the discrete kernel's normalization / partition-of-unity behaviour on the intended grid.
- [ ] Test points near the domain boundary separately; truncated support can change normalization.
- [ ] Confirm grid spacing and kernel scaling are consistent with the formulation.
- [ ] Compare against the exact kernel formula used as the reference.

## D. Force spreading and velocity interpolation

- [ ] Test spreading using a single marker and a simple known force.
- [ ] Check units and grid-spacing factors in the spread force.
- [ ] Perform the documented force-conservation check and report its error.
- [ ] Test interpolation on a constant velocity field; the interpolated velocity should reproduce the constant where the kernel is properly supported and normalized.
- [ ] Test interpolation on a simple linear field and quantify the error.
- [ ] Check marker/grid indexing and behaviour near boundaries.
- [ ] Where the formulation requires it, test the discrete work/adjoint relationship between spreading and interpolation.

## E. Pressure–velocity projection prototype

- [ ] Write down the discrete divergence operator, pressure-gradient operator, and pressure Poisson equation actually implemented.
- [ ] Check sign conventions and boundary conditions for each operator.
- [ ] Use a manufactured velocity field with known divergence to test the divergence calculation.
- [ ] Verify that the pressure correction reduces the velocity-field divergence.
- [ ] Record the divergence norm before and after projection using a defined metric and tolerance.
- [ ] Check pressure boundary conditions and the pressure reference / gauge.
- [ ] Test a simple flow with a known solution before relying on the cavity plots.
- [ ] Check whether the grid arrangement and discrete operators are compatible; a low pressure residual alone is insufficient.

## F. Oscillating-lid cavity prototype

- [ ] Confirm the geometry and lid/boundary velocities match the documented setup.
- [ ] Check no-slip or other intended boundary conditions on every wall.
- [ ] Confirm the flow responds in the expected direction when the lid motion is reversed.
- [ ] Check a zero-lid-velocity case; the expected stationary solution should be recovered under the stated conditions.
- [ ] Compare against a trusted benchmark or published reference using matching parameters.
- [ ] Check grid and time-step refinement for velocity and pressure quantities of interest.
- [ ] Confirm pressure plots use a clearly stated reference value and are not interpreted as absolute pressure unless justified.

## G. Time integration and stability

- [ ] Derive or cite the stability restrictions relevant to the implemented advection, diffusion, and forcing terms.
- [ ] Record the relevant dimensionless time-step measures (for example, advective CFL and viscous diffusion restriction where applicable).
- [ ] Run a smaller time step and compare the same quantities of interest.
- [ ] Monitor kinetic energy, maximum velocity, divergence, and pressure correction over time.
- [ ] Investigate unexplained energy growth, oscillation, drift, or dependence on the initial state.

## H. Coupling milestone — required before calling this a membrane-fluid simulation

- [ ] Define the internal membrane geometry and its motion / force law.
- [ ] Define how membrane forces are spread to the fluid grid.
- [ ] Define how fluid velocity is interpolated back to the membrane.
- [ ] Confirm the fluid update and membrane update are coupled at the intended time level.
- [ ] Test action–reaction / force balance and relevant conservation properties.
- [ ] Demonstrate a simple limiting case with a known or independently predictable response.
- [ ] Show that the membrane actually influences the internal fluid solution, rather than merely being plotted alongside a separate cavity flow.
- [ ] Verify grid/time-step sensitivity of pressure difference, flow rate, and membrane motion.
- [ ] Only after these checks, compare against experimental data or a suitable literature benchmark and discuss physical validation.

## I. Minimum evidence for each test

For every test, record:

1. **Question:** what property is being checked?
2. **Prediction:** what should happen, and why?
3. **Setup:** parameters, grid, time step, and boundary conditions.
4. **Metric:** a numerical error or diagnostic, with units and tolerance.
5. **Result:** measured value and pass/fail.
6. **Interpretation:** what this test does and does not establish.

## Test log

| Date | Test / case | Expected result or tolerance | Measured result / error | Pass? | Evidence / notes |
|---|---|---|---|---|---|
| | | | | | |


## J. Stage 2 immersed-boundary coupling demo

File: `stage2_immersed_boundary_demo.py`

This is a small integration experiment using a periodically translating line of markers. It is **not** a closed pulsatile chamber.

- [ ] Run the script from a fresh terminal and record versions, run command, and any warnings.
- [ ] Check that zero membrane amplitude produces zero prescribed marker speed and no forcing (after setting the amplitude to zero).
- [ ] Check the marker target velocity against the derivative of its prescribed displacement.
- [ ] Test interpolation on a constant velocity field and verify that it returns the same constant.
- [ ] Compare the integrated Eulerian spread force with the integrated Lagrangian marker force.
- [ ] Record the pressure Poisson residual as well as divergence before/after projection; do not rely on divergence alone.
- [ ] Increase pressure iterations and verify whether the projection diagnostics converge.
- [ ] Halve the time step and compare velocity, marker slip, and divergence.
- [ ] Vary the penalty gain and check for instability or excessive forcing.
- [ ] Refine the grid and marker spacing together; compare the same diagnostic quantities.
- [ ] Confirm all fields remain finite and that no unexplained energy growth occurs.
- [ ] Document the limitation: periodic tangential translation demonstrates the coupling mechanism but does not yet model chamber-volume change or pumping.
