# HeartSim — Simulation Verification Checklist

## Purpose

Use this file to decide whether a simulation result is trustworthy before interpreting it as engineering evidence. A plot that looks plausible, a solver that runs, or a small residual alone does **not** establish accuracy.

Work through the checks in order. Record the method, result, and evidence for each completed check. Do not mark a check complete just because AI-generated code says it passed.

## Status key

- [ ] Not started
- [x] Completed and supported by recorded evidence
- [!] Failed, blocked, or needs investigation

## 1. Define the model before testing it

- [ ] Write down the physical question the model is intended to answer.
- [ ] State what is modelled and what is deliberately excluded.
- [ ] List assumptions, including fluid properties, geometry, dimensionality, and membrane motion.
- [ ] Define all variables and units; check dimensional consistency of governing equations.
- [ ] Specify initial conditions and every boundary condition.
- [ ] Identify the outputs that will be used to make decisions (for example pressure difference, flow rate, velocity, displacement, or power).
- [ ] Separate prescribed inputs from quantities the solver is expected to predict.

**Evidence to save:** model description, governing equations, parameter table, and a sketch of the domain and boundaries.

## 2. Verify the implementation

- [ ] Check that the implemented equations match the documented equations and sign conventions.
- [ ] Check array dimensions, coordinate indexing, units, and boundary-condition application.
- [ ] Check that pressure and velocity are not accidentally mixed between time steps or overwritten by later notebook/code sections.
- [ ] Run the smallest case possible and confirm all computed values are finite (no unexpected NaN or infinity).
- [ ] Add deterministic tests for key numerical functions.
- [ ] Re-run tests from a clean Python process, not only in an already-used notebook session.

## 3. Conservation and physical consistency

- [ ] Measure mass/volume conservation over the domain.
- [ ] Report a normalized conservation error and define an acceptable tolerance before running the test.
- [ ] Check the sign and magnitude of pressure gradients, velocities, forces, and fluxes against physical expectations.
- [ ] Check that units are consistent in force, pressure, flow, energy, and power calculations.
- [ ] Investigate any unphysical growth, unexplained drift, or asymmetric result in a symmetric setup.

## 4. Verification against known solutions

- [ ] Select a simple case with an analytical solution or trusted benchmark.
- [ ] Match the benchmark's geometry, assumptions, boundary conditions, and parameter values.
- [ ] Compare quantitative outputs, not just plots.
- [ ] Calculate absolute and relative errors; handle near-zero reference values carefully.
- [ ] Record whether the test passes the tolerance set in advance.
- [ ] If no suitable benchmark exists, label the result as unverified rather than accurate.

## 5. Numerical convergence and stability

- [ ] Run at least three systematically refined spatial resolutions where feasible.
- [ ] Reduce the time step independently of spatial resolution.
- [ ] Compare selected quantities of interest across refinements.
- [ ] Check whether errors decrease at the expected rate for the method and test case.
- [ ] Test the solver's stability limits and report the chosen time step and relevant dimensionless numbers.
- [ ] Distinguish a converged iterative solve from a grid/time-step-converged physical prediction.

## 6. Sensitivity and uncertainty

- [ ] Identify parameters that are uncertain or weakly known.
- [ ] Vary important parameters over justified ranges and report their effect on the outputs.
- [ ] Separate numerical error, parameter uncertainty, and model-form limitations.
- [ ] Avoid reporting more precision than the assumptions and evidence support.

## 7. Reproducibility and reporting

- [ ] Record Python and dependency versions, run command, configuration, and random seeds (if any).
- [ ] Save the exact parameter values and solver settings used to generate each result.
- [ ] Generate plots with labelled axes, units, legends, and clear time/position information.
- [ ] Save numerical diagnostics alongside figures.
- [ ] Make a clean, repeatable run possible from the repository instructions.
- [ ] Record failed tests and limitations; do not silently delete or bypass failures.

## 8. Before making an engineering claim

- [ ] State which verification checks passed and which remain incomplete.
- [ ] State whether results are benchmarked, numerically converged, and/or validated against physical measurements.
- [ ] Do not call a result **validated against reality** unless suitable experimental or observational data have been compared.
- [ ] State the intended range of validity and known limitations.

## First task

Start with **one small, measurable test case**. Write down the expected result before running the code, compare the computed value against it, and document the error. Do not begin by adding more model complexity.

## Test log

| Date | Test / case | Expected result or tolerance | Measured result / error | Pass? | Evidence / notes |
|---|---|---|---|---|---|
| | | | | | |
