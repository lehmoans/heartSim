"""
Stage 2: immersed-boundary coupling demonstration.

A prescribed, periodically translating line of Lagrangian markers is coupled to a
2D incompressible periodic fluid grid. A regularised Peskin kernel interpolates
fluid velocity to markers and spreads marker forcing back to the Eulerian grid.

This is a numerical-methods demonstration, NOT a pulsatile chamber or physiological
heart model. The periodic line translates tangentially to avoid changing the
enclosed volume; a closed deforming membrane and chamber walls are later milestones.

Dependencies: numpy, matplotlib
Run: python stage2_immersed_boundary_demo.py
"""
import numpy as np
import matplotlib.pyplot as plt


def peskin_phi(r):
    """Four-point Peskin kernel, support |r| <= 2."""
    r = np.abs(np.asarray(r, dtype=float))
    out = np.zeros_like(r)
    m1 = r <= 1.0
    out[m1] = (
        3.0 - 2.0 * r[m1]
        + np.sqrt(np.maximum(0.0, 1.0 + 4.0 * r[m1] - 4.0 * r[m1] ** 2))
    ) / 8.0
    m2 = (r > 1.0) & (r <= 2.0)
    out[m2] = (
        5.0 - 2.0 * r[m2]
        - np.sqrt(np.maximum(0.0, -7.0 + 12.0 * r[m2] - 4.0 * r[m2] ** 2))
    ) / 8.0
    return out


def periodic_distance(a, b, length):
    """Shortest signed distance between periodic coordinates."""
    return (a - b + 0.5 * length) % length - 0.5 * length


def interpolate_marker_velocity(u, v, xm, ym, X, Y, dx, dy, Lx, Ly):
    """Interpolate Eulerian velocity to one Lagrangian marker."""
    rx = periodic_distance(X, xm, Lx) / dx
    ry = periodic_distance(Y, ym, Ly) / dy
    weights = peskin_phi(rx) * peskin_phi(ry)  # dimensionless kernel weights
    weights /= np.sum(weights)
    return np.sum(u * weights), np.sum(v * weights)


def spread_marker_forces(fx, fy, xm, ym, fmx, fmy, ds, X, Y, dx, dy, Lx, Ly):
    """Spread marker force densities to the periodic Eulerian grid."""
    for k in range(len(xm)):
        rx = periodic_distance(X, xm[k], Lx) / dx
        ry = periodic_distance(Y, ym[k], Ly) / dy
        delta = (peskin_phi(rx) / dx) * (peskin_phi(ry) / dy)
        fx += fmx[k] * delta * ds
        fy += fmy[k] * delta * ds
    return fx, fy


def pressure_poisson_periodic(rhs, dx, dy, iterations=200, initial=None):
    """Jacobi solve of the periodic five-point Poisson equation; pressure mean is zero."""
    p = np.zeros_like(rhs) if initial is None else initial.copy()
    p_next = np.zeros_like(rhs)
    dx2, dy2 = dx * dx, dy * dy
    denom = 2.0 * (dx2 + dy2)
    for _ in range(iterations):
        p_next[:] = (
            (np.roll(p, -1, axis=1) + np.roll(p, 1, axis=1)) * dy2
            + (np.roll(p, -1, axis=0) + np.roll(p, 1, axis=0)) * dx2
            - rhs * dx2 * dy2
        ) / denom
        p_next -= np.mean(p_next)  # remove the arbitrary pressure constant
        p, p_next = p_next, p
    return p


def forward_divergence(u, v, dx, dy):
    """Divergence operator compatible with the backward pressure gradient."""
    return (
        (np.roll(u, -1, axis=1) - u) / dx
        + (np.roll(v, -1, axis=0) - v) / dy
    )


def backward_gradient(p, dx, dy):
    """Backward pressure gradient, compatible with forward divergence."""
    gx = (p - np.roll(p, 1, axis=1)) / dx
    gy = (p - np.roll(p, 1, axis=0)) / dy
    return gx, gy


def laplacian(q, dx, dy):
    return (
        (np.roll(q, -1, axis=1) - 2.0 * q + np.roll(q, 1, axis=1)) / dx**2
        + (np.roll(q, -1, axis=0) - 2.0 * q + np.roll(q, 1, axis=0)) / dy**2
    )


def main():
    # Nondimensional demonstration parameters.
    nx = ny = 32
    Lx = Ly = 1.0
    dx, dy = Lx / nx, Ly / ny
    rho = 1.0
    nu = 0.01
    dt = 0.001
    n_steps = 200
    pressure_iterations = 200

    # Prescribed line motion: tangential translation, not chamber pumping.
    n_markers = 32
    ds = Lx / n_markers
    X0 = np.linspace(0.0, Lx, n_markers, endpoint=False)
    Y0 = np.full(n_markers, 0.5 * Ly)
    amplitude = 0.03
    frequency = 1.0
    penalty = 2.0  # illustrative feedback gain; requires stability/sensitivity study

    x = np.arange(nx) * dx
    y = np.arange(ny) * dy
    X, Y = np.meshgrid(x, y, indexing="xy")
    u = np.zeros((ny, nx))
    v = np.zeros((ny, nx))
    p = np.zeros((ny, nx))

    times, target_speeds, marker_speeds = [], [], []
    div_rms = []
    poisson_residual_rms = []
    final_fx = np.zeros_like(u)
    final_fy = np.zeros_like(v)

    for step in range(n_steps):
        t = step * dt
        displacement = amplitude * np.sin(2.0 * np.pi * frequency * t)
        target_ux = 2.0 * np.pi * frequency * amplitude * np.cos(
            2.0 * np.pi * frequency * t
        )
        xm = (X0 + displacement) % Lx
        ym = Y0.copy()

        # Interpolate current fluid velocity to the moving markers.
        um = np.zeros(n_markers)
        vm = np.zeros(n_markers)
        for k in range(n_markers):
            um[k], vm[k] = interpolate_marker_velocity(
                u, v, xm[k], ym[k], X, Y, dx, dy, Lx, Ly
            )

        # Penalty feedback: marker forcing attempts to reduce no-slip mismatch.
        fmx = penalty * (target_ux - um)
        fmy = penalty * (0.0 - vm)
        fx = np.zeros_like(u)
        fy = np.zeros_like(v)
        fx, fy = spread_marker_forces(
            fx, fy, xm, ym, fmx, fmy, ds, X, Y, dx, dy, Lx, Ly
        )

        # Explicit incompressible Navier-Stokes predictor: advection + viscosity + IB force.
        dudx = (np.roll(u, -1, axis=1) - np.roll(u, 1, axis=1)) / (2.0 * dx)
        dudy = (np.roll(u, -1, axis=0) - np.roll(u, 1, axis=0)) / (2.0 * dy)
        dvdx = (np.roll(v, -1, axis=1) - np.roll(v, 1, axis=1)) / (2.0 * dx)
        dvdy = (np.roll(v, -1, axis=0) - np.roll(v, 1, axis=0)) / (2.0 * dy)

        u_star = u + dt * (
            -u * dudx - v * dudy + nu * laplacian(u, dx, dy) + fx / rho
        )
        v_star = v + dt * (
            -u * dvdx - v * dvdy + nu * laplacian(v, dx, dy) + fy / rho
        )

        # Pressure projection. The discrete divergence and gradient are paired.
        rhs = (rho / dt) * forward_divergence(u_star, v_star, dx, dy)
        p = pressure_poisson_periodic(rhs, dx, dy, pressure_iterations, initial=p)
        poisson_residual = laplacian(p, dx, dy) - rhs
        poisson_residual_rms.append(float(np.sqrt(np.mean(poisson_residual**2))))
        grad_px, grad_py = backward_gradient(p, dx, dy)
        u = u_star - (dt / rho) * grad_px
        v = v_star - (dt / rho) * grad_py

        times.append(t)
        target_speeds.append(target_ux)
        marker_speeds.append(float(np.mean(um)))
        div_rms.append(float(np.sqrt(np.mean(forward_divergence(u, v, dx, dy) ** 2))))
        final_fx, final_fy = fx, fy

    speed = np.sqrt(u**2 + v**2)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    im0 = axes[0].contourf(X, Y, speed, levels=18)
    fig.colorbar(im0, ax=axes[0], label="Speed")
    axes[0].quiver(X[::2, ::2], Y[::2, ::2], u[::2, ::2], v[::2, ::2])
    axes[0].set_title("Final fluid velocity")

    im1 = axes[1].contourf(X, Y, final_fx, levels=18)
    fig.colorbar(im1, ax=axes[1], label="Spread x-force density")
    axes[1].set_title("Immersed-boundary forcing")

    axes[2].plot(times, target_speeds, label="Prescribed marker speed")
    axes[2].plot(times, marker_speeds, label="Interpolated fluid speed")
    axes[2].set_title("No-slip coupling response")
    axes[2].set_xlabel("Time")
    axes[2].legend()
    for ax in axes:
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_aspect("equal") if ax is not axes[2] else None
    plt.tight_layout()
    plt.show()

    print(f"Grid: {nx} x {ny}; steps: {n_steps}; dt: {dt:g}")
    print(f"Final RMS discrete divergence: {div_rms[-1]:.4e}")
    print(f"Final RMS pressure-Poisson residual: {poisson_residual_rms[-1]:.4e}")
    print(f"Final maximum speed: {np.max(speed):.4e}")
    print("Prototype only: inspect time-step, pressure-iteration, and penalty-gain sensitivity.")
    print("Not yet a closed chamber, validated solver, or physiological pump model.")


if __name__ == "__main__":
    main()
