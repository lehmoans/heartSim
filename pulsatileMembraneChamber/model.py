"""
Pulsatile Membrane Chamber — notebook code consolidated into one script.

This file preserves the code and code comments from the notebook in execution order.
The early cells are component tests; the later cells run the oscillating-lid prototype
and produce spatial velocity/pressure plots plus numerical diagnostics.

Dependencies: numpy, matplotlib

Important: the oscillating-lid solver is an educational prototype, not a validated
internal-membrane or production-quality immersed-boundary solver.
"""

import numpy as np
import matplotlib.pyplot as plt

# Illustrative, dimensionally consistent parameters.
y0 = 0.50       # mean membrane position [m]
amplitude = 0.05  # displacement amplitude [m]
frequency = 1.0   # oscillation frequency [Hz]

period = 1.0 / frequency
t = np.linspace(0.0, 2.0 * period, 500)

y_membrane = y0 + amplitude * np.sin(2.0 * np.pi * frequency * t)
v_membrane = 2.0 * np.pi * frequency * amplitude * np.cos(
    2.0 * np.pi * frequency * t
)

fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
axes[0].plot(t, y_membrane)
axes[0].set_ylabel("Membrane position y [m]")
axes[0].grid(True)
axes[1].plot(t, v_membrane)
axes[1].set_xlabel("Time [s]")
axes[1].set_ylabel("Membrane velocity [m/s]")
axes[1].grid(True)
plt.tight_layout()
plt.show()

# Section 3 — Spectral pressure-projection test
# Spectral projection test on a periodic square domain.
# This cell tests pressure correction only; it is not a full Navier–Stokes time step.

nx = ny = 64 #nodes
Lx = Ly = 1.0
dx, dy = Lx / nx, Ly / ny
dt = 0.001
rho = 1000.0

x = np.arange(nx) * dx #x points
y = np.arange(ny) * dy #y points
Xg, Yg = np.meshgrid(x, y, indexing="xy") #mesh

# Arbitrary smooth periodic predictor field, not necessarily divergence-free.
u_star = np.sin(2 * np.pi * Xg / Lx) + 0.25 * np.cos(2 * np.pi * Yg / Ly)
v_star = 0.4 * np.cos(2 * np.pi * Xg / Lx) + np.sin(2 * np.pi * Yg / Ly)

kx = 2 * np.pi * np.fft.fftfreq(nx, d=dx)
ky = 2 * np.pi * np.fft.fftfreq(ny, d=dy)
KX, KY = np.meshgrid(kx, ky, indexing="xy")
K2 = KX**2 + KY**2

u_hat = np.fft.fft2(u_star)
v_hat = np.fft.fft2(v_star)

# Fourier representation of divergence: i*kx*u_hat + i*ky*v_hat.
div_hat = 1j * KX * u_hat + 1j * KY * v_hat

# Solve -|k|^2 p_hat = (rho/dt) div_hat.
p_hat = np.zeros_like(div_hat, dtype=complex)
nonzero = K2 > 0
p_hat[nonzero] = -(rho / dt) * div_hat[nonzero] / K2[nonzero]
p_hat[0, 0] = 0.0  # pressure's arbitrary constant is set to zero

# Correct velocity in Fourier space.
u_new_hat = u_hat - (dt / rho) * (1j * KX * p_hat)
v_new_hat = v_hat - (dt / rho) * (1j * KY * p_hat)
u_new = np.fft.ifft2(u_new_hat).real
v_new = np.fft.ifft2(v_new_hat).real

div_before = np.fft.ifft2(div_hat).real
div_after_hat = 1j * KX * u_new_hat + 1j * KY * v_new_hat
div_after = np.fft.ifft2(div_after_hat).real

print(f"Maximum |divergence| before projection: {np.max(np.abs(div_before)):.3e}")
print(f"Maximum |divergence| after projection:  {np.max(np.abs(div_after)):.3e}")

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
im0 = axes[0].imshow(div_before, origin="lower", extent=[0, Lx, 0, Ly])
axes[0].set_title("Divergence before projection")
plt.colorbar(im0, ax=axes[0], shrink=0.8)
im1 = axes[1].imshow(div_after, origin="lower", extent=[0, Lx, 0, Ly])
axes[1].set_title("Divergence after projection")
plt.colorbar(im1, ax=axes[1], shrink=0.8)
for ax in axes:
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
plt.tight_layout()
plt.show()

# Section 4 — Regularised Peskin delta kernel
def peskin_phi(r):
    """Four-point Peskin 1D kernel phi(r), with support |r| <= 2."""
    r = np.abs(np.asarray(r, dtype=float))
    out = np.zeros_like(r)

    m1 = r <= 1.0
    out[m1] = (
        3.0 - 2.0 * r[m1]
        + np.sqrt(np.maximum(0.0, 1.0 + 4.0 * r[m1] - 4.0 * r[m1]**2))
    ) / 8.0

    m2 = (r > 1.0) & (r <= 2.0)
    out[m2] = (
        5.0 - 2.0 * r[m2]
        - np.sqrt(np.maximum(0.0, -7.0 + 12.0 * r[m2] - 4.0 * r[m2]**2))
    ) / 8.0

    return out

r = np.linspace(-2.5, 2.5, 500)
plt.figure(figsize=(7, 3))
plt.plot(r, peskin_phi(r))
plt.xlabel("r")
plt.ylabel("phi(r)")
plt.title("Four-point Peskin regularised delta kernel")
plt.grid(True)
plt.show()

# Section 5 — Force spreading and conservation check
# Demonstration: spread vertical force density from markers on a periodic horizontal line.
# The marker force values below are illustrative, not physiological.

nx = ny = 64
Lx = Ly = 1.0
dx, dy = Lx / nx, Ly / ny
x = np.arange(nx) * dx
y = np.arange(ny) * dy
Xgrid, Ygrid = np.meshgrid(x, y, indexing="xy")

n_markers = 32
Xm = np.linspace(0.0, Lx, n_markers, endpoint=False)
Ym = np.full(n_markers, 0.5)
ds = Lx / n_markers
f_marker_y = np.ones(n_markers) * 2.0  # N/m, illustrative force density

Fgrid_y = np.zeros((ny, nx))

def periodic_distance(a, b, length):
    """Signed shortest periodic distance from grid coordinate a to marker b."""
    d = a - b
    return (d + 0.5 * length) % length - 0.5 * length

for xm, ym, fm in zip(Xm, Ym, f_marker_y):
    rx = periodic_distance(Xgrid, xm, Lx) / dx
    ry = periodic_distance(Ygrid, ym, Ly) / dy
    delta_x = peskin_phi(rx) / dx
    delta_y = peskin_phi(ry) / dy
    Fgrid_y += fm * delta_x * delta_y * ds

force_from_markers = np.sum(f_marker_y * ds)
force_on_grid = np.sum(Fgrid_y) * dx * dy

print(f"Integrated marker force: {force_from_markers:.8f} N")
print(f"Integrated grid force:   {force_on_grid:.8f} N")
print(f"Absolute difference:     {abs(force_from_markers-force_on_grid):.3e} N")

plt.figure(figsize=(7, 4))
plt.imshow(Fgrid_y, origin="lower", extent=[0, Lx, 0, Ly], aspect="equal")
plt.colorbar(label="Spread vertical force density [N/m²]")
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.title("Boundary force spread onto the fixed grid")
plt.show()

# Section 6 — Interpolate velocity to membrane markers
# Interpolate the projected test velocity to marker points using the same kernel.
U_marker_x = np.zeros(n_markers)
U_marker_y = np.zeros(n_markers)

for k, (xm, ym) in enumerate(zip(Xm, Ym)):
    rx = periodic_distance(Xgrid, xm, Lx) / dx
    ry = periodic_distance(Ygrid, ym, Ly) / dy
    weights = (peskin_phi(rx) / dx) * (peskin_phi(ry) / dy) * dx * dy
    U_marker_x[k] = np.sum(u_new * weights)
    U_marker_y[k] = np.sum(v_new * weights)

print("Interpolated marker velocities (first five markers):")
for k in range(min(5, n_markers)):
    print(f"marker {k:2d}: Ux={U_marker_x[k]: .5e} m/s, Uy={U_marker_y[k]: .5e} m/s")

# Section 7 — Time-step estimates
# Example time-step estimates. Replace these values with the parameters for your case.
U_max = 0.5       # maximum speed [m/s], assumed for estimate
h = 1.0 / 64      # grid spacing [m] for a unit-length domain
nu = 1.0e-6       # kinematic viscosity [m^2/s]
C_adv = 0.5

dt_adv = C_adv * h / U_max
dt_visc = h**2 / (4.0 * nu)
dt_estimate = min(dt_adv, dt_visc)

print(f"Advection estimate: {dt_adv:.4e} s")
print(f"Explicit diffusion estimate: {dt_visc:.4e} s")
print(f"Minimum of these two estimates: {dt_estimate:.4e} s")
print("This is not yet a final time step: check membrane coupling and accuracy too.")

# Section 10 — End-to-end oscillating-lid cavity prototype
# End-to-end oscillating-lid cavity prototype.
# Dependencies: NumPy and Matplotlib only.
# All values are illustrative and dimensionless / consistently scaled.

import numpy as np
import matplotlib.pyplot as plt

# Grid and fluid parameters
nx = ny = 41
Lx = Ly = 1.0
dx = Lx / (nx - 1)
dy = Ly / (ny - 1)
rho = 1.0
nu = 0.05
dt_flow = 0.001
n_steps = 600
pressure_iterations = 60

# Oscillating top-wall speed: U_lid(t) = U0 sin(2 pi f t)
U0 = 0.1
f_lid = 1.0
period_lid = 1.0 / f_lid

x_flow = np.linspace(0.0, Lx, nx)
y_flow = np.linspace(0.0, Ly, ny)
X_flow, Y_flow = np.meshgrid(x_flow, y_flow)

u_flow = np.zeros((ny, nx))
v_flow = np.zeros((ny, nx))
p_flow = np.zeros((ny, nx))
times = []
max_speeds = []
divergence_rms = []
lid_speeds = []

def apply_velocity_bc(u, v, lid_speed):
    # No-slip stationary side and bottom walls; prescribed oscillating top wall.
    u[:, 0] = 0.0
    u[:, -1] = 0.0
    u[0, :] = 0.0
    u[-1, :] = lid_speed
    v[:, 0] = 0.0
    v[:, -1] = 0.0
    v[0, :] = 0.0
    v[-1, :] = 0.0

def pressure_poisson(p, b, iterations):
    # Jacobi iteration for a teaching pressure-Poisson solve.
    pn = np.empty_like(p)
    for _ in range(iterations):
        pn[:] = p
        p[1:-1, 1:-1] = (
            (pn[1:-1, 2:] + pn[1:-1, :-2]) * dy**2
            + (pn[2:, 1:-1] + pn[:-2, 1:-1]) * dx**2
            - b[1:-1, 1:-1] * dx**2 * dy**2
        ) / (2.0 * (dx**2 + dy**2))

        # Approximate homogeneous Neumann pressure boundaries, with a gauge fix.
        p[:, 0] = p[:, 1]
        p[:, -1] = p[:, -2]
        p[-1, :] = p[-2, :]
        p[0, :] = 0.0
    return p

for step in range(n_steps):
    t_now = step * dt_flow
    lid_speed = U0 * np.sin(2.0 * np.pi * f_lid * t_now)

    # Enforce current wall values before computing derivatives.
    apply_velocity_bc(u_flow, v_flow, lid_speed)
    un = u_flow.copy()
    vn = v_flow.copy()

    dudx = (un[1:-1, 2:] - un[1:-1, :-2]) / (2.0 * dx)
    dudy = (un[2:, 1:-1] - un[:-2, 1:-1]) / (2.0 * dy)
    dvdx = (vn[1:-1, 2:] - vn[1:-1, :-2]) / (2.0 * dx)
    dvdy = (vn[2:, 1:-1] - vn[:-2, 1:-1]) / (2.0 * dy)

    b = np.zeros_like(p_flow)
    b[1:-1, 1:-1] = rho * (
        (dudx + dvdy) / dt_flow
        - dudx**2
        - 2.0 * dudy * dvdx
        - dvdy**2
    )
    p_flow = pressure_poisson(p_flow, b, pressure_iterations)

    # Explicit advection + diffusion + pressure gradient.
    u_flow[1:-1, 1:-1] = (
        un[1:-1, 1:-1]
        - un[1:-1, 1:-1] * dt_flow * dudx
        - vn[1:-1, 1:-1] * dt_flow * dudy
        - dt_flow / (2.0 * rho * dx) * (p_flow[1:-1, 2:] - p_flow[1:-1, :-2])
        + nu * dt_flow * (
            (un[1:-1, 2:] - 2.0 * un[1:-1, 1:-1] + un[1:-1, :-2]) / dx**2
            + (un[2:, 1:-1] - 2.0 * un[1:-1, 1:-1] + un[:-2, 1:-1]) / dy**2
        )
    )
    v_flow[1:-1, 1:-1] = (
        vn[1:-1, 1:-1]
        - un[1:-1, 1:-1] * dt_flow * dvdx
        - vn[1:-1, 1:-1] * dt_flow * dvdy
        - dt_flow / (2.0 * rho * dy) * (p_flow[2:, 1:-1] - p_flow[:-2, 1:-1])
        + nu * dt_flow * (
            (vn[1:-1, 2:] - 2.0 * vn[1:-1, 1:-1] + vn[1:-1, :-2]) / dx**2
            + (vn[2:, 1:-1] - 2.0 * vn[1:-1, 1:-1] + vn[:-2, 1:-1]) / dy**2
        )
    )

    apply_velocity_bc(u_flow, v_flow, lid_speed)

    # Save diagnostics at each step. Boundary derivatives are excluded from RMS.
    div_flow = (
        (u_flow[1:-1, 2:] - u_flow[1:-1, :-2]) / (2.0 * dx)
        + (v_flow[2:, 1:-1] - v_flow[:-2, 1:-1]) / (2.0 * dy)
    )
    times.append(t_now)
    lid_speeds.append(lid_speed)
    max_speeds.append(np.max(np.sqrt(u_flow**2 + v_flow**2)))
    divergence_rms.append(np.sqrt(np.mean(div_flow**2)))

print(f"Simulated duration: {n_steps * dt_flow:.3f} s")
print(f"Maximum speed at final step: {max_speeds[-1]:.4e}")
print(f"Interior RMS divergence at final step: {divergence_rms[-1]:.4e}")

# Section 11 — Flow-field plots and diagnostics
# Plot 1: final velocity magnitude and velocity vectors.
speed_flow = np.sqrt(u_flow**2 + v_flow**2)
fig, ax = plt.subplots(figsize=(7, 5))
im = ax.contourf(X_flow, Y_flow, speed_flow, levels=20)
fig.colorbar(im, ax=ax, label="Speed magnitude")
skip = 2
ax.quiver(
    X_flow[::skip, ::skip], Y_flow[::skip, ::skip],
    u_flow[::skip, ::skip], v_flow[::skip, ::skip],
    color="white", scale=1.5
)
ax.set_title(f"Oscillating-lid cavity: speed at t={times[-1]:.3f}")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_aspect("equal")
plt.tight_layout()
plt.show()

# Plot 2: streamlines and pressure.
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
axes[0].streamplot(x_flow, y_flow, u_flow, v_flow, density=1.2)
axes[0].set_title("Final velocity streamlines")
axes[0].set_xlabel("x")
axes[0].set_ylabel("y")
axes[0].set_aspect("equal")

im2 = axes[1].contourf(X_flow, Y_flow, p_flow, levels=20)
fig.colorbar(im2, ax=axes[1], label="Pressure (relative)")
axes[1].set_title("Final relative pressure")
axes[1].set_xlabel("x")
axes[1].set_ylabel("y")
axes[1].set_aspect("equal")
plt.tight_layout()
plt.show()

# Plot 3: spatial velocity components u and v across the domain.
# These distinguish horizontal and vertical flow, which speed magnitude alone hides.
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
im_u = axes[0].contourf(X_flow, Y_flow, u_flow, levels=21)
fig.colorbar(im_u, ax=axes[0], label="Horizontal velocity u")
axes[0].set_title("Horizontal velocity field u(x, y)")
axes[0].set_xlabel("x")
axes[0].set_ylabel("y")
axes[0].set_aspect("equal")

im_v = axes[1].contourf(X_flow, Y_flow, v_flow, levels=21)
fig.colorbar(im_v, ax=axes[1], label="Vertical velocity v")
axes[1].set_title("Vertical velocity field v(x, y)")
axes[1].set_xlabel("x")
axes[1].set_ylabel("y")
axes[1].set_aspect("equal")
plt.show()

# Plot 4: pressure field with the velocity field overlaid.
# Pressure is relative (gauge-dependent); arrows show the local flow direction.
fig, ax = plt.subplots(figsize=(7, 5))
im_p = ax.contourf(X_flow, Y_flow, p_flow, levels=21)
fig.colorbar(im_p, ax=ax, label="Relative pressure")
skip = 2
ax.quiver(
    X_flow[::skip, ::skip], Y_flow[::skip, ::skip],
    u_flow[::skip, ::skip], v_flow[::skip, ::skip],
    color="white", scale=1.5
)
ax.set_title("Pressure field with velocity vectors")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_aspect("equal")
plt.tight_layout()
plt.show()

# Plot 5: velocity vectors coloured by speed magnitude.
fig, ax = plt.subplots(figsize=(7, 5))
im_q = ax.quiver(
    X_flow[::skip, ::skip], Y_flow[::skip, ::skip],
    u_flow[::skip, ::skip], v_flow[::skip, ::skip],
    speed_flow[::skip, ::skip], scale=1.5
)
fig.colorbar(im_q, ax=ax, label="Speed magnitude")
ax.set_title("Velocity vector field coloured by speed")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_aspect("equal")
plt.tight_layout()
plt.show()

# Plot 6: imposed wall speed and solver diagnostics over time.
fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
axes[0].plot(times, lid_speeds)
axes[0].set_ylabel("Top-wall speed")
axes[0].set_title("Prescribed pulsation")
axes[0].grid(True)

axes[1].plot(times, divergence_rms, label="Interior RMS divergence")
axes[1].set_xlabel("Time")
axes[1].set_ylabel("RMS divergence")
axes[1].set_yscale("log")
axes[1].grid(True)
axes[1].legend()
plt.tight_layout()
plt.show()

# Plot 7: maximum speed through time.
plt.figure(figsize=(8, 3.5))
plt.plot(times, max_speeds)
plt.xlabel("Time")
plt.ylabel("Maximum speed")
plt.title("Flow response over the simulation")
plt.grid(True)
plt.tight_layout()
plt.show()


