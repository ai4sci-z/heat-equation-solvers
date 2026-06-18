"""
2D Heat equation: u_t = alpha * (u_xx + u_yy)
Domain: [0,1]^2, Dirichlet BC u=0
IC: sin(pi*x)*sin(pi*y)
Exact: exp(-2*alpha*pi^2*t) * sin(pi*x)*sin(pi*y)
Method: explicit FTCS
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# ── parameters ───────────────────────────────────────────────────────────────
alpha = 1.0
N     = 50
T     = 0.1
x     = np.linspace(0, 1, N + 2)
dx    = x[1] - x[0]
r     = 0.2                          # r = alpha*dt/dx^2 <= 0.25 for 2D
dt    = r * dx**2 / alpha
steps = int(T / dt) + 1

X, Y = np.meshgrid(x, x, indexing='ij')

# ── initial condition ─────────────────────────────────────────────────────────
u = np.sin(np.pi * X) * np.sin(np.pi * Y)
u[0, :]  = u[-1, :] = 0.0
u[:, 0]  = u[:, -1] = 0.0

# ── time march ────────────────────────────────────────────────────────────────
snap_steps = [0, steps // 2, steps - 1]
snaps = {}

for n in range(steps):
    if n in snap_steps:
        snaps[n * dt] = u.copy()
    lap = (u[2:, 1:-1] - 2*u[1:-1, 1:-1] + u[:-2, 1:-1]) / dx**2 \
        + (u[1:-1, 2:] - 2*u[1:-1, 1:-1] + u[1:-1, :-2]) / dx**2
    u_new = u.copy()
    u_new[1:-1, 1:-1] = u[1:-1, 1:-1] + alpha * dt * lap
    u = u_new
    u[0, :] = u[-1, :] = u[:, 0] = u[:, -1] = 0.0

# ── exact & error ─────────────────────────────────────────────────────────────
def exact2d(t):
    return np.exp(-2 * alpha * np.pi**2 * t) * np.sin(np.pi * X) * np.sin(np.pi * Y)

t_final = (steps - 1) * dt
u_ex    = exact2d(t_final)
l2_err  = np.sqrt(dx**2 * np.sum((u - u_ex)**2))
print(f"2D Heat  N={N}  T={T:.2f}  L2 error = {l2_err:.3e}")

# ── plot ──────────────────────────────────────────────────────────────────────
os.makedirs("../results", exist_ok=True)
fig, axes = plt.subplots(2, 3, figsize=(14, 8))
times_list = sorted(snaps.keys())

for col, t_val in enumerate(times_list):
    u_s  = snaps[t_val]
    u_ex_s = exact2d(t_val)
    err  = u_s - u_ex_s

    vmax = u_s.max()
    axes[0, col].contourf(X, Y, u_s, levels=20, cmap='hot')
    axes[0, col].set_title(f'Numerical  t={t_val:.3f}')
    axes[0, col].set_xlabel('x'); axes[0, col].set_ylabel('y')
    im = axes[1, col].contourf(X, Y, err, levels=20, cmap='RdBu_r')
    axes[1, col].set_title(f'Error  t={t_val:.3f}')
    axes[1, col].set_xlabel('x'); axes[1, col].set_ylabel('y')
    plt.colorbar(im, ax=axes[1, col])

plt.suptitle(f'2D Heat Equation (FTCS)  N={N}  L2={l2_err:.2e}', fontsize=12)
plt.tight_layout()
plt.savefig('../results/python_heat2d.png', dpi=150)
plt.show()
print("Saved ../results/python_heat2d.png")
