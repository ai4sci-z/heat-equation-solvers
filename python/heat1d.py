"""
1D Heat equation: u_t = alpha * u_xx
Domain: x in [0,1], Dirichlet BC u=0, IC u=sin(pi*x)
Exact solution: exp(-alpha*pi^2*t) * sin(pi*x)
Method: explicit FTCS finite difference
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# ── parameters ──────────────────────────────────────────────────────────────
alpha = 1.0
N     = 50
T     = 0.1
x     = np.linspace(0, 1, N + 2)
dx    = x[1] - x[0]
r     = 0.4                        # CFL number: r = alpha*dt/dx^2 <= 0.5
dt    = r * dx**2 / alpha
steps = int(T / dt) + 1

# ── initial condition ────────────────────────────────────────────────────────
u = np.sin(np.pi * x)
u[0] = u[-1] = 0.0                 # enforce BC

# ── time march ───────────────────────────────────────────────────────────────
save_times = {0: None, int(steps * 0.5): None, steps - 1: None}
u_save = {}

for n in range(steps):
    if n in save_times:
        u_save[n * dt] = u.copy()
    u_new        = u.copy()
    u_new[1:-1]  = u[1:-1] + r * (u[2:] - 2*u[1:-1] + u[:-2])
    u_new[0]     = u_new[-1] = 0.0
    u             = u_new

# ── exact solution ────────────────────────────────────────────────────────────
def exact(t):
    return np.exp(-alpha * np.pi**2 * t) * np.sin(np.pi * x)

# ── L2 error at final time ────────────────────────────────────────────────────
t_final   = (steps - 1) * dt
u_ex      = exact(t_final)
l2_err    = np.sqrt(dx * np.sum((u - u_ex)**2))
print(f"1D Heat  N={N}  T={T:.2f}  L2 error = {l2_err:.3e}")

# ── plot ──────────────────────────────────────────────────────────────────────
os.makedirs("../results", exist_ok=True)
fig, axes = plt.subplots(1, 2, figsize=(11, 4))

ax = axes[0]
colors = plt.cm.plasma(np.linspace(0.2, 0.9, len(u_save)))
for (t_val, u_val), color in zip(u_save.items(), colors):
    ax.plot(x, u_val, color=color, lw=2, label=f't={t_val:.3f}')
    ax.plot(x, exact(t_val), '--', color=color, lw=1, alpha=0.7)
ax.set_xlabel('x'); ax.set_ylabel('u(x,t)')
ax.set_title('1D Heat Equation (solid=numerical, dashed=exact)')
ax.legend(fontsize=8); ax.grid(alpha=0.3)

ax = axes[1]
ax.plot(x, u - u_ex, 'r-', lw=1.5)
ax.set_xlabel('x'); ax.set_ylabel('error')
ax.set_title(f'Error at t={t_final:.3f}  L2={l2_err:.2e}')
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('../results/python_heat1d.png', dpi=150)
plt.show()
print("Saved ../results/python_heat1d.png")
