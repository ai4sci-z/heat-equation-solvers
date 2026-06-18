"""
3D Heat equation: u_t = alpha*(u_xx + u_yy + u_zz)
Domain: [0,1]^3, Dirichlet BC u=0
IC: sin(pi*x)*sin(pi*y)*sin(pi*z)
Exact: exp(-3*alpha*pi^2*t) * sin(pi*x)*sin(pi*y)*sin(pi*z)
Method: explicit FTCS
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# ── parameters ───────────────────────────────────────────────────────────────
alpha = 1.0
N     = 20                           # smaller grid for 3D memory
T     = 0.05
x     = np.linspace(0, 1, N + 2)
dx    = x[1] - x[0]
r     = 0.13                         # r = alpha*dt/dx^2 <= 1/6 for 3D
dt    = r * dx**2 / alpha
steps = int(T / dt) + 1

X, Y, Z = np.meshgrid(x, x, x, indexing='ij')

# ── initial condition ─────────────────────────────────────────────────────────
u = np.sin(np.pi * X) * np.sin(np.pi * Y) * np.sin(np.pi * Z)
u[0,:,:]=u[-1,:,:]=u[:,0,:]=u[:,-1,:]=u[:,:,0]=u[:,:,-1]=0.0

def apply_bc(v):
    v[0,:,:]=v[-1,:,:]=v[:,0,:]=v[:,-1,:]=v[:,:,0]=v[:,:,-1]=0.0

# ── time march ────────────────────────────────────────────────────────────────
snap_steps = [0, steps // 2, steps - 1]
snaps = {}

for n in range(steps):
    if n in snap_steps:
        snaps[n * dt] = u.copy()
    lap = ((u[2:,1:-1,1:-1] - 2*u[1:-1,1:-1,1:-1] + u[:-2,1:-1,1:-1])
         + (u[1:-1,2:,1:-1] - 2*u[1:-1,1:-1,1:-1] + u[1:-1,:-2,1:-1])
         + (u[1:-1,1:-1,2:] - 2*u[1:-1,1:-1,1:-1] + u[1:-1,1:-1,:-2])) / dx**2
    u_new = u.copy()
    u_new[1:-1,1:-1,1:-1] = u[1:-1,1:-1,1:-1] + alpha * dt * lap
    apply_bc(u_new)
    u = u_new

# ── exact & error ─────────────────────────────────────────────────────────────
def exact3d(t):
    return np.exp(-3*alpha*np.pi**2*t) * np.sin(np.pi*X)*np.sin(np.pi*Y)*np.sin(np.pi*Z)

t_final = (steps - 1) * dt
u_ex    = exact3d(t_final)
l2_err  = np.sqrt(dx**3 * np.sum((u - u_ex)**2))
print(f"3D Heat  N={N}  T={T:.3f}  L2 error = {l2_err:.3e}")

# ── plot: z=0.5 midplane slices ────────────────────────────────────────────────
os.makedirs("../results", exist_ok=True)
mid = N // 2 + 1  # index of z=0.5

fig, axes = plt.subplots(2, 3, figsize=(14, 8))
times_list = sorted(snaps.keys())

for col, t_val in enumerate(times_list):
    u_s     = snaps[t_val][:, :, mid]
    u_ex_s  = exact3d(t_val)[:, :, mid]
    err     = u_s - u_ex_s

    axes[0, col].contourf(X[:,:,mid], Y[:,:,mid], u_s, levels=20, cmap='hot')
    axes[0, col].set_title(f'Numerical z=0.5  t={t_val:.3f}')
    axes[0, col].set_xlabel('x'); axes[0, col].set_ylabel('y')
    im = axes[1, col].contourf(X[:,:,mid], Y[:,:,mid], err, levels=20, cmap='RdBu_r')
    axes[1, col].set_title(f'Error z=0.5  t={t_val:.3f}')
    axes[1, col].set_xlabel('x'); axes[1, col].set_ylabel('y')
    plt.colorbar(im, ax=axes[1, col])

plt.suptitle(f'3D Heat Equation (FTCS) — midplane z=0.5  N={N}  L2={l2_err:.2e}', fontsize=11)
plt.tight_layout()
plt.savefig('../results/python_heat3d.png', dpi=150)
plt.show()
print("Saved ../results/python_heat3d.png")
