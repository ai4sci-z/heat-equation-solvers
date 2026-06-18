# 1D Heat equation: u_t = alpha*u_xx
# Domain: [0,1], BC: u=0, IC: sin(pi*x)
# Exact: exp(-alpha*pi^2*t)*sin(pi*x)
# Method: explicit FTCS

using LinearAlgebra, Printf
using Plots; gr()

# ── parameters ────────────────────────────────────────────────────────────────
const alpha = 1.0
const N     = 50
const T     = 0.1
x  = LinRange(0, 1, N+2)
dx = x[2] - x[1]
r  = 0.4
dt = r * dx^2 / alpha
steps = Int(floor(T / dt)) + 1

# ── initial condition ─────────────────────────────────────────────────────────
u     = sin.(π .* x)
u[1]  = u[end] = 0.0
u_new = copy(u)

# ── snapshots ─────────────────────────────────────────────────────────────────
snap_steps = Set([1, div(steps,2), steps])
snaps = Dict{Float64, Vector{Float64}}()

# ── time march ────────────────────────────────────────────────────────────────
for n in 1:steps
    n in snap_steps && (snaps[n*dt] = copy(u))
    @inbounds for i in 2:N+1
        u_new[i] = u[i] + r*(u[i+1] - 2u[i] + u[i-1])
    end
    u_new[1] = u_new[end] = 0.0
    u, u_new = u_new, u
end

# ── exact & error ─────────────────────────────────────────────────────────────
exact(t) = exp(-alpha * π^2 * t) .* sin.(π .* x)
t_final  = (steps-1) * dt
l2_err   = sqrt(dx * sum((u .- exact(t_final)).^2))
@printf("1D Heat  N=%d  T=%.2f  L2 error = %.3e\n", N, T, l2_err)

# ── plot ──────────────────────────────────────────────────────────────────────
mkpath("../results")
p1 = plot(title="1D Heat (solid=num, dashed=exact)", xlabel="x", ylabel="u",
          legend=:topright)
colors = [:blue, :red, :green]
t_keys = sort(collect(keys(snaps)))
for (k, t_val) in enumerate(t_keys)
    plot!(p1, x, snaps[t_val], lw=2, color=colors[k], label="num t=$(round(t_val,digits=3))")
    plot!(p1, x, exact(t_val), lw=1, ls=:dash, color=colors[k], label=false)
end

p2 = plot(x, u .- exact(t_final), lw=1.5, color=:red,
          title=@sprintf("Error t=%.3f  L2=%.2e", t_final, l2_err),
          xlabel="x", ylabel="error", legend=false)

plt = plot(p1, p2, layout=(1,2), size=(1000,400),
           plot_title="Julia — 1D Heat Equation (FTCS)  N=$N")
savefig(plt, "../results/julia_heat1d.png")
println("Saved ../results/julia_heat1d.png")
