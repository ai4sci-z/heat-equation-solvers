# 2D Heat equation: u_t = alpha*(u_xx + u_yy)
# Domain: [0,1]^2, Dirichlet BC, IC: sin(pi*x)*sin(pi*y)
# Exact: exp(-2*alpha*pi^2*t)*sin(pi*x)*sin(pi*y)
# Method: explicit FTCS

using Printf
using Plots; gr()

const alpha = 1.0
const N     = 50
const T     = 0.1
x  = LinRange(0, 1, N+2)
dx = x[2] - x[1]
r  = 0.2
dt = r * dx^2 / alpha
steps = Int(floor(T / dt)) + 1

X = [xi for xi in x, _ in x]
Y = [yi for _  in x, yi in x]

u     = sin.(π .* X) .* sin.(π .* Y)
u[1,:].=0; u[end,:].=0; u[:,1].=0; u[:,end].=0
u_new = copy(u)

exact2d(t) = exp(-2*alpha*π^2*t) .* sin.(π.*X) .* sin.(π.*Y)

snap_steps = Set([1, div(steps,2), steps])
snaps = Dict{Float64, Matrix{Float64}}()

for n in 1:steps
    n in snap_steps && (snaps[n*dt] = copy(u))
    @inbounds for j in 2:N+1, i in 2:N+1
        u_new[i,j] = u[i,j] + r*(u[i+1,j]+u[i-1,j]+u[i,j+1]+u[i,j-1]-4u[i,j])
    end
    u_new[1,:].=0; u_new[end,:].=0; u_new[:,1].=0; u_new[:,end].=0
    u, u_new = u_new, u
end

t_final = (steps-1)*dt
l2_err  = sqrt(dx^2 * sum((u .- exact2d(t_final)).^2))
@printf("2D Heat  N=%d  T=%.2f  L2 error = %.3e\n", N, T, l2_err)

mkpath("../results")
t_keys = sort(collect(keys(snaps)))
plots_list = []
for t_val in t_keys
    u_s = snaps[t_val]
    p = contourf(x, x, u_s', levels=20, cmap=:hot, colorbar=true,
                 title="t=$(round(t_val,digits=3))", xlabel="x", ylabel="y")
    push!(plots_list, p)
end
for t_val in t_keys
    err = snaps[t_val] .- exact2d(t_val)
    p = contourf(x, x, err', levels=20, cmap=:RdBu, colorbar=true,
                 title="Error t=$(round(t_val,digits=3))", xlabel="x", ylabel="y")
    push!(plots_list, p)
end

plt = plot(plots_list..., layout=(2,3), size=(1400,800),
           plot_title="Julia — 2D Heat (FTCS)  N=$N  L2=$(round(l2_err,sigdigits=2))")
savefig(plt, "../results/julia_heat2d.png")
println("Saved ../results/julia_heat2d.png")
