# 3D Heat equation: u_t = alpha*(u_xx+u_yy+u_zz)
# Domain: [0,1]^3, Dirichlet BC, IC: sin(pi*x)*sin(pi*y)*sin(pi*z)
# Exact: exp(-3*alpha*pi^2*t)*sin(pi*x)*sin(pi*y)*sin(pi*z)
# Method: explicit FTCS

using Printf
using Plots; gr()

const alpha = 1.0
const N     = 20
const T     = 0.05
x  = LinRange(0, 1, N+2)
dx = x[2] - x[1]
r  = 0.13
dt = r * dx^2 / alpha
steps = Int(floor(T / dt)) + 1

X = [xi for xi in x, _ in x, __ in x]
Y = [yi for _  in x, yi in x, __ in x]
Z = [zi for _  in x, __ in x, zi in x]

u     = sin.(π.*X) .* sin.(π.*Y) .* sin.(π.*Z)
u[1,:,:].=0; u[end,:,:].=0; u[:,1,:].=0; u[:,end,:].=0
u[:,:,1].=0; u[:,:,end].=0
u_new = copy(u)

exact3d(t) = exp(-3*alpha*π^2*t) .* sin.(π.*X) .* sin.(π.*Y) .* sin.(π.*Z)

snap_steps = Set([1, div(steps,2), steps])
snaps = Dict{Float64, Array{Float64,3}}()

for n in 1:steps
    n in snap_steps && (snaps[n*dt] = copy(u))
    @inbounds for k in 2:N+1, j in 2:N+1, i in 2:N+1
        u_new[i,j,k] = u[i,j,k] + r*(
            u[i+1,j,k]+u[i-1,j,k]+u[i,j+1,k]+u[i,j-1,k]+
            u[i,j,k+1]+u[i,j,k-1]-6u[i,j,k])
    end
    u_new[1,:,:].=0; u_new[end,:,:].=0; u_new[:,1,:].=0; u_new[:,end,:].=0
    u_new[:,:,1].=0; u_new[:,:,end].=0
    u, u_new = u_new, u
end

t_final = (steps-1)*dt
l2_err  = sqrt(dx^3 * sum((u .- exact3d(t_final)).^2))
@printf("3D Heat  N=%d  T=%.3f  L2 error = %.3e\n", N, T, l2_err)

mkpath("../results")
mid = div(N,2) + 1
t_keys = sort(collect(keys(snaps)))
plots_list = []
for t_val in t_keys
    u_s = snaps[t_val][:,:,mid]
    p = contourf(x, x, u_s', levels=20, cmap=:hot, colorbar=true,
                 title="z=0.5  t=$(round(t_val,digits=3))", xlabel="x", ylabel="y")
    push!(plots_list, p)
end
for t_val in t_keys
    err = (snaps[t_val] .- exact3d(t_val))[:,:,mid]
    p = contourf(x, x, err', levels=20, cmap=:RdBu, colorbar=true,
                 title="Error z=0.5  t=$(round(t_val,digits=3))", xlabel="x", ylabel="y")
    push!(plots_list, p)
end
plt = plot(plots_list..., layout=(2,3), size=(1400,800),
           plot_title="Julia — 3D Heat (FTCS) z=0.5  N=$N  L2=$(round(l2_err,sigdigits=2))")
savefig(plt, "../results/julia_heat3d.png")
println("Saved ../results/julia_heat3d.png")
