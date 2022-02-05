using Plots
using DelimitedFiles
include("stepbl.jl")
include("kellersBox.jl")

# Settings
L = 2
H = 0.1

uinf = 10
rho = 1.225
nu = 0.0002
mu = nu*rho
re = uinf*L/nu

#--------------------------------------------------------
# Explicit
nx = 2000
xs = collect(range(0,stop=L,length=nx))
dx = diff(xs)[1]

ny = 100
ys = collect(range(0,stop=H,length=ny))
dy = diff(ys)[1]

ue = 10.5 .- xs./2
ue[1:Int(nx/2)] .= uinf

# bondary and initial conditions
u_init = zeros(ny)
u_init[end] = ue[1]
v_init = zeros(ny)

u, v, tws = bl_exp(xs, ys, ue, nu, mu, u_init, v_init)

Cf = sum(tws.*dx)/(0.5*rho*uinf^2*L)
println("Cf=$Cf")

p = plot(u[end,:]./uinf, ys,label="explicit")

bl_thicks = zeros(nx)
for i = 1:nx
    bl_thicks[i] = getBlThick(ys,u[i,:],ue[i])
end
p2 = plot(xs,bl_thicks, label="explicit")

cfs = tws./(0.5*rho*uinf^2)
p3 = plot(xs, cfs, label="explicit")


# --------------------------------------------------------------------------
# Implicit
nx = 1000
xs = collect(range(0,stop=L,length=nx))
dx = diff(xs)[1]

ny = 100
ys = collect(range(0,stop=H,length=ny))
dy = diff(ys)[1]

u = zeros(nx,ny)
v = zeros(nx,ny)

ue = 10.5 .- xs./2
ue[1:Int(nx/2)] .= uinf

# bondary and initial conditions
u_init = zeros(ny)
u_init[end] = ue[1]
v_init = zeros(ny)

u, v, tws = bl_imp(xs, ys, ue, nu, mu, u_init, v_init)

Cf = sum(tws.*dx)/(0.5*rho*uinf^2*L)
println("Cf=$Cf")

p = plot!(p, u[end,:]./uinf, ys,label="implicit")

bl_thicks = zeros(nx)
for i = 1:nx
    bl_thicks[i] = getBlThick(ys,u[i,:],ue[i])
end
p2 = plot!(p2, xs,bl_thicks, label="implicit")

cfs = tws./(0.5*rho*uinf^2)
p3 = plot!(p3, xs, cfs, label="implicit")


# --------------------------------------------------------------------------
# Keller's Box Method

np = 100 # number of net points, including at wall
nx = 500

xs = collect(LinRange(0,L,nx))
dx = diff(xs)

ue = 10.5 .- xs./2
ue[1:Int(nx/2)] .= uinf

ηe = 8
eta = collect(LinRange(0,ηe, np))

# untransformed y coordinates
ys = zeros(nx, np)
for i = 1:nx
    ys[i,:] .= eta/sqrt(ue[i]/(nu*xs[i]))
end

uvel, tws, f, f1, f2 = bl_kbox(xs, eta, ue, nu)
Cf = sum(tws[2:end].*dx)/(0.5*rho*uinf^2*L)
println("Cf=$Cf")

p = plot!(p, uvel[end,:]./uinf, ys[end,:],label="Kbox")

bl_thicks = zeros(nx)
for i = 1:nx
    bl_thicks[i] = getBlThick(ys[i,:],uvel[i,:],ue[i])
end
p2 = plot!(p2, xs,bl_thicks, label="Kbox")

cfs = tws./(0.5*rho*uinf^2)
p3 = plot!(p3, xs, cfs, label="Kbox")

# --------------------------------------------
# Plotting

# Data from Schetz
howarthdata = readdlm("Howarth Data.csv",',')
p = plot!(p,howarthdata[:,1], howarthdata[:,2],label="Schetz Pg. 111")

p = plot!(p, xlims=(0,1), ylims=(0,0.06))
p = plot!(p, xlabel="U/Uinf", ylabel="y")
p2 = plot!(p2, xlabel="x", ylabel="δ")
p3 = plot!(p3, xlabel="x", ylabel="cf(x)")
p3 = plot!(p3, xlims=(xs[2],L))
display(p)
display(p2)
display(p3)