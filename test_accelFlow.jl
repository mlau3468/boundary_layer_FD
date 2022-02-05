using Plots
using DelimitedFiles
include("stepbl.jl")

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
dx = L/nx

ny = 100
ys = collect(range(0,stop=H,length=ny))
dy = H/ny

u = zeros(nx,ny)
v = zeros(nx,ny)

ue = 9.5 .+ xs./2
ue[1:Int(nx/2)] .= uinf

# boundary condition
u[1,2:end] .= uinf
u[:,end] .= ue
u[:,1] .= 0
v[1,:] .= 0
v[:,1] .= 0


Cf = 0
sep = false
i = 1
while (i <= nx-1) & !sep
    @views stepBL_exp!(u[i,:], v[i,:], u[i+1,:], v[i+1,:],ue[i], ue[i+1], dx, dy, nu)
    tw = mu./(2*dy).*(4 .*u[i+1,2] .- u[i+1,3])
    cf = tw./(0.5*rho*uinf^2)
    if cf >0
        global i = i + 1
        global Cf = Cf + cf * dx
    else
        global sep = true
        println("separated at i=$i")
    end
end
Cf = Cf/2
println("Cf=$Cf")

p = plot(u[end,:]./uinf, ys,label="explicit")

bl_thicks = zeros(nx)
for i = 1:nx
    bl_thicks[i] = getBlThick(ys,u[i,:],ue[i])
end
p2 = plot(xs,bl_thicks, label="explicit")

# --------------------------------------------------------------------------

# Implicit

nx = 1000
xs = collect(range(0,stop=L,length=nx))
dx = L/nx

ny = 100
ys = collect(range(0,stop=H,length=ny))
dy = H/ny

u = zeros(nx,ny)
v = zeros(nx,ny)

ue = 9.5 .+ xs./2
ue[1:Int(nx/2)] .= uinf


# boundary condition
u[1,2:end] .= uinf
u[:,end] .= ue
u[:,1] .= 0
v[1,:] .= 0
v[:,1] .= 0

A = zeros(ny-2, ny-2)
B = zeros(ny-2)

Cf = 0
sep = false
i = 1
while (i <= nx-1) & !sep
    @views stepBL_imp!(u[i,:], v[i,:], u[i+1,:], v[i+1,:], ue[i], ue[i+1], nu, dx, dy, ny, A, B)
    tw = mu./(2*dy).*(4 .*u[i+1,2] .- u[i+1,3])
    cf = tw./(0.5*rho*uinf^2)
    if cf >0
        global i = i + 1
        global Cf = Cf + cf * dx
    else
        global sep = true
        println("separated at i=$i")
    end
end
Cf = Cf/2
println("Cf=$Cf")

p = plot!(p, u[end,:]./uinf, ys,label="implicit")

bl_thicks = zeros(nx)
for i = 1:nx
    bl_thicks[i] = getBlThick(ys,u[i,:],ue[i])
end
p2 = plot!(p2, xs,bl_thicks, label="implicit")


# --------------------------------------------
# Plotting
p = plot!(p, xlims=(0,1.2), ylims=(0,0.06))
p = plot!(p, xlabel="U/Uinf", ylabel="y")
p2 = plot!(p2, xlabel="x", ylabel="δ")
display(p)
display(p2)