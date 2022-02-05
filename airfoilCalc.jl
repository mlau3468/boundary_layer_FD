include("airfoilTools.jl")
include("stepbl.jl")
include("kellersBox.jl")
using LinearAlgebra
using DelimitedFiles


function changeframe(tan1, norm1, tan2, norm2, vec)
    vec = vec[1].*tan1 .+ vec[2].*norm1 # to global

    x = tan2[1]*vec[1] + tan2[2]*vec[2] # to local
    z = norm2[1]*vec[1] + norm2[2]*vec[2]

    return [x;z]
end

function linVortPan(p1, p2, p, normal, tangent, panLen)
    # project point onto panel coordinates
    # warning: not valid when p lies on p1-p2
    # unit strengths at both panel points assumed
    @views x = tangent[1]*(p[1]-p1[1]) + tangent[2]*(p[2]-p1[2])
    @views z = normal[1]*(p[1]-p1[1]) + normal[2]*(p[2]-p1[2])

    theta1 = atan(z,x)
    theta2 = atan(z, x-panLen)
    r1 = sqrt(x^2+z^2)
    r2 = sqrt((x-panLen)^2+z^2)

    # first panel point
    up_a = -z/2/pi/panLen*log(r2/r1) + (panLen-x)/2/pi/panLen * (theta2-theta1)
    wp_a = -(panLen-x) /2/pi/panLen * log(r1/r2) -z/2/pi/panLen * (panLen/z-(theta2-theta1))

    # second panel point
    up_b = z/2/pi/panLen*log(r2/r1) + x/2/pi/panLen * (theta2-theta1)
    wp_b = -x/2/pi/(panLen)*log(r1/r2) + z/2/pi/(panLen) * (panLen/z-(theta2-theta1))

    # global frame
    uw_a = up_a.*tangent .+ wp_a.*normal
    uw_b = up_b.*tangent .+ wp_b.*normal

    return [uw_a'; uw_b']

end

function main()
    pts = genNACA([0,0,1,2], 150, true)
    npan = size(pts,2)-1
    c_pts, thetas, norms, tangents, panLen = procPanels(pts)

    # conditions
    U = 1
    chord = 1
    rho = 1.225
    alf = 0
    mu = 1.81*10^-5
    nu = mu/rho
    re = rho*U*chord/mu

    # Initialize solver matrix
    A = zeros(npan+1, npan+1) # normal matrix
    B = zeros(npan, npan+1) # tangent matrix
    RHS = zeros(npan+1)
    u_vec = U .* [cosd(alf), sind(alf)]

    # build RHS matrix
    for i = 1:npan
        RHS[i] = dot(-u_vec,norms[:, i])
    end

    # build influence coefficient matrix
    for i = 1:npan
        for j = 1:npan
            if i == j
                A[i,j] += -0.5/pi #wp_a
                B[i,j] += 0.25 #up_a
                A[i,j+1] += 0.5/pi #wp_b
                B[i,j+1] += 0.25 #up_b
            else
                res = linVortPan(pts[:,j],pts[:,j+1], c_pts[:,i], norms[:,j],tangents[:,j], panLen[j])
                uw_a = res[1,:]
                uw_b = res[2,:]

                A[i,j] += uw_a'norms[:,i]
                B[i,j] += uw_a'tangents[:,i]
                A[i,j+1] += uw_b'norms[:,i]
                B[i,j+1] += uw_b'tangents[:,i]
            end
        end
    end

    # Kutta condition
    A[end, 1] = 1
    A[end, end] = 1
    RHS[end] = 0

    # solve for circulation strength
    gam = A\RHS

    # lift by kutta jokowski
    lift = 0.0
    for i = 1:npan
        lift += rho*U.*(gam[i] + gam[i+1])/2*panLen[i]
    end
    # lift coefficient
    cl = lift/0.5/rho/U^2/chord

    # get pressure coefficient and velocities
    vs = zeros(npan)
    cp = zeros(npan)
    for i = 1:npan
        vs[i] = B[i,:]'gam + dot(u_vec,tangents[:,i])
        cp[i] = 1-vs[i]^2/U^2
    end

    # find stagnation point
    idx = 0
    found = false
    while !found && idx < npan
        idx += 1
        if vs[idx+1] >0 
            found = true
        end
    end

    botsurf = collect(1:idx)
    upsurf = collect(idx+1:npan)
    upsurf = upsurf[4:end-5]

    ny = 100
    H = 0.8

    xs = cumsum(panLen[upsurf])
    xs = xs .- xs[1]
    ue = vs[upsurf]

    #ys = collect(range(0,stop=H,length=ny))

    # boundary conditions
    u_init = zeros(ny)
    u_init[2:end] .= ue[1]
    v_init = zeros(ny)

    #u, v, tws = bl_imp(xs, ys, ue, nu, mu, u_init, v_init)

    ηe = 10
    eta = collect(LinRange(0,ηe, ny))

    # untransformed y coordinates
    nx = length(xs)
    ys = zeros(nx, ny)
    for i = 1:nx
        ys[i,:] .= eta/sqrt(ue[i]/(nu*xs[i]))
    end

    u, tws, f, f1, f2 = bl_kbox(xs, eta, ue, nu, mu)

    # get BL thickness
    nx = length(upsurf)
    bl_thicks = zeros(nx)
    for i = 1:nx
        bl_thicks[i] = getBlThick(ys[i,:],u[i,:],ue[i])
        bl_thicks[i] = getBlThick(ys[i,:],u[i,:],ue[i])
    end

    # calculate viscous force vector for each panel
    cfs = zeros(npan)
    tws_pan = zeros(npan)
    tws_pan[upsurf] .= tws
    dfs_vec = zeros(2,npan)
    for i = 1:npan
        dfs_vec[:,i] .= tws_pan[i].*panLen[i].*tangents[:,i]
        cfs[i] = tws_pan[i]/(0.5*rho*U^2*chord)
    end

    # sum viscous forces
    f_visc = sum(dfs_vec, dims=2)

    # drag coefficient
    cd = f_visc[1]/0.5/rho/U^2/chord

    println("CL: $cl")
    println("CD: $cd")

    cpPlot = plot(c_pts[1,:], cp, yflip=true, label="calc")

    xf_cp_i = readdlm("xfoil_ivisc_cp.txt", skipstart=3)
    cpPlot = plot!(cpPlot, xf_cp_i[:,1], xf_cp_i[:,3], label="xfoil")
    display(cpPlot)

    xf_data_i = readdlm("xfoil_ivisc_res.txt", skipstart=1)
    xf_data_v = readdlm("xfoil_visc_res.txt", skipstart=1)

    p1 = plot(c_pts[1,upsurf], bl_thicks)
    p2 = plot(c_pts[1,upsurf], cfs[upsurf], label="calc")
    p2 = plot!(p2, xf_data_v[:,2], xf_data_v[:,7], label="xfoil")
    p2 = plot!(p2, xlabel="x", ylabel="cf")

    p3 = plot(c_pts[1,:], vs, label="calc")
    p3 = plot!(p3, xf_data_i[:,2], xf_data_i[:,4], label="xfoil")
    p3 = plot!(p3, xlabel="x", ylabel="ue")

    display(p1)
    display(p2)
    display(p3)
end

main()