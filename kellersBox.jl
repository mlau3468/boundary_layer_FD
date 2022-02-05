# Source: An Engineering Approach to the Calculation of Aeroynamic Flows, Cebeci, 1999

function lagrangeInterp(x1, x2, x3, y1, y2, y3, x0)
    return (x0-x2)*(x0-x3)/(x1-x2)/(x1-x3)*y1 + 
    (x0-x1)*(x0-x3)/(x2-x1)/(x2-x3)*y2+(x0-x1)*(x0-x2)/
    (x3-x1)/(x3-x2)*y3
end

function ivpl(η)
    # initial guess for falkner skan
    N = length(η)
    f = zeros(N) # f
    u = zeros(N) # f'
    v = zeros(N) # f''
    b = zeros(N) # (1+νt/ν) from eddy viscosity
    for j = 1:N
        etab = η[j]/η[N]
        etab2 = etab^2
        f[j] = 0.25*η[N]*etab*(3-0.5*etab2)
        u[j] = 0.5*etab*(3-etab2)
        v[j] = 1.5*(1-etab2)/η[N]
        b[j] = 1
    end
    return f, u, v, b
end

function fskvel(f, u ,v, ue)
    # f = f
    # u = f'
    # v = f''
    # computes velocity profile from falkner skan solution
    np = length(f)
    uvel = zeros(np)
    for n = 1:np
        uvel[n] = ue*u[n]
    end
    return uvel
end

function coef!(x, s1, s2, s3, s4, s5, s6, s7, s8, f, u, v, b, iter, nx, itv, np, deta, r1, r2, r3, p1, p2)

    cel = 0.0
    if nx > 1
        cel = 0.5*(x[nx]+x[nx-1])/(x[nx]-x[nx-1])
    end 
    λ1 = itv 
    λ2 = 1-itv
    p1p = λ1*p1[nx] + cel + 0.5*λ2

    for j = 2:np
    # present station
        usb = 0.5*(u[2,j]^2 + u[2,j-1]^2)
        ub = 0.5*(u[2,j]+u[2,j-1])
        vb = 0.5*(v[2,j]+v[2,j-1])
        fb = 0.5*(f[2,j]+f[2,j-1])
        fvb = 0.5*(f[2,j]*v[2,j]+f[2,j-1]*v[2,j-1])
        derbv = (b[2,j]*v[2,j]-b[2,j-1]*v[2,j-1])/deta[j]
        flare = 1
        if ub < 0
            flare = 0
        end
        p2p = λ1*p2[nx] + cel*flare
        if nx == 1
            crb = -p2[nx]
            cfb = 0
            cvb = 0
        else
            # previous station
            cfb = 0.5*(f[1,j]+f[1,j-1])
            cvb = 0.5*(v[1,j]+v[1,j-1])
            cusb = 0.5*(u[1,j]^2+u[1,j-1]^2)
            cfvb = 0.5*(f[1,j]*v[1,j]+f[1,j-1]*v[1,j-1])
            cderbv = (b[1,j]*v[1,j]-b[1,j-1]*v[1,j-1])/deta[j]
            clb = cderbv + λ1*(p1[nx]*cfvb+p2[nx]*(1-cusb)) + 0.5*λ2*cfvb
            crb = -clb - λ1*p2[nx] - cel*flare*cusb + cel*cfvb
        end

        # Sj
        s1[j] = b[2,j]/deta[j] + 0.5*p1p*f[2,j] - 0.5*cel*cfb
        s2[j] = -b[2,j-1]/deta[j] + 0.5*p1p*f[2,j-1] - 0.5*cel*cfb
        s3[j] = 0.5*(p1p*v[2,j] + cel*cvb)
        s4[j] = 0.5*(p1p*v[2,j-1] + cel*cvb)
        s5[j] = -p2p*u[2,j]
        s6[j] = -p2p*u[2,j-1]
        s7[j] = λ2*cel*u[2,np]
        s8[j] = λ2*cel*u[2,np]

        # Rj
        r2[j] = crb - (derbv + p1p*fvb - p2p*usb + λ2*cel*(u[2,np]^2-u[1,np]^2) + cel*(fb*cvb-vb*cfb))
        r1[j] = f[2,j-1] - f[2,j] + ub*deta[j]
        r3[j-1] = u[2,j-1] - u[2,j] + vb*deta[j]
    end

    # boundary conditions
    r1[1] = 0
    r2[np] = 0
    if itv == 1 #standard mode
        gamma1 = 0
        gamma2 = 1
        r3[np] = 0
    else # inverse mode
        gamma1 = cc[nx,nx] * sqrt(x[nx])
        gamma2 = 1 - eta[np] * gamma1
        r3[np] = gi - (gamma2*u[2,np] + gamma1*f[2,np])
    end
end

function bl_kbox(x, eta, ue, nu, mu)
    np = length(eta)
    Nx = length(x)

    dx = diff(x)
    deta = zeros(np)
    deta[2:end] = diff(eta) #deta[1] is dummy value so indexing is cosistent with Cebeci

    # array to save values
    twsave = zeros(Nx)
    usave = zeros(Nx, np)

    itermax = 50 # maximum number of newton iterations
    itv = 1 # standard mode

    # set up Cebeci
    s1 = zeros(np)
    s2 = zeros(np)
    s3 = zeros(np)
    s4 = zeros(np)
    s5 = zeros(np)
    s6 = zeros(np)
    s7 = zeros(np)
    s8 = zeros(np)
    r1 = zeros(np)
    r2 = zeros(np)
    r3 = zeros(np)
    f = zeros(2,np) #f [[prev], [current]]
    u = zeros(2, np) #f' [[prev], [current]]
    v = zeros(2, np) # f'' [[prev], [current]]
    b = zeros(2, np)
    p1 = zeros(Nx) # freestream velocity parameter m
    p2 = zeros(Nx) # freestream velocity parameter (m+1)/2

    J = np-1
    A = zeros(3*J+3, 3*J+3)
    R = zeros(3*J+3)

    # calculate freestream velocity parameters
    p2[1] = 0
    p1[1] = (p2[1]+1)/2
    for i = 2:Nx
       #duedx = lagrangeInterp(x[i-1],x[i],x[i+1],ue[i-1],ue[i],ue[i+1],x[i])
       duedx = (ue[i]-ue[i-1])/(x[i]-x[i-1])
       p2[i] = x[i]/ue[i]*duedx
       p1[i] = (p2[i]+1)/2
    end

    # initialize solution: present f, u, v, b = initial values
    f_init,u_init,v_init,b_init = ivpl(eta)

    f[2,:] .= f_init[:]
    u[2,:] .= u_init[:]
    v[2,:] .= v_init[:]
    b[2,:] .= b_init[:]

    go = true
    nx = 1
    while go && nx <= Nx

        # prev station is the current station
        f[1,:] .= f[2,:]
        u[1,:] .= u[2,:]
        v[1,:] .= v[2,:]
        b[1,:] .= b[2,:]

        # Newton Iteration
        converged = false
        iter = 1 # interation counter
        while (~converged) && (iter<itermax)

            coef!(x, s1, s2, s3, s4, s5, s6, s7, s8, f, u, v, b, iter, nx, itv, np, deta, r1, r2, r3, p1, p2)

            # first block j=0
            A[1:3,1:3] = [1 0 0; 0 1 0; 0 -1 -deta[2]/2] #A0
            A[1:3, 4:6] = [0 0 0; 0 0 0; 0 1 -deta[2]/2] #C0
            R[1:3] = [0;0;r3[1]] #R0

            # 1<j<J
            for j = 2:np
                if j == np
                    aj = [1 -deta[j]/2 0; s3[j] s5[j] s1[j]; 0 1 0]
                    bj = [-1 -deta[j]/2 0; s4[j] s6[j] s2[j]; 0 0 0]
                    rj = [r1[j]; r2[j];0]

                    ida = 3*(j-1)+1
                    A[ida:ida+2, ida:ida+2] = aj
                    A[ida:ida+2, ida-3:ida-1] = bj
                    R[ida:ida+2] = rj
                else
                    aj = [1 -deta[j]/2 0; s3[j] s5[j] s1[j]; 0 -1 -deta[j+1]/2]
                    bj = [-1 -deta[j]/2 0; s4[j] s6[j] s2[j]; 0 0 0]
                    cj = [0 0 0; 0 0 0; 0 1 -deta[j+1]/2]
                    rj = [r1[j];r2[j];r3[j]]

                    ida = 3*(j-1)+1
                    A[ida:ida+2, ida:ida+2] = aj
                    A[ida:ida+2, ida-3:ida-1] = bj
                    A[ida:ida+2, ida+3:ida+5] = cj
                    R[ida:ida+2] = rj
                end

            end

            del = A\R

            # update iterates
            for j = 1:np
                f[2,j] = f[2,j] + del[(j-1)*3+1]
                u[2,j] = u[2,j] + del[(j-1)*3+2]
                v[2,j] = v[2,j] + del[(j-1)*3+3]
            end

            # test for convergence
            if maximum(abs.(del)) < 0.00001
                converged = true
            end
            iter += 1
        end

        if converged
            # save values
            for n = 1:np
                usave[nx,n] = ue[nx]*u[2,n]
            end

            if nx > 1
                tw =  mu*sqrt(ue[nx]^3/(nu*x[nx]))*v[2,1]
                if tw < 0
                    failx = x[nx]
                    println("Negative shear stress at x=$failx.")
                    go = false
                else
                    twsave[nx] = tw
                end
            end
        else
            failx = x[nx]
            println("Not converged after $iter iterations at x=$failx.")
            go = false
        end

        nx += 1
    end
    return usave, twsave, f, u, v
end