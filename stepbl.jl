function stepBL_exp!(u, v, u_new, v_new, ue1, ue2, dx, dy, nu)
    for j=2:ny-1
        Q = nu*dx/(u[j]*dy*dy)
        u_new[j] = Q*(u[j+1] + u[j-1]) - 
        (2*Q - 1)*u[j] + 1/u[j]*(ue2^2-ue1^2)/2 - 
        v[j]/u[j]*dx./dy*0.5*(u[j+1] - u[j-1])
    end
    for j =2:ny
        v_new[j] = v_new[j-1] - dy/(2.0*dx)*(u_new[j] - u[j] + u_new[j-1] - u[j-1])
    end
end

function stepBL_imp!(u,v,u_new,v_new,ue1,ue2,nu,dx,dy,ny, A, B)
    for j = 2:ny-1
        alf = nu*dx/(u[j]*dy^2)
        bet = v[j]*dx/(2*u[j]*dy)
        A[j-1,j-1] = 1+2*alf
        if j > 2
            A[j-1,j-2] = -alf
        end
        if j < ny-1
            A[j-1,j] = -alf
        end
        B[j-1] = u[j] - bet*(u[j+1]-u[j-1])+ (ue2^2-ue1^2)/(2*u[j])
        if j == ny-1
            B[j-1] = B[j-1] + alf*ue2
        end
    end

    sol = A\B
    u_new[2:end-1] .= sol
    u_new[end] = ue2
    for j = 2:ny
        v_new[j] = v_new[j-1] - dy/(2*dx)*(u_new[j]-u[j]+u_new[j-1]-u[j-1])
    end
end

function getBlThick(ys,us, ue)
    found = false
    idx = 1
    while !found
        if us[idx]/ue >= 0.99
            found = true
            # linearly interpolate
            x1 = us[idx-1]/ue
            x2 = us[idx]/ue
            return (ys[idx]-ys[idx-1])/(x2-x1)*(0.99-x1) + ys[idx-1]
        elseif idx == length(ys)
            return NaN
        else
            idx += 1
        end
    end
end

function bl_exp(xs, ys, ue, nu, mu, u_init, v_init)
    nx = length(xs)
    ny = length(ys)

    dy = diff(ys)[1]
    dx = diff(xs)
 
    u = zeros(nx,ny) # horizontal velocity
    v = zeros(nx,ny) # vertical velocity
    tws = zeros(nx) # wall 2D shear stress
    
    # boundary condition
    u[1,:] = u_init
    v[1,:] = v_init
    u[1,2:end] .= ue[1]
    u[:,end] .= ue
    
    sep = false
    i = 1
    while (i <= nx-1) & !sep
        @views stepBL_exp!(u[i,:], v[i,:], u[i+1,:], v[i+1,:],ue[i], ue[i+1], dx[i], dy, nu)
        tw = mu/(2*dy)*(4*u[i+1,2] - u[i+1,3])
        if tw >0
            i = i + 1
            tws[i] = tw
        else
            sep = true
            println("separated at i=$i")
        end
    end

    return u, v, tws
end

function bl_imp(xs, ys, ue, nu, mu, u_init, v_init)
    nx = length(xs)
    ny = length(ys)    

    dx = diff(xs)
    dy = diff(ys)[1]

    u = zeros(nx,ny) # horizontal velocity
    v = zeros(nx,ny) # vertical velocity
    tws = zeros(nx) # wall 2D shear stress

    # boundary condition
    u[1,:] = u_init
    v[1,:] = v_init
    u[1,2:end] .= ue[1]
    u[:,end] .= ue

    A = zeros(ny-2, ny-2)
    B = zeros(ny-2)
   
    sep = false
    i = 1
    while (i <= nx-1) & !sep
        @views stepBL_imp!(u[i,:], v[i,:], u[i+1,:], v[i+1,:], ue[i], ue[i+1], nu, dx[i], dy, ny, A, B)
        tw = mu./(2*dy)*(4 *u[i+1,2] - u[i+1,3])
        if tw > 0
            i = i + 1
            tws[i] = tw
        else
            sep = true
            println("separated at i=$i")
        end
    end

    return u, v, tws
end