
function averaged_reduced_state_dynamics_state_partial(
    p, f, g, h, k, m, t, δ, α1, α2, α3, μ, tMax, Isp, g0, τs, ws)

    # Construct thrust direction vector
    α       = SVector(α1, α2, α3)

    # Allocate storage for integrand partials
    dIdx    = zeros(5,6)

    # Approximate integral with quadrature
    L0      = -π + 1e-6
    Lf      =  π - 1e-6
    ddxdx   = zeros(6,6)
    @inbounds for i in eachindex(τs)
        # Compute true longitude for current point in quadrature
        Li      = 0.5*(Lf - L0)*τs[i] + 0.5*(Lf + L0)

        # Compute averaged reduced dynamics integrand state partials
        xi      = SVector(p,f,g,h,k,Li,m)
        averaged_reduced_state_dynamics_integrand_state_partials!(
            dIdx, xi, t, δ, α, μ, tMax, Isp, g0)

        # Compute integrand
        dIdx  .*= 0.5*(Lf - L0)

        # Update Jacobian approximation
        for col in axes(ddxdx,2)
            for row in 1:5
                ddxdx[row,col] += ws[i]*dIdx[row,col]
            end
        end
    end

    # Place in single tuple and return
    return Tuple(ddxdx)
end

function averaged_reduced_state_dynamics_control_partial(
    p, f, g, h, k, m, t, δ, α1, α2, α3, μ, tMax, Isp, g0, τs, ws)

    # Construct thrust direction vector
    α       = SVector(α1, α2, α3)

    # Allocate storage for integrand partials
    dIdu    = zeros(5,4)

    # Approximate integral with quadrature
    L0      = -π
    Lf      =  π
    ddxdu   = zeros(6,4)
    @inbounds for i in eachindex(τs)
        # Compute true longitude for current point in quadrature
        Li      = 0.5*(Lf - L0)*τs[i] + 0.5*(Lf + L0)

        # Compute averaged reduced dynamics integrand state partials
        xi      = SVector(p,f,g,h,k,Li,m)
        averaged_reduced_state_dynamics_integrand_control_partials!(
            dIdu, xi, t, δ, α, μ, tMax, Isp, g0)

        # Compute integrand
        dIdu  .*= 0.5*(Lf - L0)

        # Update Jacobian approximation
        for col in axes(ddxdu,2)
            for row in 1:5
                ddxdu[row,col] += ws[i]*dIdu[row,col]
            end
        end
    end

    # Compute partials for mass dynamics
    ddxdu[6,1]  = -tMax / (Isp*g0)

    # Place in single tuple and return
    return Tuple(ddxdu)
end

function averaged_reduced_state_dynamics_integrand_state_partials!(dIdx, xi, t, δ, α, μ, tMax, Isp, g0)
    # Grab states
    p,f,g,h,k,L,m = xi

    # Compute orbital period for current osculating orbit
    a       = p / (1.0 - f*f - g*g)
    T0      = 2.0*π*sqrt(a*a*a / μ)

    # Compute MEE dynamics
    dmee    = mee_dynamics(xi, t, δ, α, μ, tMax)
    dx      = view(dmee, 1:5)

    # Compute MEE dynamics partial
    ddmeedx = mee_dynamics_partials(xi, t, δ, α, μ, tMax)
    ddoedx  = view(ddmeedx, 1:5, 1:7 .!= 6)

    # Compute s
    w       = 1.0 + f*cos(L) + g*sin(L)
    dLdtp   = sqrt(μ*p)*(w / p)^2
    s       = 1.0 / (T0 * dLdtp)

    # Compute s partials
    dsdx    = view(s_partials(p, f, g, h, k, L, m, μ), 1:7 .!= 6)

    # Compute integrand partial
    #dIdx    = s*ddoedx + dmee*transpose(dsdx)
    dIdx   .= s*ddoedx
    mul!(dIdx, dx, transpose(dsdx), 1.0, 1.0)
    return nothing
end

function averaged_reduced_state_dynamics_integrand_control_partials!(dIdu, xi, t, δ, α, μ, tMax, Isp, g0)
    # Grab states
    p,f,g,h,k,L,m = xi

    # Compute orbital period for current osculating orbit
    a       = p / (1.0 - f*f - g*g)
    T0      = 2.0*π*sqrt(a*a*a / μ)

    # Compute MEE dynamics
    dmee    = mee_dynamics(xi, t, δ, α, μ, tMax)
    dx      = view(dmee, 1:5)

    # Compute s
    w       = 1.0 + f*cos(L) + g*sin(L)
    dLdtp   = sqrt(μ*p)*(w / p)^2
    s       = 1.0 / (T0 * dLdtp)

    # Compute partial wrt scale factor
    dIdu[:,1]  .= SVector(s*dx[1] / δ,
                          s*dx[2] / δ,
                          s*dx[3] / δ,
                          s*dx[4] / δ,
                          s*dx[5] / δ) 

    # Compute partial wrt thrust direction
    Boe         = compute_Boe(xi, μ)
    dIdu[:,2:4].= (δ*tMax*s / m) * Boe

    # Return nothing
    return nothing
end

function compute_Boe(x,μ)
    # Grab states
    p,f,g,h,k,L,m = x

    # Construct Boe matrix
    sqrtpmu = sqrt(p / μ)
    w       = 1.0 + f*cos(L) + g*sin(L)
    κ       = h*sin(L) - k*cos(L)
    ss      = 1.0 + h^2 + k^2
    wInv    = 1.0 / w
    SMatrix{5,3}(0.0, sqrtpmu*sin(L), -sqrtpmu*cos(L), 0.0, 0.0,
                 2.0*p*wInv*sqrtpmu, sqrtpmu*wInv*((w + 1.0)*cos(L) + f),
                 sqrtpmu*wInv*((w + 1.0)*sin(L) + g), 0.0, 0.0,
                 0.0, -sqrtpmu*g*wInv*κ, sqrtpmu*f*wInv*κ, 
                 0.5*sqrtpmu*ss*cos(L)*wInv, 0.5*sqrtpmu*ss*sin(L)*wInv)
end