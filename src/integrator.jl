import Zygote: gradient

export Leapfrog

abstract type Integrator end

struct Leapfrog{T <: Terminator} <: Integrator
    ϵ::Float64
    tc::T
end

# CONSIDER: Use pullback and grab ll when generating path
function (x::Leapfrog)(s, θ)
    s = copy(s)
    path = [copy(s)]
    df = ∇(θ)
    a = df(q(s))
    while true
        p(s) .-= 0.5 * x.ϵ * a
        q(s) .+= x.ϵ * (θ.mass\p(s))
        a = df(q(s))
        p(s) .-= 0.5 * x.ϵ * a
        push!(path, copy(s))
        x.tc(path) && break
    end
    return path
end