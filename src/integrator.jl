import Zygote: gradient

export Leapfrog

abstract type Integrator end

struct Leapfrog{T <: Terminator} <: Integrator
    ϵ::Float64
    tc::T
end

# TODO: Incorporate mass into hamiltonian equations
# CONSIDER: Use pullback and grab ll when generating path
function (x::Leapfrog)(q, p, f)
    path = [(q, p)]
    p -= 0.5 * x.ϵ * f(q)
    while !x.tc(path)
        q += x.ϵ * p
        p -= x.ϵ * f(q)
        push!(path, (q, p))
    end
    q += x.ϵ * p
    p -= 0.5 * x.ϵ * f(q)
    push!(path, (q, p))
    return path
end

# function (x::Leapfrog)(q, p, f::Function)
#     p -= 0.5 * x.ϵ * f(q)
#     for i in range(length = x.L - 1)
#         q += x.ϵ * p
#         p -= x.ϵ * f(q)
#     end
#     q += x.ϵ * p
#     p -= 0.5 * x.ϵ * f(q)
#     return q, p
# end