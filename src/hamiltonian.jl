import LinearAlgebra: ⋅
import Zygote: pullback

export Hamiltonian, ∇, mass

struct Hamiltonian{P, M <: Mass}
    density::P # Unnormalized log posterior
    mass::M # LLᵀ decomposition
end

function gradient(H::Hamiltonian, s::AbstractArray) 
    ll, pb = pullback(H.density, s)
    return ll, pb(1.0)[1]
end

energy(H::Hamiltonian, s) = kinetic(H, s) + potential(H, s)

kinetic(H::Hamiltonian, s) = 0.5 * p(s)⋅(H.mass\p(s))

potential(H::Hamiltonian, s) = -H.density(q(s))

mass(H::Hamiltonian) = H.mass