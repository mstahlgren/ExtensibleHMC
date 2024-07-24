import LinearAlgebra: ⋅
import Zygote: withgradient

export Hamiltonian, ∇, mass

struct Hamiltonian{P, M <: Mass}
    density::P # Unnormalized log posterior
    mass::M # LLᵀ decomposition
end

(H::Hamiltonian)(s) = withgradient(H.density, s)

energy(H::Hamiltonian, s) = kinetic(H, s) + potential(H, s)

kinetic(H::Hamiltonian, s) = 0.5 * p(s)⋅(H.mass\p(s))

potential(::Hamiltonian, s) = -ll(s)

mass(H::Hamiltonian) = H.mass