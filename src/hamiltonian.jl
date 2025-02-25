struct Hamiltonian{P, M <: AbstractMass}
    density::P # Unnormalized log likelihood, gradient
    mass::M
end

(H::Hamiltonian)(q, a = similar(q)) = H.density(q, a)

v(H::Hamiltonian, p) = H.mass\p

kinetic(H::Hamiltonian, p) = 0.5 * dot(p, v(H, p)) + logabsdet(H.mass)

refresh(H::Hamiltonian) = rand(H.mass)