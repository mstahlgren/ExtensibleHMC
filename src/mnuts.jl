# Multinominal No U-Turn Sampler with robust generalized U-turn criteria.

import LinearAlgebra: ⋅
import LogExpFunctions: logaddexp

struct MNUTS <: Sampler
    ϵ::Float64
    max_depth::Int
    max_ΔE::Float64
end

MNUTS(ϵ, max_depth = 10) = MNUTS(ϵ, max_depth, 1000)

buffsize(ϕ::MNUTS) = 3 * (2^(ϕ.max_depth+1)-1) + 5 # 5 => 63 => 193, 6 => 127 => 385, 7 => 255 => 769

struct BinaryTree{T}
    left::State{T}
    right::State{T}
    prop::State{T}
    msum::T
    psum::Float64
    esum::Float64
    steps::Int
end

function BinaryTree(θ::Hamiltonian, q₀, buffer)
    s = State(θ, q₀, buffer)
    display(q(s))
    display(p(s))
    display(copy!(pop!(buffer), a(s)))
    ls = State(q(s), p(s), copy!(pop!(buffer), a(s)), ll(s), ke(s))
    BinaryTree(ls, s, s, p(s), 0.0, 0.0, 0)
end

function BinaryTree(prop::State, l::BinaryTree, r::BinaryTree, buffer, esum = logaddexp(l.esum, r.esum))
    BinaryTree(l.left, r.right, prop, pop!(buffer) .= l.msum .+ r.msum, l.psum + r.psum, esum, l.steps + r.steps)
end

function BinaryTree(l::BinaryTree, r::BinaryTree, buffer)
    esum = logaddexp(l.esum, r.esum)
    BinaryTree(mh(esum, l.esum) ? l.prop : r.prop, l, r, buffer, esum)
end

mh(a, p) = rand() < exp(p - a)

uturn(msum, vₗ, vᵣ) = vₗ ⋅ msum < 0 || vᵣ ⋅ msum < 0

function uturn(θ::Hamiltonian, t::T, l::T, r::T, buffer) where T <: BinaryTree
    outer = uturn(t.msum, v(θ, p(t.left)), v(θ, p(t.right)))
    left  = uturn(peek(buffer) .= l.msum .+ p(r.left), v(θ, p(l.left)), v(θ, p(r.left)))
    right = uturn(peek(buffer) .= r.msum .+ p(l.right), v(θ, p(l.right)), v(θ, p(r.right)))
    return outer || left || right
end

function sample(ϕ::MNUTS, θ::Hamiltonian, q₀, buffer::Buffer = Buffer(ϕ, length(q₀)))
    tree = BinaryTree(θ, q₀, buffer)
    E₀, turned, div = energy(tree.prop), false, false
    for j = 0:ϕ.max_depth
        v = rand((-1, 1))
        if isone(v)
            tree′, turned, div = buildtree(ϕ, θ, tree.right, 1, j, E₀, buffer)
            ltree, rtree = tree, tree′
        else
            tree′, turned, div = buildtree(ϕ, θ, tree.left, -1, j, E₀, buffer)
            ltree, rtree = tree′, tree
        end
        accepted = !div && !turned && mh(tree.esum, tree′.esum)
        tree = BinaryTree(accepted ? tree′.prop : tree.prop, ltree, rtree, buffer)
        turned = turned || uturn(θ, tree, ltree, rtree, buffer)
        if div || turned break end
    end
    return Sample(copy(q(tree.prop)), ll(tree.prop), tree.steps, tree.psum / tree.steps, div, !div && !turned)
end

function buildtree(ϕ::MNUTS, θ, s, v, j, E₀, buffer)
    if iszero(j) return buildleaf(ϕ, θ, s, v * ϕ.ϵ, E₀, buffer) end
    tree₁, uturn₁, div₁ = buildtree(ϕ, θ, s, v, j - 1, E₀, buffer)
    if uturn₁ || div₁ return tree₁, uturn₁, div₁ end
    if isone(v)
        tree₂, uturn₂, div₂ = buildtree(ϕ, θ, tree₁.right, 1, j - 1, E₀, buffer)
        ltree, rtree = tree₁, tree₂
    else 
        tree₂, uturn₂, div₂ = buildtree(ϕ, θ, tree₁.left, -1, j - 1, E₀, buffer)
        ltree, rtree = tree₂, tree₁
    end
    tr = BinaryTree(ltree, rtree, buffer)
    return tr, uturn₂ || uturn(θ, tr, ltree, rtree, buffer), div₂
end

function buildleaf(ϕ::MNUTS, θ, s, ϵ, E₀, buffer)
    s′ = leapfrog(θ, s, ϵ, buffer)
    E′ = energy(s′)
    ΔE = E′ - E₀
    P = min(1.0, exp(-ΔE))
    d = ΔE > ϕ.max_ΔE
    return BinaryTree(s′, s′, s′, p(s′), P, -ΔE, 1), false, d
end
