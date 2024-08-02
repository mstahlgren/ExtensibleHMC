import LinearAlgebra: ⋅
import LogExpFunctions: logaddexp

export MNUTS

# Multinominal No U-Turn Sampler with robust generalized U-turn criteria.

struct MNUTS <: Sampler
    ϵ::Float64
    max_depth::Int
    max_ΔE::Float64
end

MNUTS(ϵ) = MNUTS(ϵ, 12, 1000)

struct BinaryTree{T}
    left::State{T}
    right::State{T}
    prop::State{T}
    msum::T
    esum::Float64
    steps::Int
end

function BinaryTree(θ, q)
    s = State(θ, q)
    BinaryTree(s, s, s, p(s), energy(s), 1)
end

function BinaryTree(prop::State, l::BinaryTree, r::BinaryTree, esum = logaddexp(l.esum, r.esum))
    BinaryTree(l.left, r.right, prop, l.msum + r.msum, esum, l.steps + r.steps)
end

function BinaryTree(l::BinaryTree, r::BinaryTree)
    esum = logaddexp(l.esum, r.esum)
    BinaryTree(mh(esum, l.esum) ? l.prop : r.prop, l, r, esum)
end

mh(a, p) = rand() < exp(p - a)

uturn(msum, vₗ, vᵣ) = vₗ ⋅ msum < 0 || vᵣ ⋅ msum < 0

function uturn(θ::Hamiltonian, t::T, l::T, r::T) where T <: BinaryTree
    outer = uturn(t.msum, v(θ, t.left), v(θ, t.right))
    left  = uturn(l.msum .+ p(r.left), v(θ, l.left), v(θ, r.left))
    right = uturn(r.msum .+ p(l.right), v(θ, l.right), v(θ, r.right))
    return outer || left || right
end

function StatsBase.sample(ϕ::MNUTS, θ::Hamiltonian, q₀)
    tree, j, turned, div = BinaryTree(θ, q₀), 0, false, false
    while j < ϕ.max_depth && !turned && !div
        v = rand((-1, 1))
        if isone(v)
            tree′, turned, div = buildtree(ϕ, θ, tree.right, 1, j)
            ltree, rtree = tree, tree′
        else
            tree′, turned, div = buildtree(ϕ, θ, tree.left, -1, j)
            ltree, rtree = tree′, tree
        end
        accepted = !div && !turned && mh(energy(tree.prop), energy(tree′.prop))
        tree = BinaryTree(accepted ? tree′.prop : tree.prop, ltree, rtree)
        turned = turned || uturn(θ, tree, ltree, rtree)
        if !div && !turned j = j + 1 end
    end
    return Sample(q(tree.prop), ll(tree.prop), tree.steps, q₀ != q(tree.prop), div, j == ϕ.max_depth)
end

function buildtree(ϕ::MNUTS, θ, s, v, j)
    if iszero(j) return buildleaf(ϕ, θ, s, v * ϕ.ϵ) end
    tree₁, uturn₁, div₁ = buildtree(ϕ, θ, s, v, j - 1)
    if uturn₁ || div₁ return tree₁, uturn₁, div₁ end
    if isone(v)
        tree₂, uturn₂, div₂ = buildtree(ϕ, θ, tree₁.right, 1, j - 1)
        ltree, rtree = tree₁, tree₂
    else 
        tree₂, uturn₂, div₂ = buildtree(ϕ, θ, tree₁.left, -1, j - 1)
        ltree, rtree = tree₂, tree₁
    end
    tr = BinaryTree(ltree, rtree)
    return tr, uturn₂ || uturn(θ, tr, ltree, rtree), div₂
end

function buildleaf(ϕ::MNUTS, θ, s, ϵ)
    s′ = leapfrog(θ, s, ϵ)
    E′ = energy(s′)
    d = E′ - energy(s) > ϕ.max_ΔE
    return BinaryTree(s′, s′, s′, p(s′), E′, 1), false, d
end