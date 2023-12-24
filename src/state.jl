export State, p, q

struct State{T <: AbstractVecOrMat}
    q::Vector{T}
    p::Vector{T}
end

Base.copy(s::State) = State(q(s) |> copy, p(s) |> copy)

q(s::State) = s.q

p(s::State) = s.p