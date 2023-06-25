export State, p, q

struct State
    q::Vector{Float64}
    p::Vector{Float64}
end

Base.copy(s::State) = State(q(s) |> copy, p(s) |> copy)

q(s::State) = s.q

p(s::State) = s.p