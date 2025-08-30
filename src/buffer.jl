mutable struct Buffer
    const values  ::Matrix{Float64}
    current       ::Int
end

Buffer(m::Int, n) = Buffer(Matrix{Float64}(undef, m, n), 1)

reset!(b::Buffer) = b.current = 1

function Base.pop!(b::Buffer)
    if b.current > size(b.values, 2) BoundsError(b, b.current) |> throw end
    b.current += 1
    return @inbounds view(b.values, :, b.current - 1)
end

function peek(b::Buffer)
    if b.current > size(b.values, 2) BoundsError(b, b.current) |> throw end
    return @inbounds view(b.values, :, b.current)
end