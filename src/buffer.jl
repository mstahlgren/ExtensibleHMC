buffsize(sampler) = 3 * (2^(sampler.max_depth+1)-1) + 5 # 5 => 63 => 193, 6 => 127 => 385, 7 => 255 => 769

mutable struct Buffer
    values      ::Matrix{Float64}
    step        ::Int
    current     ::Int
end

Buffer(sampler, ex) = Buffer(Matrix{Float64}(undef, size(ex, 1), size(ex, 2) * buffsize(sampler)), size(ex, 2), 1)

reset!(b::Buffer) = b.current = 1

function Base.pop!(b::Buffer)
    if b.current > size(b.values, 2) throw(BoundsError(b, b.current)) end
    retrange = b.current:b.current+b.step-1
    b.current += b.step
    #println("$((b.current-1) / b.step), out of $(size(b.values, 2)/b.step)")
    return @inbounds view(b.values, :, retrange)
end

function peek(b::Buffer)
    if b.current > size(b.values, 2) throw(BoundsError(b, b.current)) end
    return @inbounds view(b.values, :, b.current:b.current+b.step-1)
end