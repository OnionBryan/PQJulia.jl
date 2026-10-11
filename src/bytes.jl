"""
Byte-string concatenation for hash inputs and encodings. Base.vcat on byte vectors carries an
assertion whose message path is a runtime dispatch on Julia 1.12, which the dispatch gate
(test/timing/dispatch.jl) reports on every secret path; catbytes preallocates and copies.
"""
module Bytes

export catbytes

@inline nbytes(::UInt8) = 1
@inline nbytes(x::AbstractVector{UInt8}) = length(x)
@inline totalbytes() = 0
@inline totalbytes(x, xs...) = nbytes(x) + totalbytes(xs...)

@inline function put(r::Vector{UInt8}, o::Int, x::UInt8)
    @inbounds r[o+1] = x
    o + 1
end
@inline function put(r::Vector{UInt8}, o::Int, x::AbstractVector{UInt8})
    n = length(x)
    n == 0 || copyto!(r, o + 1, x, firstindex(x), n)
    o + n
end
@inline putall(r::Vector{UInt8}, o::Int) = r
@inline putall(r::Vector{UInt8}, o::Int, x, xs...) = putall(r, put(r, o, x), xs...)

"""
    catbytes(xs...) -> Vector{UInt8}

The concatenation of `xs`, each a `UInt8` or an `AbstractVector{UInt8}`, as a new `Vector{UInt8}`.
"""
catbytes(xs::Union{UInt8, AbstractVector{UInt8}}...) =
    putall(Vector{UInt8}(undef, totalbytes(xs...)), 0, xs...)

end # module Bytes
