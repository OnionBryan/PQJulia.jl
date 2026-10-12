"""Zeroing of secret buffers the package owns; see SECURITY.md, "Secrets in memory"."""
module Wipe

export wipe!

# Empty asm that reads the buffer: the zeroing cannot be removed as a dead store.
barrier(p::Ptr{Cvoid}) = Base.llvmcall("""call void asm sideeffect "", "r,~{memory}"(ptr %0)
ret void""", Cvoid, Tuple{Ptr{Cvoid}}, p)

function wipe!(x::Array{T}) where {T}
    isbitstype(T) || throw(ArgumentError("wipe! needs plain-data elements, got $T"))
    GC.@preserve x begin
        p = Ptr{Cvoid}(pointer(x))
        ccall(:memset, Ptr{Cvoid}, (Ptr{Cvoid}, Cint, Csize_t), p, 0, sizeof(x))
        barrier(p)
    end
    x
end
wipe!(x::Vector{<:Array}) = (foreach(wipe!, x); x)

# The mutable storage under a wrapper (view, reshape, reinterpret, transpose, PermutedDimsArray):
# follow `parent` to the end. Only Array, Memory and BitArray storage is written to; anything
# else (CodeUnits of an immutable String, ranges, unknown array types) is refused.
storage(x::AbstractArray) = (p = parent(x); p === x ? x : storage(p))
writable(::Union{Array, Memory, BitArray}) = true
writable(_) = false
base_pointer(s::Union{Array, Memory}) = Ptr{Cvoid}(pointer(s))
base_pointer(s::BitArray) = Ptr{Cvoid}(pointer(s.chunks))

# Contiguous memory is cleared with memset over its element stride (Base.elsize, which exceeds
# sizeof(T) for padded types); any other array of numbers is zeroed element by element. In both
# cases the barrier receives the storage pointer, so the stores cannot be dropped.
function wipe!(x::AbstractArray{T}) where {T}
    isbitstype(T) || throw(ArgumentError("wipe! needs plain-data elements, got $T"))
    s = storage(x)
    writable(s) || throw(ArgumentError("wipe! cannot write to the storage of a $(typeof(x))"))
    if x isa StridedArray && strides(x) == Base.size_to_strides(1, size(x)...)
        GC.@preserve x s begin
            ccall(:memset, Ptr{Cvoid}, (Ptr{Cvoid}, Cint, Csize_t), pointer(x), 0, length(x) * Base.elsize(x))
            barrier(base_pointer(s))
        end
    elseif T <: Number
        fill!(x, zero(T))
        GC.@preserve s barrier(base_pointer(s))
    else
        throw(ArgumentError("wipe! cannot zero a non-contiguous $(typeof(x))"))
    end
    x
end
wipe!(x) = throw(ArgumentError("wipe! zeroes arrays of plain data, got $(typeof(x))"))
wipe!() = nothing
# Vararg{Any, N} and a recursive helper: compiled for each argument list, so every buffer is wiped
# through a static call. A plain `xs...` method is compiled for a widened Vararg{Any} signature,
# and `foreach` over it dispatches at run time (which JET does not report).
wipe!(a, b, xs::Vararg{Any, N}) where {N} = (wipe!(a); wipe!(b); wipe_each(xs...); nothing)
@inline wipe_each() = nothing
@inline wipe_each(x, xs::Vararg{Any, N}) where {N} = (wipe!(x); wipe_each(xs...))

end # module Wipe
