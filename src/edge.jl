import Base: eltype, Tuple, ==
import Graphs: AbstractEdge, src, dst, reverse
import Polymer: block

struct PolymerBlockEdge{T<:Real} <: AbstractEdge{Int}
    src::Int
    dst::Int
    block::PolymerBlock{T}

    function PolymerBlockEdge(src::Integer, dst::Integer, block::PolymerBlock{T}) where T
        return new{T}(src, dst, block)
    end
end

eltype(e::PolymerBlockEdge) = eltype(src(e))

# Accessors
src(e::PolymerBlockEdge) = e.src
dst(e::PolymerBlockEdge) = e.dst
weight(e::PolymerBlockEdge) = e.block.f  # block_length(e.block)
block(e::PolymerBlockEdge) = e.block

# I/O
show(io::IO, e::PolymerBlockEdge) = print(io, "Edge $(e.src) => $(e.dst) mapping to $(e.block)")

# Conversions
Tuple(e::PolymerBlockEdge) = (src(e), dst(e), weight(e))

(::Type{PolymerBlockEdge{T}})(e::PolymerBlockEdge) where T = PolymerBlockEdge{T}(e.src, e.dst, e.block)

# Convenience functions - note that these do not use weight.
reverse(e::T) where {T <: PolymerBlockEdge} = T(dst(e), src(e), block(e))
==(e1::PolymerBlockEdge, e2::PolymerBlockEdge) = (src(e1) == src(e2) && dst(e1) == dst(e2))
==(e1::PolymerBlockEdge, e2::AbstractEdge) = (src(e1) == src(e2) && dst(e1) == dst(e2))
==(e1::AbstractEdge, e2::PolymerBlockEdge) = (src(e1) == src(e2) && dst(e1) == dst(e2))

# IO
Base.show(io::IO, e::PolymerBlockEdge) = print(io, "Edge $(e.src) => $(e.dst) mapping to $(e.block)")
