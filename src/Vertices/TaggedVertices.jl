struct TaggedVertices{Tag,V<:Integer} <: AbstractVertices{V}
        bitset::BitSet
end

vertices(vertex_set::TaggedVertices) = vertex_set.bitset

Base.sort(vertex_set::TaggedVertices) = vertex_set
Base.collect(vertex_set::TaggedVertices) = collect(vertex_set.bitset)
Base.hash(vertex_set::TaggedVertices, h::UInt) = hash(vertex_set.bitset, h)

Base.intersect(vertex_set1::T, vertex_set2::T) where {T<:TaggedVertices} = T(intersect(vertex_set1.bitset, vertex_set2.bitset))
Base.setdiff(vertex_set1::T, vertex_set2::T) where {T<:TaggedVertices} = T(setdiff(vertex_set1.bitset, vertex_set2.bitset))
Base.union(vertex_set1::T, vertex_set2::T) where {T<:TaggedVertices} = T(union(vertex_set1.bitset, vertex_set2.bitset))
function Base.union(vertex_set::T, itr) where {T<:TaggedVertices}
        bs = copy(vertex_set.bitset)
        for x in itr
                union!(bs, x.bitset)
        end
        T(bs)
end

Base.in(v::Integer, vertex_set::TaggedVertices) = v in vertex_set.bitset
Base.eltype(::TaggedVertices{Tag,V}) where {Tag,V} = V

TaggedVertices{Tag,V}() where {Tag,V<:Integer} = TaggedVertices{Tag,V}(BitSet())
TaggedVertices{Tag,V}(i::Int) where {Tag,V<:Integer} = TaggedVertices{Tag,V}(BitSet(i))
TaggedVertices{Tag,V}(itr::AbstractVector{<:Integer}) where {Tag,V<:Integer} = TaggedVertices{Tag,V}(BitSet(itr))

"""
    ExpansionVertices{V}

A set of expansion-lattice vertices (`TaggedVertices` tagged `:expansion`), indexing
sites of a cluster-expansion lattice's expansion unit cell.
"""
const ExpansionVertices{V} = TaggedVertices{:expansion,V}

ExpansionVertices() = ExpansionVertices{Int}(BitSet())
ExpansionVertices(i::Int) = ExpansionVertices{Int}(BitSet(i))
ExpansionVertices(itr::AbstractVector{<:Integer}) = ExpansionVertices{Int}(BitSet(itr))
ExpansionVertices(bs::BitSet) = ExpansionVertices{Int}(bs)

"""
    LatticeVertices{V}

A set of lattice vertices (`TaggedVertices` tagged `:lattice`), indexing physical
sites of a lattice.
"""
const LatticeVertices{V} = TaggedVertices{:lattice,V}

LatticeVertices() = LatticeVertices{Int}(BitSet())
LatticeVertices(i::Int) = LatticeVertices{Int}(BitSet(i))
LatticeVertices(itr::AbstractVector{<:Integer}) = LatticeVertices{Int}(BitSet(itr))
LatticeVertices(bs::BitSet) = LatticeVertices{Int}(bs)
