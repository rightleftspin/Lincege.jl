
struct BondCluster{V<:AbstractSet{<:Unsigned},H<:Unsigned} <: AbstractCluster
        expansion_vertices::V
        ghash::H
end
