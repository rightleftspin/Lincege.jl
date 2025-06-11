struct Tiling{B<:AbstractVector{<:AbstractVector{<:AbstractVector{<:Float}}},T<:AbstractVector{AbstractVector{<:Float}},L<:AbstractVector{AbstractVector{<:Integer}},N<:AbstractVector{<:Float},F}
        tiling_unit::B
        translation_vectors::T
        translation_labels::L
        expansion_neighbors::N
        real_space_neighbors::N

        # Optional
        labels::L
        neighbor_fn::F
end

function Tiling(
        basis::AbstractVector{<:AbstractVector{<:Float}},
        primitive_vectors::AbstractVector{AbstractVector{<:Float}},
        neighbors::AbstractVector{<:Float};
        labels::AbstractVector{<:Integer}=repeat([1], length(basis)),
        neighbor_fn::Function=((c1, c2, d, dir) -> nothing),
        expand_by_basis::Bool=false,
)

        tiling_unit = Vector{Vector{Vector{Float64}}}[]
        if expand_by_basis
                tiling_unit = [basis]
        else
                for basis_elem in basis
                        push!(tiling_unit, [basis_elem])
                end
        end

        translation_vectors = primitive_vectors
        translation_labels = [collect(1:length(basis))]
        opt_labels = [labels]

        Tiling(
                tiling_unit,
                translation_vectors,
                translation_labels,
                neighbors,
                neighbors,
                opt_labels,
                neighbor_fn
        )
end



real_space_neighbors(tiling::Tiling) = tiling.real_space_neighbors







