struct Bond
    site1::Int
    site2::Int
    direction::Vector{Int}
    bond_type::Int
end

# coordinates are written [x1 x2 x3 ...; y1 y2 y3 ...; z1 z2 z3 ...]
struct UnitCell
    basis::Matrix{Float64}
    primitive_vectors::Matrix{Float64}
    site_colors::Vector{Int}
    bonds::Vector{Bond}
end

function UnitCell(basis::Vector{<:Vector{Float64}}, primitive_vectors::Vector{<:Vector{Float64}}, site_colors::Vector{Int}, bonds::Vector{Bond})
    UnitCell(hcat(basis...), hcat(primitive_vectors...), site_colors, bonds)
end

shift_unit_cell(unit_cell::UnitCell, shift_vector::Vector{Int}) = unit_cell.primitive_vectors * shift_vector[1:end-1] + unit_cell.basis[:, x[end]]
shift_unit_cell(unit_cell::UnitCell, shift_vectors::Matrix{Int}) = stack(v -> shift_unit_cell(unit_cell, v), eachcol(shift_vectors))
image_unit_cell(unit_cell::UnitCell) = _NI("image_unit_cell")
