abstract type AbstractUnitCell end

basis_size(unit_cell::AbstractUnitCell) = _NI("basis_size")
dimension(unit_cell::AbstractUnitCell) = _NI("dimension")
shift_unit_cell(unit_cell::AbstractUnitCell, shift_vector::AbstractVector{<:Int}) = _NI("shift_unit_cell")
# matrix input is of the form [shift_x shift_y shift_z; ...] where each row is a shift vector for the unit cell
shift_unit_cell(unit_cell::AbstractUnitCell, shift_vectors::AbstractMatrix{<:Int}) = _NI("shift_unit_cell")
"""
    image_unit_cell(unit_cell)
    image_unit_cell(expansion_unit_cell; show_expansion_bonds=false)

Return a Plots object showing the sites, physical bonds, and unlabeled primitive
cell boundary in two or three dimensions. Load `Plots` and `ColorSchemes` first.

For an `ExpansionUnitCell`, enable `show_expansion_bonds` to draw dashed
connections between expansion-element centers. Coincident sites with the same
color and duplicate physical bonds of the same type are drawn once, matching
coordinates rounded to ten decimal places. Neighboring physical bond endpoints
are also marked as sites in expansion-cell plots.
"""
function image_unit_cell end

"""
    image_lattice(unit_cell; cells=nothing)
    image_lattice(expansion_unit_cell; cells=-2:2, show_expansion_bonds=false)

Return a Plots object showing a repeated lattice in two or three dimensions.
Load `Plots` and `ColorSchemes` first. `cells` supplies the integer translations
along each primitive vector; `cells=0:1` repeats two cells per direction.
Ordinary `UnitCell` plots default to `-6:6` in 2D and `-2:2` in 3D.
`ExpansionUnitCell` plots default to `-2:2` in both dimensions.

For expansion cells, `show_expansion_bonds=true` adds dashed connections between
element centers. Physical bonds retain their types, shared sites retain their
colors, and coincident sites and duplicate physical bonds are drawn once as in
[`image_unit_cell`](@ref). Outgoing physical bonds include neighboring endpoints.
The returned object can be saved with `Plots.savefig`.
"""
function image_lattice end

include("UnitCell.jl")
include("ExpansionUnitCell.jl")
include("util.jl")
