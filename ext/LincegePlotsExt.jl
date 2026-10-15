module LincegePlotsExt

using Plots
using ColorSchemes
using Lincege

function Lincege.image_unit_cell(unit_cell::UnitCell)
        dim = Lincege.dimension(unit_cell)
        if dim == 2
                return _image_unit_cell_2d(unit_cell)
        elseif dim == 3
                gr()
                return _image_unit_cell_3d(unit_cell)
        else
                error("Only 2D and 3D unit cells can be visualized. Got dimension: $dim")
        end
end

function _image_unit_cell_2d(unit_cell::UnitCell)
        sites_info = Vector{Float64}[]

        max_x = 0.1
        min_x = -0.1
        max_y = 0.1
        min_y = -0.1

        pv1 = unit_cell.primitive_vectors[:, 1]
        pv2 = unit_cell.primitive_vectors[:, 2]

        unit_cell_lines = [
                [[0, 0], pv1],
                [[0, 0], pv2],
                [pv2, pv1 .+ pv2],
                [pv1, pv1 .+ pv2]
        ]

        unit_cell_plot = plot(legend=:best)

        for i in 1:Lincege.basis_size(unit_cell)
                basis_pos = unit_cell.basis[:, i]
                push!(sites_info, [basis_pos[1], basis_pos[2], unit_cell.site_colors[i]])
        end

        for bond in unit_cell.bonds
                dir_x = bond.direction[1]
                dir_y = bond.direction[2]

                basis1 = unit_cell.basis[:, bond.site1]
                basis2 = unit_cell.basis[:, bond.site2]

                site1_x = basis1[1]
                site1_y = basis1[2]
                site2_x = basis2[1] + dir_x * pv1[1] + dir_y * pv2[1]
                site2_y = basis2[2] + dir_x * pv1[2] + dir_y * pv2[2]

                max_x = max(max_x, site1_x, site2_x)
                min_x = min(min_x, site1_x, site2_x)
                max_y = max(max_y, site1_y, site2_y)
                min_y = min(min_y, site1_y, site2_y)

                bond_palette = colorschemes[:seaborn_bright]
                plot!(unit_cell_plot, [site1_x, site2_x], [site1_y, site2_y],
                        label="", color=bond_palette[bond.bond_type], lw=2, dpi=1000)
        end

        min_x = min(min_x, pv1[1], pv2[1])
        max_x = max(max_x, pv1[1] + pv2[1], pv1[1], pv2[1])
        min_y = min(min_y, pv1[2], pv2[2])
        max_y = max(max_y, pv1[2] + pv2[2])

        x_positions = getindex.(sites_info, 1)
        y_positions = getindex.(sites_info, 2)
        colors = getindex.(sites_info, 3)
        atom_palette = colorschemes[:Accent_8]

        for c in unique(colors)
                mask = colors .== c
                color_index = Int(c)
                plot!(unit_cell_plot, x_positions[mask], y_positions[mask],
                        seriestype=:scatter,
                        aspect_ratio=1,
                        xlims=(min_x * 1.1, max_x * 1.1),
                        ylims=(min_y * 1.1, max_y * 1.1),
                        color=atom_palette[(color_index+1)],
                        label="Type $color_index",
                        markerstrokewidth=2,
                        markersize=6,
                        dpi=1000)
        end

        for line in unit_cell_lines
                x_points = [line[1][1], line[2][1]]
                y_points = [line[1][2], line[2][2]]
                plot!(unit_cell_plot, x_points, y_points, linestyle=:dash, color="black",
                        label="")
        end

        display(unit_cell_plot)
        return unit_cell_plot
end

function _image_unit_cell_3d(unit_cell::UnitCell)
        sites_info = Vector{Float64}[]
        max_x = 0.1
        min_x = -0.1
        max_y = 0.1
        min_y = -0.1
        max_z = 0.1
        min_z = -0.1

        pv1 = unit_cell.primitive_vectors[:, 1]
        pv2 = unit_cell.primitive_vectors[:, 2]
        pv3 = unit_cell.primitive_vectors[:, 3]

        unit_cell_lines = [
                [[0, 0, 0], pv1], [[0, 0, 0], pv2], [pv1, pv1 .+ pv2], [pv2, pv1 .+ pv2],
                [pv3, pv3 .+ pv1], [pv3, pv3 .+ pv2], [pv3 .+ pv1, pv3 .+ pv1 .+ pv2],
                [pv3 .+ pv2, pv3 .+ pv1 .+ pv2], [[0, 0, 0], pv3], [pv1, pv1 .+ pv3],
                [pv2, pv2 .+ pv3], [pv1 .+ pv2, pv1 .+ pv2 .+ pv3]
        ]

        unit_cell_plot = plot3d(legend=:outertopright)

        for i in 1:Lincege.basis_size(unit_cell)
                basis_pos = unit_cell.basis[:, i]
                push!(sites_info, [basis_pos[1], basis_pos[2], basis_pos[3], unit_cell.site_colors[i]])
        end

        for bond in unit_cell.bonds
                basis1 = unit_cell.basis[:, bond.site1]
                basis2 = unit_cell.basis[:, bond.site2]
                dir_x, dir_y, dir_z = bond.direction[1], bond.direction[2], bond.direction[3]

                site1_x, site1_y, site1_z = basis1[1], basis1[2], basis1[3]
                site2_x = basis2[1] + dir_x * pv1[1] + dir_y * pv2[1] + dir_z * pv3[1]
                site2_y = basis2[2] + dir_x * pv1[2] + dir_y * pv2[2] + dir_z * pv3[2]
                site2_z = basis2[3] + dir_x * pv1[3] + dir_y * pv2[3] + dir_z * pv3[3]

                max_x = max(max_x, site1_x, site2_x)
                min_x = min(min_x, site1_x, site2_x)
                max_y = max(max_y, site1_y, site2_y)
                min_y = min(min_y, site1_y, site2_y)
                max_z = max(max_z, site1_z, site2_z)
                min_z = min(min_z, site1_z, site2_z)

                bond_palette = colorschemes[:seaborn_bright]
                plot3d!(unit_cell_plot, [site1_x, site2_x], [site1_y, site2_y], [site1_z, site2_z],
                        label="", color=bond_palette[bond.bond_type], lw=2)
        end

        for i in 1:3
                pv = unit_cell.primitive_vectors[:, i]
                min_x = min(min_x, pv[1])
                max_x = max(max_x, pv[1])
                min_y = min(min_y, pv[2])
                max_y = max(max_y, pv[2])
                min_z = min(min_z, pv[3])
                max_z = max(max_z, pv[3])
        end

        for line in unit_cell_lines
                plot3d!(unit_cell_plot, [line[1][1], line[2][1]], [line[1][2], line[2][2]], [line[1][3], line[2][3]],
                        linestyle=:dash, color="black", label="")
        end

        plot3d!(unit_cell_plot,
                xlims=(min_x * 1.1, max_x * 1.1),
                ylims=(min_y * 1.1, max_y * 1.1),
                zlims=(min_z * 1.1, max_z * 1.1),
                aspect_ratio=:equal, xlabel="X", ylabel="Y", zlabel="Z")

        x_positions = getindex.(sites_info, 1)
        y_positions = getindex.(sites_info, 2)
        z_positions = getindex.(sites_info, 3)
        colors = getindex.(sites_info, 4)
        atom_palette = colorschemes[:Accent_8]

        for c in unique(colors)
                mask = colors .== c
                color_index = Int(c)
                plot3d!(unit_cell_plot, x_positions[mask], y_positions[mask], z_positions[mask],
                        seriestype=:scatter3d, markershape=:circle, aspect_ratio=:equal,
                        color=atom_palette[color_index], label="Type $color_index",
                        markerstrokewidth=2, markersize=4)
        end

        display(unit_cell_plot)
        return unit_cell_plot
end

function Lincege.image_lattice(unit_cell::UnitCell; cells=nothing)
        dim = Lincege.dimension(unit_cell)
        try
                if dim == 2
                        return _image_lattice_2d(unit_cell; cells=something(cells, -6:6))
                elseif dim == 3
                        gr()
                        return _image_lattice_3d(unit_cell; cells=something(cells, -2:2))
                else
                        error("Only 2D and 3D unit cells can be visualized. Got dimension: $dim")
                end
        catch e
                @warn "Visualization failed: $e"
                return 1
        end
end

function _image_lattice_2d(unit_cell::UnitCell; cells=-6:6)
        sitesInfo = []

        pv1 = unit_cell.primitive_vectors[:, 1]
        pv2 = unit_cell.primitive_vectors[:, 2]
        basis_1 = unit_cell.basis[:, 1]
        xLim = basis_1[1]+pv1[1]*2
        latticePlot = plot(legend=:best)
        for x in cells
                for y in cells
                        deltaX = pv1[1]*x + pv2[1]*y
                        deltaY = pv1[2]*x+pv2[2]*y
                        for i in 1:Lincege.basis_size(unit_cell)
                                basis_pos = unit_cell.basis[:, i]
                                push!(sitesInfo, [basis_pos[1]+deltaX, basis_pos[2]+deltaY, unit_cell.site_colors[i]])
                        end
        for bond in unit_cell.bonds
                prX = bond.direction[1]
                prY = bond.direction[2]

                basis1 = unit_cell.basis[:, bond.site1]
                basis2 = unit_cell.basis[:, bond.site2]

                site1X = basis1[1] + deltaX
                site1Y = basis1[2] + deltaY
                site2X = basis2[1] + prX * pv1[1] + prY * pv2[1] + deltaX
                site2Y = basis2[2] + prX * pv1[2] + prY * pv2[2] + deltaY

                my_palette = colorschemes[:seaborn_bright]
                plot!(latticePlot, [site1X, site2X], [site1Y, site2Y],
                        label="", color=my_palette[bond.bond_type], lw=2, dpi=500)
        end
end
end

        xPositions = getindex.(sitesInfo, 1)
        yPositions = getindex.(sitesInfo, 2)
        colors = getindex.(sitesInfo, 3)
        paletteAtoms = colorschemes[:Accent_8]

        for c in unique(colors)
                mask = colors .== c
                colorIndex = Int(c)
                plot!(latticePlot, xPositions[mask], yPositions[mask],
                        seriestype=:scatter,
                        aspect_ratio=1,
                        color=paletteAtoms[(colorIndex)],
                        label="Type $colorIndex",
                        markerstrokewidth=2,
                        markersize=6,
                        xlims=(-xLim,xLim),
                        ylims=(-xLim,xLim),
                        dpi=500)
        end

        #display(latticePlot)
        return latticePlot
end
function _image_lattice_3d(unit_cell::UnitCell; cells=-2:2)
        bonds = [
                ExpansionBond([1, bond.site1], [1, bond.site2],
                        bond.direction, bond.bond_type)
                for bond in unit_cell.bonds
        ]
        expansion_cell = ExpansionUnitCell(
                [unit_cell.basis], unit_cell.primitive_vectors,
                [unit_cell.site_colors], [collect(1:Lincege.basis_size(unit_cell))],
                bonds, Bond[])
        return Lincege.image_lattice(expansion_cell; cells)
end

function Lincege.image_unit_cell(
        unit_cell::ExpansionUnitCell; show_expansion_bonds=false)
        dim = Lincege.dimension(unit_cell)
        dim in (2, 3) || throw(ArgumentError("Only 2D and 3D unit cells are supported"))
        dim == 3 && gr()
        return _image_expansion(unit_cell, [zeros(Int, dim)];
                show_expansion_bonds, show_boundary=true)
end

function Lincege.image_lattice(
        unit_cell::ExpansionUnitCell;
        cells=-2:2, show_expansion_bonds=false)
        dim = Lincege.dimension(unit_cell)
        dim in (2, 3) ||
                throw(ArgumentError("Only 2D and 3D lattices are supported"))
        dim == 3 && gr()
        shifts = [collect(shift) for shift in
                Iterators.product(ntuple(_ -> cells, dim)...)]
        isempty(shifts) && throw(ArgumentError("cells must not be empty"))
        return _image_expansion(unit_cell, shifts;
                show_expansion_bonds, show_boundary=false)
end

function _image_expansion(unit_cell, shifts;
        show_expansion_bonds, show_boundary)
        dim = Lincege.dimension(unit_cell)
        positions = Vector{Float64}[]
        site_colors = Int[]
        site_indices = Dict{Tuple,Int}()
        physical_bonds = Set{Tuple{Int,Int,Int}}()

        function add_site(position, color)
                key = (Tuple(round.(position; digits=10)), color)
                return get!(site_indices, key) do
                        push!(positions, position)
                        push!(site_colors, color)
                        length(positions)
                end
        end

        for shift in shifts
                offset = unit_cell.primitive_vectors * shift
                for (element, basis) in enumerate(unit_cell.basis)
                        for site in axes(basis, 2)
                                add_site(basis[:, site] + offset,
                                        unit_cell.site_colors[element][site])
                        end
                end
                for bond in unit_cell.bonds
                        element1, site1 = bond.site1
                        element2, site2 = bond.site2
                        position1 = unit_cell.basis[element1][:, site1] + offset
                        position2 = unit_cell.basis[element2][:, site2] +
                                unit_cell.primitive_vectors * bond.direction + offset
                        index1 = add_site(position1,
                                unit_cell.site_colors[element1][site1])
                        index2 = add_site(position2,
                                unit_cell.site_colors[element2][site2])
                        push!(physical_bonds,
                                (minmax(index1, index2)..., bond.bond_type))
                end
        end

        figure = dim == 2 ?
                plot(legend=:best, aspect_ratio=:equal) :
                plot3d(legend=:outertopright, aspect_ratio=:equal)

        function draw_line(position1, position2; kwargs...)
                coordinates = ([position1[d], position2[d]] for d in 1:dim)
                plot!(figure, coordinates...; kwargs...)
        end

        bond_types = sort!(unique([bond[3] for bond in physical_bonds]))
        bond_palette = palette(:seaborn_bright, max(1, length(bond_types)))
        labeled_types = Set{Int}()
        for (index1, index2, bond_type) in sort!(collect(physical_bonds))
                label = bond_type in labeled_types ? "" : "Bond type $bond_type"
                push!(labeled_types, bond_type)
                draw_line(positions[index1], positions[index2];
                        color=bond_palette[findfirst(==(bond_type), bond_types)],
                        label, lw=2)
        end

        if show_expansion_bonds
                centers = [vec(sum(basis; dims=2)) / size(basis, 2)
                        for basis in unit_cell.basis]
                drawn = Set{Tuple}()
                first_line = true
                for shift in shifts
                        offset = unit_cell.primitive_vectors * shift
                        for bond in unit_cell.expansion_bonds
                                position1 = centers[bond.site1] + offset
                                position2 = centers[bond.site2] +
                                        unit_cell.primitive_vectors * bond.direction +
                                        offset
                                endpoints = minmax(
                                        Tuple(round.(position1; digits=10)),
                                        Tuple(round.(position2; digits=10)))
                                endpoints in drawn && continue
                                push!(drawn, endpoints)
                                draw_line(position1, position2;
                                        linestyle=:dash, color=:gray,
                                        label=first_line ? "Expansion connectivity" : "")
                                first_line = false
                        end
                end
        end

        if show_boundary
                for corner in Iterators.product(ntuple(_ -> 0:1, dim)...)
                        for axis in 1:dim
                                corner[axis] == 0 || continue
                                start = collect(corner)
                                stop = copy(start)
                                stop[axis] = 1
                                draw_line(unit_cell.primitive_vectors * start,
                                        unit_cell.primitive_vectors * stop;
                                        linestyle=:dash, color=:black,
                                        label="")
                        end
                end
        end

        colors = sort!(unique(site_colors))
        atom_palette = palette(:Accent_8, max(1, length(colors)))
        for (index, color) in enumerate(colors)
                selected = findall(==(color), site_colors)
                coordinates = ([positions[i][d] for i in selected] for d in 1:dim)
                scatter!(figure, coordinates...;
                        color=atom_palette[index], label="Type $color",
                        markerstrokewidth=2, markersize=dim == 2 ? 6 : 4)
        end
        return figure
end
end
