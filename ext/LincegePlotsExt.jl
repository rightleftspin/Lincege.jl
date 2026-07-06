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
        first_line = true

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
                        label=first_line ? "Unit Cell" : "")
                first_line = false
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
        first_line = true

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
                        linestyle=:dash, color="black", label=first_line ? "Unit Cell" : "")
                first_line = false
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

end
