module ex8

using NLCE

basis = [[0, 0], [1, 0]]

colors = [1, 2]

primitive_vec = [[1, -1], [1, 1]]

neighborhood = [1]

max_order = 13

nlce_clusters = NLCE.site_color_NLCE(basis, colors, primitive_vec, neighborhood, max_order)

filepath = "examples/outputs/ex-8/ssl_dimer"
mkpath(filepath)
filename = filepath * "/ssl_dimer"

NLCE.write_to_file_colors(nlce_clusters, filename)

end
