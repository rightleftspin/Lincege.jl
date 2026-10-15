using Lincege

basis = [[0.0, 0.0]]
primitive_vectors = [[1.0, 0.0], [0.0, 1.0]]
bonds = [Bond(1, 1, [1, 0], 1), Bond(1, 1, [0, 1], 1)]
unit_cell = UnitCell(basis, primitive_vectors, bonds, [1])

max_order = 4
lattice = SiteExpansionLattice(max_order, unit_cell)

symmetries = Vector{Matrix{Float64}}([
        [1 0; 0 1],
        [0 1; -1 0],
        [-1 0; 0 -1],
        [0 -1; 1 0],
        [0 1; 1 0],
        [0 -1; -1 0],
        [1 0; 0 -1],
        [-1 0; 0 1],
])

translation_clusters = TranslationClusterSet(lattice)
clusters_from_lattice!(translation_clusters, lattice)

symmetric_clusters = SymmetricClusterSet(lattice, symmetries)
reduce_clusters!(symmetric_clusters, translation_clusters)

isomorphic_clusters = IsomorphicClusterSet(lattice)
reduce_clusters!(isomorphic_clusters, symmetric_clusters)

symmetric_expansion = Expansion(symmetric_clusters, lattice)
summation!(symmetric_expansion, max_order)

isomorphic_expansion = Expansion(isomorphic_clusters, lattice)
summation!(isomorphic_expansion, max_order)

json_path = joinpath(mktempdir(), "square_symmetries.json")
write_to_json(
        isomorphic_expansion, isomorphic_clusters,
        symmetric_expansion, symmetric_clusters,
        lattice, json_path,
)

txt_path = joinpath(dirname(json_path), "square_symmetries.txt")
write_to_txt_file(
        isomorphic_expansion, isomorphic_clusters,
        symmetric_expansion, symmetric_clusters,
        lattice, txt_path,
)
