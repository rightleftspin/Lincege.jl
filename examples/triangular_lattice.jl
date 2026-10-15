using Lincege

triangular_basis = [[0.0, 0.0]]
triangular_pvecs = [[1.0, 0.0], [1.0/2, sqrt(3)/2]]
triangular_bonds = [Bond(1,1,[1,0],1),Bond(1,1,[0,1],1),Bond(1,1,[-1,1],1)]
triangular_uc = UnitCell(triangular_basis, triangular_pvecs, triangular_bonds, [1])

using Plots, ColorSchemes
image_unit_cell(triangular_uc)

image_lattice(triangular_uc)

m_order = 3
lattice = SiteExpansionLattice(m_order, triangular_uc)

trans_clusters = TranslationClusterSet(lattice)
clusters_from_lattice!(trans_clusters, lattice)

iso_clusters = IsomorphicClusterSet(lattice)
reduce_clusters!(iso_clusters, trans_clusters)

expansion = Expansion(iso_clusters, lattice)
summation!(expansion, m_order)

write_to_json(expansion, lattice, "triangular_lattice.json")
