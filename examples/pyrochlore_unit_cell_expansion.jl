using Lincege

pyro_exp_uc_basis = [
        [
                [1 / 2, 1 / 2, 1 / 2],
                [1 / 2, -1 / 2, -1 / 2],
                [-1 / 2, 1 / 2, -1 / 2],
                [-1 / 2, -1 / 2, 1 / 2],
        ],
]
pyro_exp_uc_pvecs = [[2.0, 2.0, 0.0], [2.0, 0.0, 2.0], [0.0, 2.0, 2.0]]
pyro_exp_uc_lbonds = [
        # Inside a Unit Cell
        ExpansionBond([1, 1], [1, 2], [0, 0, 0], 1),
        ExpansionBond([1, 1], [1, 3], [0, 0, 0], 1),
        ExpansionBond([1, 1], [1, 4], [0, 0, 0], 1),
        ExpansionBond([1, 2], [1, 3], [0, 0, 0], 1),
        ExpansionBond([1, 2], [1, 4], [0, 0, 0], 1),
        ExpansionBond([1, 3], [1, 4], [0, 0, 0], 1),
        # Between Unit Cells
        ExpansionBond([1, 1], [1, 2], [0, 0, 1], 1),
        ExpansionBond([1, 1], [1, 3], [0, 1, 0], 1),
        ExpansionBond([1, 1], [1, 4], [1, 0, 0], 1),
        ExpansionBond([1, 2], [1, 3], [0, 1, -1], 1),
        ExpansionBond([1, 2], [1, 4], [1, 0, -1], 1),
        ExpansionBond([1, 3], [1, 4], [1, -1, 0], 1),
]
pyro_exp_uc_ebonds = [
        Bond(1, 1, [1, 0, 0], 1),
        Bond(1, 1, [0, 1, 0], 1),
        Bond(1, 1, [0, 0, 1], 1),
        Bond(1, 1, [0, 1, -1], 1),
        Bond(1, 1, [1, 0, -1], 1),
        Bond(1, 1, [1, -1, 0], 1),
]
pyro_exp_uc_uc = ExpansionUnitCell(
        pyro_exp_uc_basis,
        pyro_exp_uc_pvecs,
        pyro_exp_uc_lbonds,
        pyro_exp_uc_ebonds,
        [[1, 1, 1, 1]]
)

using Plots, ColorSchemes
image_unit_cell(pyro_exp_uc_uc)

image_lattice(pyro_exp_uc_uc; cells=0:1, show_expansion_bonds=true)

m_order = 3
lattice = StrongClusterExpansionLattice(m_order, pyro_exp_uc_uc)

trans_clusters = TranslationClusterSet(lattice)
clusters_from_lattice!(trans_clusters, lattice)

iso_clusters = IsomorphicClusterSet(lattice)
reduce_clusters!(iso_clusters, trans_clusters)

expansion = Expansion(iso_clusters, lattice)
summation!(expansion, m_order)

write_to_json(expansion, lattice, "pyrochlore_lattice_uc_exp.json")
