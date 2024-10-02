module NLCE

using NautyGraphs

# Add the relevant helper functions
include("helpers/pruning.jl")
include("helpers/util.jl")

# Add the relevant structs
include("nlce_structs/NLCEClusters.jl")
include("nlce_structs/NLCELattices.jl")

# Add the basic pipeline
include("pipeline/grow.jl")
include("pipeline/prune.jl")
include("pipeline/combine.jl")

# Add the relevant helper functions
include("helpers/generation.jl")


# Initializing wrappers for Directed (Di), Edge Weighted (Ew),
# and Vertex Labelled (Vl) lattices
const DiEwVlLattice = NLCELattice{true,true,true}

const DiLattice = NLCELattice{true,false,false}
const DiEwLattice = NLCELattice{true,true,false}

const EwLattice = NLCELattice{false,true,false}
const EwVlLattice = NLCELattice{false,true,true}

const VlLattice = NLCELattice{false,false,true}
const DiVlLattice = NLCELattice{true,true,false}

const Lattice = NLCELattice{false,false,false}


function test_stuff()
    basis = [[0, 0]]
    primitive_vec = [[1, 0], [0, 1]]
    neighborhood = [1]
    max_order = 17
    msq = max_order ^ 2
    
    lattice = generate_lattice(basis, primitive_vec, neighborhood, max_order)
    
    grow(lattice, max_order - 6, [div(msq, 2)])

    #clus::Vector{Int64} = rand(1:msq, max_order)
    #for i in 1:1000
    #    ind_clus = cluster(lattice, clus)
    #end
    
  #  println(ind_clus.adj_list)
  #  println(ind_clus.adj_matrix)
  #  for (l, m) in zip(ind_clus.adj_list_weights, ind_clus.adj_matrix_weights)
  #      println(l)
  #      println(m)
  #  end
end

using BenchmarkTools

@time test_stuff()

export 
    isomorphic_tagging
    generate_lattice
     
# Below is an example of the NLCE process for a square lattice with only nearest neighbors
# up till order 4 it uses a wrapper function for speed, check the helper functions for 
# more info

# add write to file helper
#using AlgebraicNumbers: AlgebraicNumber as AN
#using Plots; gr()
#
### Tasks before next meeting
## do square, triangle and kagome lattice ising models
#
### Setting up the lattice geometry
#basis = [[AN(0), AN(0)], [AN(1), AN(0)], [AN(1)/2, sqrt(AN(3))/2]] 
#primitive_vectors = [[AN(2), AN(0)], [AN(1), sqrt(AN(3))]]
#final_order = 7
#
#for order in 1:final_order
#    println(length(simple_nlce(basis, primitive_vectors, order)))
#end

#xmax = 2
## Initialize the temperature grid
#total_temp = collect(0.05:0.05:xmax)
#all_energies = []
#
## Loop over each desired order
#for order in 5:final_order
#    # Initialize the energies
#    total_energy = zeros(length(total_temp))
#    # Perform the NLCE sum and find resulting weights
#    nlce_weights = simple_nlce(basis, primitive_vectors, order)
#    for (cluster, mult) in nlce_weights
#        # Find the energies for each cluster and add it to the total sum
#        ising_energy = ising_energies(cluster)
#        total_energy += mult .* energy_solver(total_temp, ising_energy) 
#    end
#    push!(all_energies, total_energy)
#end
#
## Plot
#plot(total_temp, all_energies, label = collect(7:final_order))
#xlims!(0, xmax)
#ylims!(-1, 0)
#title!("Ising Energy for a triangular Lattice")
#xlabel!("Temperature")
#ylabel!("Energy")
#savefig("ising_energy.png")
#
end
