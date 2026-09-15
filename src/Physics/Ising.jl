"""
    IsingSolver <: AbstractPhysicsSolver

Concrete implementation of a physics solver for the Ising model.
"""
struct IsingSolver <: AbstractPhysicsSolver
        eigenvalues::Vector{Float64}
        eigenvectors::Matrix{Float64}
        weights::Vector{Float64}
end

function IsingSolver(cluster::ImportedCluster; J::Float64=1.0)
        n = cluster.n_sites
        n_configs = 1 << n
        energy_eigenvalues = Vector{Float64}(undef, n_configs)
        energy_eigenvectors = Matrix{Float64}(undef, n, n_configs)

        for config in 0:(n_configs-1)
                spins = [((config >> (k - 1)) & 1) == 1 ? 1.0 : -1.0 for k in 1:n]
                E = sum(b -> -J * spins[b[1]] * spins[b[2]], cluster.bonds; init=0.0)
                energy_eigenvalues[config+1] = E
                energy_eigenvectors[:, config+1] = spins
        end

        IsingSolver(energy_eigenvalues, energy_eigenvectors, cluster.weights)
end

eigenvalues(s::IsingSolver) = s.eigenvalues
eigenvectors(s::IsingSolver) = s.eigenvectors
cluster_weights(s::IsingSolver) = s.weights
