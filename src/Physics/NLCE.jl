"""
    NLCEResult

Stores the results of an NLCE computation, can be plotted by accessing temperatures and the corresponding observable.
"""
struct NLCEResult
        temperatures::Vector{Float64}
        energy::Vector{Vector{Float64}}
        specific_heat::Vector{Vector{Float64}}
        entropy::Vector{Vector{Float64}}
        n_orders::Int
end

"""
    perform_nlce(solvers, temperatures) -> NLCEResult

Calls the given NLCE solvers and computes their corresponding observables.
"""
function perform_nlce(solvers::Vector{<:AbstractPhysicsSolver}, temperatures::AbstractVector{Float64})
        betas = 1.0 ./ temperatures
        n_T = length(temperatures)
        max_order = maximum(length(cluster_weights(s)) for s in solvers)

        energy = [zeros(n_T) for _ in 1:max_order]
        specific_heat = [zeros(n_T) for _ in 1:max_order]
        entropy = [zeros(n_T) for _ in 1:max_order]

        for solver in solvers
                wts = cluster_weights(solver)
                obs = observables(solver, betas)
                for (order, w) in enumerate(wts)
                        iszero(w) && continue
                        energy[order] .+= w .* obs.energy
                        specific_heat[order] .+= w .* obs.specific_heat
                        entropy[order] .+= w .* obs.entropy
                end
        end

        NLCEResult(collect(temperatures), energy, specific_heat, entropy, max_order)
end
