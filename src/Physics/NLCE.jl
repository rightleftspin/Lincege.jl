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

function _bincoeff(n::Int, k::Int)
        k > n && return 0.0
        r = 1.0
        for d in 1:k
                r = r * n / d
                n -= 1
        end
        r
end

function _euler_resum_observable(orders::Vector{Vector{Float64}}, sum_start::Int)
        n_T = length(orders[1])
        baseline = sum_start > 1 ? reduce(.+, orders[1:sum_start-1]) : zeros(n_T)
        partial = orders[sum_start:end]
        out = zeros(n_T)
        for i in 1:length(partial)
                delta = zeros(n_T)
                for j in 1:i
                        delta .+= (-1)^j * _bincoeff(i, j) .* abs.(partial[j])
                end
                out .+= (0.5^(i + 1)) .* delta
        end
        baseline .+ out
end

function _wynn_eps!(memo::Dict{Tuple{Int,Int},Vector{Float64}}, ps::Vector{Vector{Float64}}, n_T::Int, k::Int, n::Int)
        haskey(memo, (k, n)) && return memo[(k, n)]
        val = if k == 0
                n == 0 ? copy(ps[end]) : copy(ps[n])
        elseif k == -1
                zeros(n_T)
        else
                f = _wynn_eps!(memo, ps, n_T, k - 2, n + 1)
                d = _wynn_eps!(memo, ps, n_T, k - 1, n + 1) .- _wynn_eps!(memo, ps, n_T, k - 1, n)
                f .+ 1.0 ./ d
        end
        memo[(k, n)] = val
        val
end

function _wynn_resum_observable(orders::Vector{Vector{Float64}}, n_cycles::Int)
        n_T = length(orders[1])
        n_orders = length(orders)
        ps = Vector{Vector{Float64}}(undef, n_orders)
        ps[1] = copy(orders[1])
        for k in 2:n_orders
                ps[k] = ps[k-1] .+ orders[k]
        end
        n0 = n_orders - 2 * n_cycles
        _wynn_eps!(Dict{Tuple{Int,Int},Vector{Float64}}(), ps, n_T, 2 * n_cycles, n0)
end

"""
    euler_resummation(result, sum_start)

Performs euler's resummation technique for alternating series to an NLCE result. sum_start is the order index at which the Euler sum begins; the early orders (below sum_start) are left as bare sums.
"""
function euler_resummation(result::NLCEResult, sum_start::Int=1)
        @assert sum_start >= 1 && sum_start <= result.n_orders "sum_start ($sum_start) must be between 1 and n_orders ($(result.n_orders))"
        NLCEResult(
                result.temperatures,
                vcat(result.energy, [_euler_resum_observable(result.energy[1:result.n_orders], sum_start)]),
                vcat(result.specific_heat, [_euler_resum_observable(result.specific_heat[1:result.n_orders], sum_start)]),
                vcat(result.entropy, [_euler_resum_observable(result.entropy[1:result.n_orders], sum_start)]),
                result.n_orders,
        )
end

"""
    wynn_resummation(result, n_cycles)

Performs wynn's resummation technique based on pade approximants to an NLCE result. The number of cycles of improvement should generally be maximized, this is the default.
"""
function wynn_resummation(result::NLCEResult, n_cycles::Int=fld(result.n_orders, 2))
        @assert 2 * n_cycles <= result.n_orders "2 * n_cycles ($(2 * n_cycles)) exceeds the number of orders ($(result.n_orders))"
        NLCEResult(
                result.temperatures,
                vcat(result.energy, [_wynn_resum_observable(result.energy[1:result.n_orders], n_cycles)]),
                vcat(result.specific_heat, [_wynn_resum_observable(result.specific_heat[1:result.n_orders], n_cycles)]),
                vcat(result.entropy, [_wynn_resum_observable(result.entropy[1:result.n_orders], n_cycles)]),
                result.n_orders,
        )
end

"""
    apply_resummations(result, resummations)

Applies resummations to a given NLCEResult in order of the given resummation vector. The vector is of the form [(:Euler, 3), (:Wynn, 4), ...] which translates to euler summation starting at order 3, and wynn's resummation with 4 cycles of improvement.
"""
function apply_resummations(result::NLCEResult, resummations::Vector{Tuple{Symbol,Int}})
        for (method, param) in resummations
                if method == :Euler
                        result = euler_resummation(result, param)
                elseif method == :Wynn
                        result = wynn_resummation(result, param)
                else
                        error("Unknown resummation method :$method; expected :Euler or :Wynn")
                end
        end
        result
end
