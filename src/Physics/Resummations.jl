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

function _eps_wynn(k::Int, n::Int, orders::Vector{Vector{Float64}})
        if k == 0
                n == 0 ? zero(orders[1]) : orders[n]
        elseif k == -1
                zero(orders[1])
        else
                first = _eps_wynn(k - 2, n + 1, orders)
                second = _eps_wynn(k - 1, n + 1, orders) .- _eps_wynn(k - 1, n, orders)
                result = copy(first)
                for i in eachindex(result)
                        abs(second[i]) > 1e-10 && (result[i] += 1.0 / second[i])
                end
                result
        end
end

function _wynn_resum_observable(orders::Vector{Vector{Float64}}, n_cycles::Int)
        final_order = length(orders)
        orders_cumsum = [sum(orders[1:n]) for n in 1:final_order]
        _eps_wynn((2 * n_cycles), final_order - (2 * n_cycles), orders_cumsum)
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
