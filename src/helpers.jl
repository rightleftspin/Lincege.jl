using AlgebraicNumbers

function Base.isless(a::AlgebraicNumber, b::AlgebraicNumber)
    return (real(a.apprx) < real(b.apprx)) 
end

"""
    find_bonds(coordinates, neighborhood)

Given a set of coordinates, this function returns a bond array linking each
coordinate with the corresponding coordinates that are within a specified 
neighborhood distance from the coordinate, listed by index

# Examples
```julia-repl
julia> find_bonds([[0, 0], [1, 1], [0, 1], [1, 0]], 2)

```
"""

function find_bonds(coordinates::Vector{Vector{T}}, neighbors::Int) where {T}
    magnitude = coord -> sqrt(sum(coord.^2))

    bonds = [Vector{Int}(undef, 0) for i in 1:length(coordinates)]
    bond_length = [Vector{Int}(undef, 0) for i in 1:length(coordinates)]
    bond_direction = [Vector{Int}(undef, 0) for i in 1:length(coordinates)]

    directions = Dict{Vector{T}, Int}()
    
    # Find the center of all the coordinates
    center = coordinates[argmin(magnitude.(coordinates .- Ref(sum(coordinates)/length(coordinates))))]
    # Find all the minimum distances, ignore 0, which there should only be one of,
    # and then retain only the distances you care about
    distances = sort(magnitude.(coordinates .- Ref(center)))[2:end]
    neighborhood = []
    for distance in distances
        if !(distance in neighborhood)
            push!(neighborhood, distance)
        end
        if length(neighborhood) == neighbors
            break
        end      
    end    

    # Add all the relevant bonds
    for (i, coordinate) in enumerate(coordinates)
        for (j, neighbor) in enumerate(coordinates)
            for (k, neighborhood_dist) in enumerate(neighborhood)
                direction = coordinate - neighbor
                if magnitude(direction) == neighborhood_dist 
                    push!(bonds[i], j)
                    push!(bond_length[i], k)
                    if direction in keys(directions)
                        push!(bond_direction[i], directions[direction])
                    else
                        directions[direction] = length(directions) + 1
                        push!(bond_direction[i], directions[direction])
                    end
                end
            end
        end
    end
    (bonds, bond_length, bond_direction)
end

#x = [[AlgebraicNumber(a), AlgebraicNumber(b)] for a in 0:2 for b in 0:2]
#x = [[Float64(a), Float64(b)] for a in 0:2 for b in 0:2]
#println(x)
#println(find_bonds(x, 2)[1])
#println(find_bonds(x, 2)[2])
