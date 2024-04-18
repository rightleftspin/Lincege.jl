"""
    Cell{
        Coordinate<:AlgebraicNumber,
        Coloring<:Integer
    }

A type to capture a single instance of a unit cell of a lattice, including constituent coordinates/vertices, internal bonds
and the coordinate/vertex colorings

# Fields
- `coordinates::Vector{Vector{Coordinate}}`: vector of all the coordinates
- `coloring::Vector{Coloring}`: coloring for each vertex, in the same order as vertices
- `internal_bonds::Vector{Vector{Integer}}`: vector which contains vectors of all bonds within the cell for each coordinate in coordinates, by index
- `external_bonds::Vector{Vector{Cell}}`: vector which contains vectors of all bonds outside the cell for each coordinate in coordinates, by index
"""
struct Cell{
            Coordinate<:AlgebraicNumber,
            Coloring<:Integer
           }
    coordinates::Vector{Vector{Coordinate}}, # Replace with static array internally since they are coordinates
    coloring::Vector{Coloring},
    internal_bonds::Vector{Vector{Integer}}
    external_bonds::Vector{Vector{Integer}}
end

## Constructors
function Cell(coordinates::Vector{Vector{AlgebraicNumber}})
    Cell(
        coordinates,
        ones(size(coordinates)[1]),
        Vector{Vector{Integer}}(undef, size(coordinates)[1]),
        Vector{Vector{Integer}}(undef, size(coordinates)[1])
        )
end

function Cell(coordinates::Vector{Vector{AlgebraicNumber}}, coloring::Vector{Integer})
    Cell(
        coordinates,
        coloring,
        Vector{Vector{Integer}}(undef, size(coordinates)[1])
        )
end

## Base Extensions
function Base.show(io::IO, cell::Cell)
    print(
        io,
        "Cell with $(nv(cell)) coordinates and $(ne(cell)) internal bonds",
    )
    return nothing
end

## Methods
function colored_coordinates(cell::Cell)
    return(Dict(cell.coordinates .=> cell.coloring))    
end     

function neighbors(cell::Cell{Coordinate, Any}, coordinate::Vector{Coordinate}) where {Coordinate}
    for i, cell_coordinate in enumerate(cell.coordinates)
        if coordinate == cell_coordinate
            return(internal_bonds[i], [cell.coordinates[c] for c in internal_bonds[i]])
        end
    end
    return Nothing
end

function neighbors(cell::Cell, coordinate::Integer)
    if coordinate <= nv(cell)
        return(internal_bonds[coordinate], [cell.coordinates[c] for c in internal_bonds[coordinate]])
    end
    return Nothing
end

function nv(cell::Cell)
    return(size(cell.coordinates)[1])
end

function ne(cell::Cell)
    bond_total = 0
    for vertex, bonds in enumerate(cell.internal_bonds)
        for bond in bonds
            if bond > vertex
                bond_total += 1
            end
        end
    end
    return(bond_total)
end
