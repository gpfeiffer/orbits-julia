#############################################################################
##
module permutation

import Base: length, inv, isless, one, ==, *, ^

export Perm
export degree, domain, cycles, transposition
export p, q  # test data

##  Perm data type
struct Perm
    list::Array{Integer}
end

## identity
one(perm::Perm) = Perm(domain(perm))
one(Perm, n::Integer) = Perm(1:n)

## degree
degree(perm::Perm) = length(perm.list)

## domain
domain(perm::Perm) = 1:degree(perm)

## equality
==(perm::Perm, other::Perm) = perm.list == other.list

## comparison for size
isless(perm::Perm, other::Perm) = perm.list < other.list

## point under perm
^(x::Integer, perm::Perm) = perm.list[x]

## cycles ('DFS' version:  each node i has two children: i^perm and i+1)
function cycles(perm::Perm, x = 1, open = trues(degree(perm)), cycle = [], ccc = Array{Integer}[])
    open != [] || return ccc
    if open[x]
        open[x] = false
        push!(cycle, x) # continue cycle with x
        cycles(perm, x^perm, open, cycle, ccc)
    else
        cycle == [] || push!(ccc, cycle)
        x < length(open) && cycles(perm, x+1, open, [], ccc)
    end
    return ccc
end

## inverse permutation
function inv(perm::Perm)
    list = similar(perm.list)
    list[perm.list] = domain(perm)
    Perm(list)
end

## product
*(perm::Perm, other::Perm) = Perm(other.list[perm.list])

## power
function ^(perm::Perm, n::Integer)
    n < 0 && return inv(perm)^(-n)
    n == 0 && return one(perm)
    n == 1 && return perm

    q, r = divrem(n, 2)
    return perm^q * perm^q * perm^r
end

## conjugates
^(perm::Perm, other::Perm) = other^(-1) * perm * other

## sign??

## order??

## cycles??

## transpositions
function transposition(n::Integer, j::Integer, k::Integer)
    list = collect(1:n)
    list[[j, k]] = [k, j]
    return Perm(list)
end

#transpositions(n::Integer) = [transposition(n, j-1, j) for j in 2:n]


## test data
p = Perm([4,2,3,1])
q = Perm([2,3,4,1])

end # module
