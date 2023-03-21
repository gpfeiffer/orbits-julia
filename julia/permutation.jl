#############################################################################
##
module permutation

import Base: length, hash, inv, isless, one, sign, rand, iterate, ==, *, /, ^

export Perm
export degree, domain, cycles, shape, order, transposition, shuffle!
export permuted, is_trivial, largest_moved_point
export p, q, r, transpositions, ttt  # test data

##  Perm data type
struct Perm
    list::Array{Int}
end

## identity
Perm(n::Int) = Perm(1:n)
one(Perm, n::Int) = Perm(1:n)
one(perm::Perm) = Perm(domain(perm))

## degree
degree(perm::Perm) = length(perm.list)

## domain
domain(perm::Perm) = 1:degree(perm)

## equality and hash
==(perm::Perm, other::Perm) = perm.list == other.list
hash(perm::Perm, h::UInt) = hash(perm.list, hash(Perm, h))

## comparison for size
isless(perm::Perm, other::Perm) = perm.list < other.list

## point under perm
^(x::Int, perm::Perm) = perm.list[x]

## cycles ('DFS' version:  each node i has two children: i^perm and i+1)
function cycles(perm::Perm, x = 1, open = trues(degree(perm)), cycle = [], ccc = Array{Int}[])
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

function Perm(n::Int, cycles::Array)
    perm = Perm(n)
    for c in cycles
        perm.list[c] = c[(1:end) .% end .+ 1]
    end
    return perm
end

## inverse permutation
function inv(perm::Perm)
    list = similar(perm.list)
    list[perm.list] = domain(perm)
    Perm(list)
end

## product and quotient
*(perm::Perm, other::Perm) = Perm(other.list[perm.list])
/(perm::Perm, other::Perm) = perm * inv(other)

## power
function ^(perm::Perm, n::Int)
    n < 0 && return inv(perm)^(-n)
    n == 0 && return one(perm)
    n == 1 && return perm

    q, r = divrem(n, 2)
    return perm^q * perm^q * perm^r
end

## conjugation
^(perm::Perm, other::Perm) = other^(-1) * perm * other

## permuted
permuted(list::Array, perm::Perm) = list[inv(perm).list]

## shape (aka cycle structure)
shape(perm::Perm) = sort(length.(cycles(perm)), rev=true)

## order
order(perm::Perm) = lcm(shape(perm))

## sign
sign(perm::Perm) = (-1)^(degree(perm) - length(cycles(perm)))

## transpositions
function transposition(n::Int, j::Int, k::Int)
    list = collect(1:n)
    list[[j, k]] = [k, j]
    return Perm(list)
end

## shuffle:  Fisher-Yates
function shuffle!(perm::Perm)
    n = degree(perm)
    for k in 1:n-1
        l = rand(k:n);
        perm.list[[l, k]] = perm.list[[k, l]]
    end
    return perm
end

rand(T::Type{X}, n::Int) where X <: Perm = shuffle!(one(Perm, n))

##  iterate (copied from number.jl)
iterate(x::Perm) = (x, nothing)
iterate(x::Perm, ::Any) = nothing

## is trivial?
is_trivial(perm::Perm) = perm.list == domain(perm)

## largest moved point
function largest_moved_point(perm::Perm)
    d = degree(perm)
    while d^perm == d
        d -= 1
    end
    return d
end

## test data
transpositions(n::Int) = [transposition(n, j-1, j) for j in 2:n]
ttt = transpositions(20)
p = Perm([4,2,3,1])
q = Perm([2,3,4,1])
r = rand(Perm, 20)

end # module
