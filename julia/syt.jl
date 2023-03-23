#############################################################################
##
#A  syt.jl                                                       orbits-julia
#B    by Götz Pfeiffer <goetz.pfeiffer@nuigalway.ie>
##
#C  Partitions and Standard Young Tableaus as orbits under suitable actions
##
module syt

export newtonDif, newtonSum, newtonDifR, newtonSumR, orbit
export takeAway, subsets, partitions, standardYTs, tableau_path

#############################################################################
##
##  newtonSum(list) and newtonDif(list)
##
##  Newton's partial sums and forward differences of a list
##
##  properties:
##    newtonDif(newtonSum(list)) = list
##    newtonSum(newtonDif(list)) = list
##

## [l_1, l_2 ..  l_n] -> [l_1, l_2 - l_1, l_3 - l_2 .. l_n -l_{n-1}]
function newtonDif(list)
    list == [] && return []
    vcat(list[1], newtonDif(list[2:end] .- list[1]))
end

## [s_1, s_2 ..  s_n] -> [s_1, s_1 + s_2, s_1 + s_2 + s_3 .. s_1 + ... + s_n]
function newtonSum(list)
    list == [] && return []
    vcat(list[1], newtonSum(list[2:end]) .+ list[1])
end

##  right handed versions:
##    newtonDifR(newtonSumR(list)) = list
##    newtonSumR(newtonDifR(list)) = list
##

## [l_1, l_2 ..  l_n] -> [l_2 - l_1 .. l_n - l_{n-1}, l_n]
function newtonDifR(list)
    list == [] && return []
    vcat(newtonDifR(list[1:end-1] .- list[end]), list[end])
end

## [s_1, s_2 ..  s_n] -> [s_1 + ... + s_n, s_2 + ... + s_n .. s_n]
function newtonSumR(list)
    list == [] && return []
    vcat(newtonSumR(list[1:end-1]) .+ list[end], list[end])
end

#############################################################################
##
##  Compositions (of n) vs Subsets (of {1..n-1})
##
##  how to turn a subset of [1..n-1] into a composition of n:
##  * find the complement cmp of L in [1..n]
##  * compute Newton differences of cmp
##  e.g. n = 9, L =    45 7  \subseteq [1..8],
##            cmp = 123  6 89,
##            com = 111  3 21.
##
compositionSubset(n, set) = newtonDif(setdiff(1:n, set))

##  how to turn a composition of n into a subset of [1..n-1]
subsetComposition(n, com) = setdiff(1:n, newtonSum(com))

##  examples of orbits:

##  the usual orbit algorithm, again.
function orbit(aaa, x, under)
    list = [x]
    for y in list
        for a in aaa
            z = under(y, a)
            if z != y && !(z in list)
                push!(list, z)
            end
        end
    end
    return list
end

##  take-away action
takeAway(set, s) = setdiff(set, s)

##  power set as orbit under take-away
subsets(set) = orbit(set, set, takeAway)

##  partitions of n as takeAway orbit of composition classes
function partitions(n)
    function takeAwayPartition(com, a)
        com = compositionSubset(n, takeAway(subsetComposition(n, com), a))
        sort(com, by=-)
    end
    reverse(orbit(1:n-1, [n], takeAwayPartition))
end

#############################################################################
##
##  SYT: a standard Young tableau is a shortest path in the Young lattice
##

##  normal case action
function remove1Hook2(lambda, k)
    new = copy(lambda)
    new[k] -= 1
    return new
end

##  special case action
function remove1Hook1(lambda, k)
    lambda[k] == 1 && return lambda[1:k-1]
    remove1Hook2(lambda, k)
end

##  Young lattice and shortest paths.  2-step process.
##  1. build graph "bottom up" from lambda; for each edge, record the row/col
##     positions it corresponds to.
##  2. find shortest paths "top down" from 0
##  (consider setting this up as starting with a list of partitions of n)
##
function youngLattice(lambda)

    ## how to add to the list
    function grow(list, y, r, c)
        pos = findfirst(==(y), list)
        if pos == nothing
            push!(list, y)
            push!(next, [])
            pos = length(list)
        end
        push!(next[i], (pos = pos, row = r, col = c))
    end

    ## orbit-with-edges
    list = [lambda]
    next = [[] for x in list]
    i = 0
    while i < length(list)
        i += 1
        x = list[i]
        l = length(x)
        if l > 0
            dif = newtonDifR(x)   # ;-)
            grow(list, remove1Hook1(x, l), l, x[l])  # dif[l] > 0
            for k in reverse(1:l-1)
                dif[k] > 0 && grow(list, remove1Hook2(x, k), k, x[k])
            end
        end
    end
    return (list = list, next = next)
end

##  DFS with a visitor that makes shortest paths if unknown.
function pathfinder(x, next, path)
    if ismissing(path[x])
        path[x] = []
        for z in next[x]
            for p in pathfinder(z.pos, next, path)
                push!(path[x], vcat(p, [(z.row, z.col)]))
            end
        end
        path[x] == [] && push!(path[x], [])   # empty path
    end
    return path[x]
end

## list of SYTs of shape lambda
function standardYTs(lambda)
    next = youngLattice(lambda).next
    return pathfinder(1, next, Any[missing for x in next])
end

##  tableau from path
function tableau_path(path)
    tab = Array{Int}[]
    for i in 1:length(path)
        r, c = path[i]
        c == 1 && push!(tab, [])
        push!(tab[r], i)
    end
    return tab
end

end # module
