#############################################################################
##
module orbits

import Base: in, isless, size, ==

export Orbit
export orbit, orbitx, onPoints, onRight, onPairs, onSets
export orbit_with_words, orbit_with_transversal, orbit_with_stabilizer
export orbit_with_edges

## orbit
function orbit(aaa, x, under)
    list = [x]
    for y in list
        for a in aaa
            z = under(y, a)
            z in list || push!(list, z)
        end
    end
    return list
end

## on points
onPoints(x, a) = x^a

## on right
onRight(x, a) = x * a

## on pairs
onPairs(pair::Pair, a::Perm) = Pair(pair.first^a, pair.second^a)

## on sets
onSets(set::Set, a ::Perm) = Set([x^a for x in set])

##  Orbit data type
struct Orbit
    group
    elts # assume sorted!
end

## size
size(o::Orbit) = length(o.elts)

## membership, equality, comparison
in(x, o::Orbit) = x in o.elts
==(o::Orbit, other::Orbit) = o.elts[1] == other.elts[1]
isless(o::Orbit, other::Orbit) = o.elts[1] < other.elts[1]

## orbit with words
function orbit_with_words(aaa, x, under)
    list = [x]
    words = Array{Int}[[]]
    i = 0
    while i < length(list)
        i += 1
        for k in 1:length(aaa)
            z = under(list[i], aaa[k])
            if !(z in list)
                push!(list, z)
                push!(words, onWords(words[i], k))
            end
        end
    end
    return Dict(:list => list, :words => words)
end

## orbitx
function orbitx(aaa, xxx, under)
    list = copy(xxx)
    for y in list
        for a in aaa
            z = under(y, a)
            z in list || push!(list, z)
        end
    end
    return list
end

## orbit with transversal
function orbit_with_transversal(aaa, x, under)
    list = [x]
    reps = [aaa[1]^0]
    i = 0
    while i < length(list)
        i += 1
        for k in 1:length(aaa)
            z = under(list[i], aaa[k])
            if !(z in list)
                push!(list, z)
                push!(reps, reps[i] * aaa[k])
            end
        end
    end
    return Dict(:list => list, :reps => reps)
end

## orbit with stabilizer
function orbit_with_stabilizer(aaa, x, under)
    list = [x]
    reps = [aaa[1]^0]
    stab = []
    i = 0
    while i < length(list)
        i += 1
        for k in 1:length(aaa)
            z = under(list[i], aaa[k])
            l = findfirst(==(z), list)
            if l == nothing
                push!(list, z)
                push!(reps, reps[i] * aaa[k])
            else   # x^(reps[i] * a) = x^reps[l]
                push!(stab, reps[i] * aaa[k] / reps[l])
            end
        end
    end
    return Dict(:list => list, :reps => reps, :stab => stab)
end

## orbit with edges
function orbit_with_edges(aaa, x, under)
    list = [x]
    edges = []
    i = 0
    while i < length(list)
        i += 1
        for k in 1:length(aaa)
            z = under(list[i], aaa[k])
            l = findfirst(==(z), list)
            if l == nothing
                push!(list, z)
                l = length(list)
            end
            push!(edges, [i, l])
        end
    end
    return Dict(:list => list, :edges => edges)
end

end # module
