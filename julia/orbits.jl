#############################################################################
##
module orbits

import Base: in, isless, size, ==

export Orbit
export orbit, orbitx, onPoints, onRight

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

end # module
