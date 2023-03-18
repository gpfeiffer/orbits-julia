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
