#############################################################################
##
module orbits

export orbit, onPoints, onRight

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

end # module
