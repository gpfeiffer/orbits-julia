#############################################################################
##
module permgroup

using permutation
using orbits

import Base: in, <=, ==

export PermGp, elements, closure

## Perm group
struct PermGp
    gens::Array{Perm}
    one::Perm
end

## elements
elements(group::PermGp) = sort(orbit(group.gens, group.one, onRight))

## membership
in(a::Perm, group::PermGp) = a in elements(group)

## subgroup
<=(other::PermGp, group::PermGp) = all(a -> a in group, other.gens)

## equality
==(group::PermGp, other::PermGp) = other <= group && group <= other

##  closure
closure(group::PermGp, a::Perm) = PermGp(union(group.gens, a), group.one)

end # module
