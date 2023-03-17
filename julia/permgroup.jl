#############################################################################
##
module permgroup

using permutation
using orbits

import Base: in, ==

export PermGp, elements, closure, subgroups

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
is_subgp(other::PermGp, gp::PermGp) = all(a -> a in other, gp.gens)

## equality
==(gp::PermGp, other::PermGp) = is_subgp(other, gp) &&  is_subgp(gp, other)

##  closure
closure(group::PermGp, a::Perm) = PermGp(union(group.gens, a), group.one)
onGroups(x, a) = closure(x, a)

##  subgroups
subgroups(group) = orbit(elements(group), PermGp([], group.one), onGroups)


end # module
