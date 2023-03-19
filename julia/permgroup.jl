#############################################################################
##
module permgroup

using permutation
using orbits

import Base: in, isless, iterate, size, ==, ^
import permutation: is_trivial, largest_moved_point

export PermGp, elements, conjClasses, closure, subgroups, subgpClasses

## Perm group
struct PermGp
    gens::Array{Perm}
    one::Perm
end

## elements
elements(group::PermGp) = sort(orbit(group.gens, group.one, onRight))

## size
size(group::PermGp) = length(elements(group))

## membership
in(a::Perm, group::PermGp) = a in elements(group)

## equality, comparison
==(group::PermGp, other::PermGp) = elements(group) == elements(other)
isless(group::PermGp, other::PermGp) = elements(group) < elements(other)

##  closure
closure(group::PermGp, a::Perm) = PermGp(union(group.gens, a), group.one)
onGroups(x, a) = closure(x, a)

##  subgroups
subgroups(group) = orbit(elements(group), PermGp([], group.one), onGroups)

## conjugation
class(group::PermGp, a::Perm) = orbit(group.gens, a, onPoints)
^(a::Perm, group::PermGp) = Orbit(group, sort(class(group, a)))

^(group::PermGp, a::Perm) = PermGp([x^a for x in group.gens], group.one)

subgpClass(gp::PermGp, subgp::PermGp) = orbit(gp.gens, subgp, onPoints)
^(subgp::PermGp, group::PermGp) = Orbit(group, sort(subgpClass(group, subgp)))

##  iterate (copied from number.jl)
iterate(x::PermGp) = (x, nothing)
iterate(x::PermGp, ::Any) = nothing

## conjugacy classes
onClasses(x, a) = (x.elts[1] * a)^(x.group)

function conjClasses(gp::PermGp)
    orbit(orbitx(gp.gens, gp.gens, onPoints), gp.one^gp, onClasses)
end

onSubgpClasses(x, a) = onGroups(x.elts[1], a)^(x.group)

function subgpClasses(gp::PermGp)
    orbit(elements(gp), PermGp([], gp.one)^gp, onSubgpClasses)
end

## is trivial?
is_trivial(group::PermGp) = all(is_trivial, group.gens)

## largest moved point
largest_moved_point(group::PermGp) = max(largest_moved_point.(group.gens)...)

end # module
