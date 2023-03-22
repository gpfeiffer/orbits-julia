#############################################################################
##
#A  sims.jl                                                      orbits-julia
#B    by Götz Pfeiffer <goetz.pfeiffer@nuigalway.ie>
##
#C  A simple minded implementation of the Schreier-Sims algorithm
##
module simsgroup

import Base: in, size, rand
import permutation: Perm, is_trivial, largest_moved_point

export SimsGp, orbit_sims
export cube

struct SimsGp
    gens::Vector{Perm}
    one::Perm
    sims::Dict{Symbol, Any}
    SimsGp(gens, one) = new(gens, one, Dict())
end

is_trivial(group::SimsGp) = all(is_trivial, group.gens)

largest_moved_point(group::SimsGp) = max(largest_moved_point.(group.gens)...)

function in(a::Perm, group::SimsGp)
    is_trivial(group) && return is_trivial(a)
    s = sims(group)
    pos = findfirst(==(s[:list][1]^a), s[:list])
    pos == nothing && return false
    return a/s[:reps][pos] in s[:stab]
end

function closure(group::SimsGp, a::Perm)
    a in group && return group
    return SimsGp(vcat(group.gens, a), group.one)
end

function orbit_sims(aaa, x, under=^)
    list = [x]
    reps = [aaa[1]^0]
    stab = SimsGp([], aaa[1]^0)
    i = 0
    while i < length(list)
        i += 1
        for a in aaa
            z = under(list[i], a)
            l = findfirst(==(z), list)
            if l == nothing
                push!(list, z)
                push!(reps, reps[i] * a)
            else   # x^(reps[i] * a) = x^reps[l]
                stab = closure(stab, reps[i] * a / reps[l])
            end
        end
    end
    return Dict(:list => list, :reps => reps, :stab => stab)
end

function sims(group::SimsGp)
    group.sims == Dict() && merge!(group.sims,
      orbit_sims(group.gens, largest_moved_point(group)))
    return group.sims
end

function size(group::SimsGp)::BigInt
    is_trivial(group) && return 1
    length(sims(group)[:list]) * size(sims(group)[:stab])
end

function rand(group::SimsGp)
    is_trivial(group) && return group.one
    rand(sims(group)[:stab]) * rand(sims(group)[:reps])
end

## test data

##                 +--------------+
##                 |              |
##                 |  1    2    3 |
##                 |              |
##                 |  4  top    5 |
##                 |              |
##                 |  6    7    8 |
##                 |              |
##  +--------------+--------------+--------------+--------------+
##  |              |              |              |              |
##  |  9   10   11 | 17   18   19 | 25   26   27 | 33   34   35 |
##  |              |              |              |              |
##  | 12  left  13 | 20 front  21 | 28 right  29 | 36  rear  37 |
##  |              |              |              |              |
##  | 14   15   16 | 22   23   24 | 30   31   32 | 38   39   40 |
##  |              |              |              |              |
##  +--------------+--------------+--------------+--------------+
##                 |              |
##                 | 41   42   43 |
##                 |              |
##                 | 44 bottom 45 |
##                 |              |
##                 | 46   47   48 |
##                 |              |
##                 +--------------+
##
##
cube = SimsGp([Perm(48, cycles) for cycles in [
  [[ 1, 3, 8, 6],[ 2, 5, 7, 4],[ 9,33,25,17],[10,34,26,18],[11,35,27,19]],
  [[ 9,11,16,14],[10,13,15,12],[ 1,17,41,40],[ 4,20,44,37],[ 6,22,46,35]],
  [[17,19,24,22],[18,21,23,20],[ 6,25,43,16],[ 7,28,42,13],[ 8,30,41,11]],
  [[25,27,32,30],[26,29,31,28],[ 3,38,43,19],[ 5,36,45,21],[ 8,33,48,24]],
  [[33,35,40,38],[34,37,39,36],[ 3, 9,46,32],[ 2,12,47,29],[ 1,14,48,27]],
  [[41,43,48,46],[42,45,47,44],[14,22,30,38],[15,23,31,39],[16,24,32,40]],
]], Perm(48))
end # module
