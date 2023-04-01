#############################################################################
##
module coxeter

using permutation
using orbits

import permutation: Perm

export coxeterGraph, cartanMat, CoxeterGp, coxeterConjugacyClasses
export coxeterLength

using LinearAlgebra

function coxeterGraph(series, rank)
    edges = [(j-1,j) for j in 2:rank]
    series >= "D" && (edges[1] = (1,3))
    series >= "E" && (edges[2] = (2,4))
    return edges
end


function cartanMat(series, rank)
    cartan = Matrix(2I, rank, rank)
    for (i,j) in coxeterGraph(series, rank)
        cartan[i,j] = -1
        cartan[j,i] = -1
    end
    series == "B" && (cartan[1,2] = -2)
    series == "C" && (cartan[2,1] = -2)
    series == "F" && (cartan[3,4] = -2)  # sic!
    return cartan
end

absRoot(root) = sign(sum(root)) * root
onRoots(x, a) = absRoot(onRight(x, a))

function orbits_with_words_and_edges(aaa, xxx, under)
    words = Array{Int}[[i] for i in eachindex(xxx)]
    list = copy(xxx)
    edges = []
    i = 0
    while i < length(list)
        i += 1
        for (k, a) in enumerate(aaa)
            z = under(list[i], a)
            l = findfirst(==(z), list)
            if l == nothing
                push!(list, z)
                push!(words, onWords(words[i], k))
                l = length(list)
            end
            push!(edges, (i, l))
        end
    end
    return (list = list, edges = edges, words = words)
end

Perm(a, xxx, under) = Perm(indexin([under(x, a) for x in xxx], xxx))

struct CoxeterGp
    gens::Vector{Perm}
    one::Perm
    data::Dict{Symbol, Any}
#    CoxeterGp(gens, one) = new(gens, one, Dict())
end

function CoxeterGp(C)
    one = C^0
    S = axes(C,1)
    mats = []
    for s in S
        mat = C^0
        mat[:,s] = one[s,:] - C[s,:]
        push!(mats, mat)
    end
    eee = [hcat(a...) for a in eachrow(one)]
    roots = orbits_with_words_and_edges(mats, eee, onRoots)
    data = Dict(:mats => mats, :roots => roots, :rank => length(S))
    data[:N] = length(roots.list)
    data[:phi] = vcat(roots.list, -roots.list);
    data[:perms] = [Perm(m, data[:phi], onRight) for m in mats]
    CoxeterGp(data[:perms], data[:perms][1]^0, data)
end

import permgroup: PermGp
PermGp(group::CoxeterGp) = PermGp(group.gens, group.one)

using permgroup
size(group::CoxeterGp) = sizeOfGroup(PermGp(group))

data(group::CoxeterGp) = group.data

function reflections(W::CoxeterGp)
    reflection(w) = W.gens[w[1]]^prod(W.gens[w[2:end]]; init=W.one)
    reflection.(data(W)[:roots].words)
end

function coxeterLength(W, w)
    N = data(W)[:N]
    count(i^w > N for i in 1:N)
end

permCoxeterWord(W, word) = prod(W.gens[word]; init = W.one)

isLeftDescent(W, w, s) = s^w > data(W)[:N]

function firstDescent(W, w)
    n, N = getindex.(Ref(data(W)), [:rank, :N])
    s = 0
    while s < n
        s += 1
        s^w > N && return s
    end
end

function coxeterWord(W, w)
    word = Int[]
    while !is_trivial(w)
        a = firstDescent(W, w)
        push!(word, a)
        w = W.gens[a] * w
    end
    return word
end

reducedWord(W, word) = coxeterWord(W, permCoxeterWord(W, word))

function longestElement(W, J)
    wJ = W.one
    N = data(W)[:N]
    J = collect(J)
    while true
        i = findfirst(s -> s^wJ <= N, J)
        i == nothing && return wJ
        wJ = W.gens[J[i]] * wJ
    end
end

function prefixes(W, w)
    onRightDescents(w, s) = isLeftDescent(W, w^-1, s) ? w * W.gens[s] : w
    return orbit(1:data(W)[:rank], w, onRightDescents)
end

function prefixes_with_edges(W, w)
    onRightDescents(w, s) = isLeftDescent(W, w^-1, s) ? w * W.gens[s] : w
    orbit_with_edges(1:data(W)[:rank], w, onRightDescents)
end

longestCosetElement(W, J, L) = longestElement(W, J) * longestElement(W, L)

parabolicTransversal(W, J) = prefixes(W, longestCosetElement(J, 1:data(W)[:rank]))

tackOn(x, s) = sort(union(x, s))

takeAway(x, s) = sort(setdiff(x, s))

onSortedTuples(tup, a) = sort([x^a for x in tup])

function shape(W, J)
    function onParabolics(K, s)
        return onSortedTuples(K, longestCosetElement(W, K, tackOn(K, s)))
    end
    sort(orbit(1:data(W)[:rank], J, onParabolics))
end

function shapes(W)
    onShapes(x, s) = shape(W, takeAway(x[1], s))
    S = 1:data(W)[:rank]
    orbit(S, shape(W, collect(S)), onShapes)
end

function shape_with_edges(W, J)
    function onParabolics(K, s)
        onSortedTuples(K, longestCosetElement(W, K, tackOn(K, s)))
    end
    orbit_with_edges(1:data(W)[:rank], J, onParabolics)
end

function shape_with_transversal(W, J)
    S = 1:data(W)[:rank]
    list = [J]
    reps = [W.one]
    i = 0
    while i < length(list)
        i += 1
        K = list[i]
        for s in setdiff(S, K)
            a = longestCosetElement(W, K, tackOn(K, s))
            L = onSortedTuples(K, a)
            if !(L in list)
                push!(list, L);
                push!(reps, reps[i] * a);
            end
        end
    end
    return (list = list, reps = reps)
end

function parabolicComplement(W, J)
    S = 1:data(W)[:rank]
    list = [J]
    reps = [W.one]
    gens = (ears = Set(Perm[]), eyes = Set(Perm[]))
    i = 0
    while i < length(list)
        i += 1
        K = list[i]
        for s in setdiff(S, K)
            a = longestCosetElement(W, K, tackOn(K, s))
            L = onSortedTuples(K, a)
            j = findfirst(==(L), list)
            if j == nothing
                push!(list, L)
                push!(reps, reps[i] * a)
            elseif j == i
                push!(gens.ears, reps[i] * a * reps[j]^-1)
            else
                push!(gens.eyes, reps[i] * a * reps[j]^-1)
            end
        end
    end
    return gens
end

function minConjugates(W, x)
    list = [x]
    lx = coxeterLength(W, x)
    i = 0
    while i < length(list)
        i += 1;
        y = list[i]
        for s in W.gens
            z = y^s
            lz = coxeterLength(W, z)
            if lz == lx
                z in list || push!(list, z)
            elseif lz < lx
                list = [z]
                lx = lz
                i = 0
                break
            end
        end
    end
    return list
end

import permutation: is_trivial, largest_moved_point
is_trivial(group::CoxeterGp) = all(is_trivial, group.gens)
largest_moved_point(group::CoxeterGp) = max(largest_moved_point.(group.gens)...)

function coxeterMinRep(W, w)
    v = minConjugates(W, w)[1]
    K = unique(sort(coxeterWord(W, v)))
    sh = shape_with_transversal(W, K)
    (L, j) = findmin(sh.list)
    minimum(minConjugates(W, v^sh.reps[j]))
end

function coxeterConjugacyClasses(W)
    onMinReps(x, a) = coxeterMinRep(W, x * a)
    orbit(reflections(W), W.one, onMinReps)
end

end # module
