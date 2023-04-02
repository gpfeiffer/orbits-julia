#############################################################################
##
##  Variants of Relations  (Lists Version.)
##
module variants

export variantsRelations

function variantsRelations(genrel)

    abs(s::Int) = s < 0 ? genrel.invr[-s] : s

    variants = [Set{Vector{Int}}() for _  in genrel.gens]
    for list in genrel.rels
        relator =  vcat(list[1], -reverse(list[2]))

        #  u s v = 1, s = (v u)^{-1},  s^-1 = v u
        for (i, s) in enumerate(relator)
            vu = vcat(relator[i+1:end], relator[1:i-1])
            push!(variants[abs(s)], abs.(-reverse(vu)))
            push!(variants[abs(-s)], abs.(vu))
        end
    end

    [sort(collect(list), by=length) for list in variants]
end

end # module
