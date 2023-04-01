##  cyclic shifts

module shifts

export cyclic_shifts, cyclic_shifts_with_edges

using orbits, coxeter

function cyclic_shifts(W::CoxeterGp, w)
    function byCyclicShift(x, s)
        y = x^s
        coxeterLength(W, x) == coxeterLength(W, y) ?  y : x
    end
    orbit(W.gens, w, byCyclicShift)
end

function cyclic_shifts_with_edges(W::CoxeterGp, w)
    function byCyclicShift(x, s)
        y = x^s
        coxeterLength(W, x) == coxeterLength(W, y) ?  y : x
    end
    orbit_with_edges(W.gens, w, byCyclicShift)
end

end # module
