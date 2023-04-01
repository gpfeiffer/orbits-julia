##  cyclic shifts

module shifts

export cyclic_shifts, cyclic_shifts_with_edges

using orbits
using coxeter

function cyclic_shifts(W, w)
    function byCyclicShift(x, s)
        y = x^s
        coxeterLength(W, x) == coxeterLength(W, y) ?  y : x
    end
    orbit(W.gens, w, byCyclicShift)
end

##  the cyclic shift class of an involution x is just {x}

function cyclic_shifts_with_edges(W, w)
    function byCyclicShift(x, s)
        y = x^s
        coxeterLength(W, x) == coxeterLength(W, y) ?  y : x
    end
    orbit_with_edges(W.gens, w, byCyclicShift)
end

end # module
