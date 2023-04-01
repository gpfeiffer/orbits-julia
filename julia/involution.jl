module involution

onInvolutions(x, s) =  x * s == s * x ? onRight(x, s) : onPoints(x, s)

involutions(W) = orbit(W.gens, W.one, onInvolutions)

function onInvolutionClasses(x, a)
    y = x.elts[1]
    y^a <> y ? x : onRight(y, a)^x.group
end

involutionClasses(W) = orbit(reflections(W), W.one^W, onInvolutionClasses)

end # module
