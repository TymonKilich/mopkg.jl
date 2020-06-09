module mopkg

include("svopt.jl")

function line_optimize(f, a, b; eps=1e-3, method::SVOptMethod=SVHillClimb())
    optimizer = method
    optimizer(f, a, b)
end

export line_optimize, SVOptMethod

end # module
