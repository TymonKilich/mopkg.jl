module mopkg

include("svopt.jl")

function line_optimize_1(f, a, b; eps=1e-3, method::SVOptMethod=SVFiboSeach())
    optimizer = method
    optimizer(f, a, b; eps=eps)
end

function line_optimize_2(f, x0; eps=1e-3, maxit=1e5, method::SVOptMethod=SVHillClimb())
    optimizer = method
    optimizer(f, x0; ϵ=eps, maxiter=maxit)
end

export line_optimize_1, line_optimize_2, SVOptMethod

end # module
