module mopkg

include("svopt.jl")

function line_optimize(f, x0; eps=1e-3, maxit=1e5, method::SVOptMethod=SVHillClimb())
    optimizer = method
    optimizer(f, x0; ϵ=eps, maxiter=maxit)
end

function line_optimize2(f, a, b; eps=1e-3, maxit=1e5, method::SVOptMethod=SVSecant())
    optimizer = method
    optimizer(f, a, b; ϵ=eps, maxiter=maxit)
end
export line_optimize, line_optimize2, SVOptMethod

end # module
