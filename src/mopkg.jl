module mopkg

include("svopt.jl")

function line_optimize(f, x0; eps=1e-3, maxit=1e5, method::SVOptMethod=SVHillClimb())
    optimizer = method
    optimizer(f, x0; ϵ=eps, maxiter=maxit)
end

function fibonacci_optimize(f, a,b; eps=1e-3, method::SVOptMethod=SVFiboMethod())
    optimizer = method
    optimizer(f,a,b,eps)
end

export line_optimize,fibonacci_optimize, SVOptMethod

end # module
