module mopkg

include("svopt.jl")

function line_optimize(f, x0; eps=1e-3, maxit=1e5, method::SVOptMethod=SVHillClimb())
    optimizer = method
    optimizer(f, x0; ϵ=eps, maxiter=maxit)
end

function secant_optimize(f, x0, x1; eps=1e-3, maxit=1e5, method::SVOptMethod=my_secant())
    optimizer = method
    
    optimizer(f, x0, x1; ϵ=eps, iter=maxit)
end

export line_optimize, secant_optimize, SVOptMethod

end # module
