module mopkg

include("svopt.jl")

function line_optimize(f, x0; eps=1e-3, maxit=1e5, method::SVOptMethod=SVFiboSearch())
    optimizer = method
    optimizer(f, x0; ϵ=eps, maxiter=maxit)  
end

export line_optimize, SVOptMethod

end # module
