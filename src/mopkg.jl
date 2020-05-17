module mopkg
include("svopt.jl")
function line_optimize(f, x1, x2, x3; eps=0.5, method::SVOptMethod=SVPowell())  
    optimizer = method   
    optimizer(f, x1, x2, x3; epsilon=eps)
end  
            
export line_optimize, SVOptMethod

end # module