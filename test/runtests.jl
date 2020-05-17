using mopkg 
using Test 
import InteractiveUtils: subtypes 
""" 
Single variable global test functions – functions with single global minimum 
key is function (lambda), value (y, x) in minimum 
""" 
svltf = Dict( 
    (x -> x^2 - 1) => (-1.0, 0.0),
    (x -> x^2/3 + 2x - sin(x)) => (-3.423528818, -3.99083),
    (x -> x^4 + 3x^3 + x^2 + sin(x))=> (-4.915216432826826,-1.971275676169661)) 
  
@testset "Single variable optimizers" begin 
    @testset "General test for SVOptMethods" begin 
        for (fun, min) in svltf
            stval = min[2] - 2
            midval = min[2] + 3
            endval = min[2] + 6
             @testset "Epsilon tests" begin 
                for optim in subtypes(SVOptMethod) 
                    for tolerance in [1e-2, 1e-4, 1e-6] 
                        @test isapprox(line_optimize(fun, stval, midval, endval; eps=tolerance, method=optim())[1], min[1], atol = tolerance) 
                    end
                end 
            end
             
        end 
    end 
end
