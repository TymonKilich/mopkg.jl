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
    (x -> (x-3)^4+5) => (5, 3))

@testset "Single variable optimizers" begin
    @testset "General test for SVOptMethods" begin
        for (fun, min) in svltf
            tval = [min[2] - 2, min[2] + 3]
            for stval in tval
                @testset "Epsilon tests" begin
                    for optim in subtypes(SVOptMethod)
                        for tolerance in [1e-2, 1e-4, 1e-6]
                            @testset "Max steps test" begin
                                for max_steps in [50, 100, 200]
                                    @test isapprox(line_optimize(fun, stval; eps=tolerance, maxit=max_steps, method=optim())[1], min[1], atol=tolerance)
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
