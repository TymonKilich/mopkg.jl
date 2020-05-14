using mopkg
import mopkg: fdc, sfdc, find_min_interval
using Test
import InteractiveUtils: subtypes
"""
Single variable global test functions – functions with single global minimum
key is function (lambda), value (y, x) in minimum
"""
svltf = Dict(
    (x -> x^2) => [0.6147540983606556 , -0.6147540983606556, -1, 1],
    (x -> x^4 + sin(x)) => [-0.8, -1.995505617977528 ,-0.8 ,-0.4],
    (x -> x^3-x^2) => [1.0341831460674156, -0.18202247191011234, -0.2, 0.2])

@testset "Single variable optimizers" begin
    @testset "Fibonacci Search tests " begin
        for (fun, min) in svltf
            @testset "Epsilon tests" begin
                for optim in subtypes(SVOptMethod)
                    for tolerance in [1e-3, 1e-4, 1e-5, 1e-6]
                        @test isapprox(line_optimize(fun, min[3], min[4]; eps=tolerance, method=optim())[1], min[1], atol=1e-2)
                    end
                end
            end
        end
    end
end
