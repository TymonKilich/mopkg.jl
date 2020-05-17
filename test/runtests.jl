using mopkg, Test
import InteractiveUtils: subtypes

svltf = Dict(
    (x -> x^2 - 2) => (0, -2),
    (x -> x^4 + 3x^3 + x^2 + sin(x)) => (-1.9712755825569674, -1.9710820798413156),
    (x -> x^2/3 + 2x - sin(x)) => (-3.990825624164318, -3.423528818)
)

@testset "Single variable optimizers" begin
    @testset "General test for SVOptMethods" begin
        for (fun, min) in svltf
            tval = [min[2] - 3, min[2] - 1]
            for stval in tval
                @testset "Epsilon tests" begin
                    for optim in subtypes(SVOptMethod)
                        for tolerance in [1e-2, 1e-3, 1e-4]
                            @testset "sigma tests" begin
                                for sigma in [1e-3, 1e-4, 1e-5]
                                    @test isapprox(cubic_optimize(fun, stval; eps=tolerance, sig=sigma, method=optim())[1], min[1], atol=tolerance)
                                end
                            end
                        end
                    end
                end
            end
        end
    end

end
