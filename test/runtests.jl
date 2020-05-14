using mopkg
import mopkg: fdc, sfdc, find_min_interval
using Test
import InteractiveUtils: subtypes
"""
Single variable global test functions – functions with single global minimum
key is function (lambda), value (y, x) in minimum
"""
svltf = Dict(
    (x -> x^2 - 1) => (-1.0, 0.0),
    (x -> x^2/3 + 2x - sin(x)) => (-3.423528818, -3.99083),
    (x -> x^4 + 3x^3 + x^2 + sin(x)) => (-4.915216432826815, -1.9712755845857126)
)

@testset "Single variable optimizers" begin
    @testset "General test for SVOptMethods" begin
        for (fun, min) in svltf
            tval = [min[2] - 2, min[2] + 3]
            for stval in tval
                @testset "Epsilon tests" begin
                    for optim in subtypes(SVOptMethod)
                        for tolerance in [1e-2, 1e-4, 1e-6]
                            @test isapprox(line_optimize(fun, stval; eps=tolerance, method=optim())[1], min[1], atol=tolerance)
                        end
                    end
                end
            end
        end
    end
    @testset "Golden section method test" begin
        @testset "General test (SVGoldenSection)" begin
            f, x0, min = x -> x^2, 5.0, 0.0
            for tolerance in [1e-2, 1e-4, 1e-6]
                @test isapprox(line_optimize(f, x0; method=mopkg.SVGoldenSection())[1], min, atol=tolerance)
            end
        end
        @testset "Precision test (SVGoldenSection)" begin
            f, x0, min = x -> x^4 + 3x + 12, 3.0, 9.9557394291
            for tolerance in [1e-2, 1e-4, 1e-6, 1e-8, 1e-10]
                @test isapprox(line_optimize(f, x0; method=mopkg.SVGoldenSection())[1], min, atol=tolerance)
            end
        end
    end
    @testset "Finite differences tests" begin
        for xval in [-π, -π/2, 0, π/2, π]
            @test isapprox(2cos(2xval), mopkg.fdc(x -> sin(2x), xval), atol=1e-5)
        end
        @test isapprox(3^3, mopkg.fdc(x -> x^3, 3), atol=1e-5)
        @test isapprox(4*0.05^3, mopkg.fdc(x -> x^4, 0.05), atol=1e-5)
        for xval in [-π, -π/2, 0, π/2, π]
            @test isapprox(-4sin(2xval), mopkg.sfdc(x -> sin(2x), xval), atol=1e-5)
        end
        @test isapprox(12*0.05^2, mopkg.sfdc(x -> x^4, 0.05), atol=1e-5)
    end
    @testset "Find minimum interval" begin
        for (fun, min) in svltf
            tval = [min[2] - 2, min[2] + 3]
            for stval in tval
                @test find_min_interval(fun, stval)[1] ≤ min[2] ≤ find_min_interval(fun, stval)[2]
            end
        end
    end
end
