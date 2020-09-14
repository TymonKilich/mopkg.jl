using mopkg
import mopkg: fdc, sfdc, find_min_interval, fw1, fw2
using Test
import InteractiveUtils: subtypes
"""
Single variable global test functions – functions with single global minimum
key is function (lambda), value (y, x) in minimum
"""
svltf = Dict(
    (x -> x^2 - 1) => (-1.0, 0.0),
    (x -> x^2/3 + 2x - sin(x)) => (-3.423528818, -3.99083),
    (x -> (x-1)^2*10 + 10) => (10, 1)
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

    @testset "Golden ratio helper functions" begin
        ab_intervals = [
            [5,6],
            [0,0.5],
            [-1000.1,-3.14],
            [-10, 200]
        ]
        for (a, b) in ab_intervals
                @test fw1(a, b) - a == b - fw2(a, b)
        end
    end
end
