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
    (x -> x^4 - sin(x)) => (-0.4351967653715638, 0.5919189453125))

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
    @testset "SVDivInHalf tests" begin
        @test line_optimize(x -> x^4 + 2sin(x), -10; eps=1e-5, method=subtypes(SVOptMethod)[1]()) == (-1.0500382539343835, -0.7214080810546877)
        @test line_optimize(x -> x/(x^2+1), -10; eps=1e-5, method=subtypes(SVOptMethod)[1]()) == (-0.4999999999994179, -1.0000015258789068)
    end
end
