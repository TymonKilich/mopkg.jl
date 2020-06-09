using mopkg
import mopkg: fdc, sfdc, find_min_interval
using Test
import InteractiveUtils: subtypes
"""
Single variable global test functions – functions with single global minimum
key is function (lambda), value (x, y) in minimum
"""
svltf = Dict(
   	(x -> x^2 - 1) => (0.0, -1.0),
	(x -> x^2 - 2x) => (1.0, -1.0),
	(x -> 2x^2 - 4x - 6) => (1.0, -8.0),
	(x -> x^2 + 2x - sin(x)) => (-0.58, -0.27)
)

@testset "Single variable optimizers" begin
    @testset "General test for SVOptMethods" begin
        for (fun, min) in svltf
	    a = min[1] - 10
	    b = min[1] + 10 
                @testset "Epsilon tests" begin
                    for optim in subtypes(SVOptMethod)
                        for tolerance in [0.5, 0.2, 0.1]
                            @test isapprox(line_optimize(fun, a, b; eps=tolerance, method=optim())[1], min[1], atol=tolerance)
                        end
                    end
                end
        end
    end

"""Test check the lengh of range"""

    @testset "Range test" begin
	for (fun, min) in svltf
	    a = min[1] - 10
	    b = min[1] + 10 
	    for optim in subtypes(SVOptMethod)
                        for tolerance in [0.5, 0.2, 0.1]

                            @test abs(line_optimize(fun, a, b; eps=tolerance, method=optim())[4])-abs(line_optimize(fun, a, b; eps=tolerance, method=optim())[3]) < tolerance 

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
end
