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
    (x -> x^2/3 + 2x - sin(x)) => (-3.423528818, -3.99083))
"""

Metoda siecznych Działająca


"""

svltf2 = Dict(
	(x-> x^2 - 1)=> (1.0000000010382313,2.076462557454306e-9),
	(x->sin(x)+2x-3)=>(1.0630731347927889,4.176037293746049e-11),
	(x->14x-30)=>(2.142857142857143, 0.0)
	)


@testset "Metoda Siecznych" begin
	for (fun, min) in svltf2
		stval=(min[2]+1)
		stval1=(min[2]+5)
		@testset "Test Epsilonu" begin
			for optim in subtypes(SVOptMethod)
				for tolerance in [1e-2, 1e-4, 1e-6]
					@test isapprox(secants_optimize(fun, stval,stval1; EPS=tolerance, method=mopkg.SVSecant())[2], min[1] , atol= tolerance)
				end
			end
		end
	end
end

@testset "Single variable optimizers" begin
    @testset "General test for SVOptMethods" begin
        for (fun, min) in svltf
            tval = [min[2] - 2, min[2] + 3]
            for stval in tval
                @testset "Epsilon tests" begin
                    for optim in subtypes(SVOptMethod)
                        for tolerance in [1e-2, 1e-4, 1e-6]
                            @test isapprox(line_optimize(fun, stval; eps=tolerance, method=mopkg.SVHillClimb())[1], min[1], atol=tolerance)
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
end
