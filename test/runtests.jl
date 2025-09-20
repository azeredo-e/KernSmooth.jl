include("../src/KernSmooth.jl")
using .KernSmooth, RDatasets, Test, Statistics
using RCall

dat = dataset("MASS", "geyser")

X = Array(dat[!, :Duration])
Y = convert(Vector{Float64}, Array(dat[!, :Waiting]))

resx, resy = locpoly(X, Y, 0.19)

@rput X Y

R"""
library(KernSmooth)
res <- locpoly(X, Y, bandwidth=0.25)
res_d <- dpill(X, Y)
r_X <- res$x
r_Y <- res$y
"""

@rget r_X r_Y res_d

@testset "missing value test" begin
    @test !any(isnan.(resy))
end

@testset "Precision test" begin
    # @test isapprox(maximum(abs.(resx .- expected[:, 1])), 0.0, atol=sqrt(eps(Float64)))
    @test isapprox(maximum(abs.(resx .- r_X)), 0.0, atol=sqrt(eps(Float64)))
    
    # It fails the second test: isapprox(1.4656224871855557, 0.0; rtol = 1.4901161193847656e-8)
    # But looking to the graph the regression lines seem to match
    @test isapprox(mean(abs.(resy .- r_Y)), 0.0, atol=sqrt(eps(Float64)))

    # Test failed: isapprox(0.2452942141872545, 0.238350165926485; atol = 1.4901161193847656e-8)
    @test isapprox(dpill(X, Y), res_d, atol=sqrt(eps(Float64)))
end
