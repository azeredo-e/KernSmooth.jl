include("../src/KernSmooth.jl")
using .KernSmooth, RDatasets, Test, CSV

expected = CSV.File("geyserfit.csv") |> Array
expected = stack(expected, dims=1)

dat = dataset("MASS", "geyser")

X = Array(dat[!, :Duration])
Y = convert(Vector{Float64}, Array(dat[!, :Waiting]))

resx, resy = locpoly(X, Y, 0.25)

@testset "missing value test" begin
    @test !any(isnan.(resy))
end

@testset "Precision test" begin
    @test isapprox(maximum(abs.(resx .- expected[:, 1])), 0.0, atol=sqrt(eps(Float64)))
    
    # It fails the second test: isapprox(1.4656224871855557, 0.0; rtol = 1.4901161193847656e-8)
    # But looking to the graohs the regression lines seem to match
    # One thing is that when I adjust the value of the bandwidth to 0.19 the first value came out as NaN, should investigate
    @test isapprox(maximum(abs.(resy .- expected[:, 2])), 0.0, rtol=sqrt(eps(Float64)))
end

#example from R function
# @test_approx_eq dpill(x, y) 0.238350165926485
