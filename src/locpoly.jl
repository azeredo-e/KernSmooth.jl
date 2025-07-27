## For computing a binned local polynomial
## regression estimator of a univariate regression
## function or its derivative.
## The data are discretised on an equally
## spaced grid. The bandwidths are discretised on a
## logarithmically spaced grid.
using LinearAlgebra
using SpecialFunctions


"""
    locpoly(x::Vector{Float64}, y::Vector{Float64}, bandwidth::Union{Float64, Vector{Float64}};
        drv::Int = 0,
        degree::Int=drv+1,
        kernel::Symbol = :normal,
        gridsize::Int = 401,
        bwdisc::Int = 25,
        range_x::Vector{Float64}=Float64[],
        binned::Bool = false,
        truncate::Bool = true
    )

Estimates a regression function, or its derivatives, using local polynomials.

Estimates the local polynomial regression of E[Y|X] for a series. Bandwidth is an obligatory parameter, the method uses
a fast binned implementation over an equally-spaced grid.  
If the passed bandwidth is a `Float64` the "width" of the bins for the regression is the same throughout all the data.  
If passed a `Vector{Float64}` it must first be of the same length as the optional keyword parameter `gridsize`.
The data is then equally split in "`gridsize`" parts and for each a corresponding value for the bandwidth is attributed.

# Args
- x::Vector{Float64} : X values of the series;
- y::Vector{Float64} : Y values of the series;
- bandwidth::Union{Float64, Vector{Float64}} : bandwidth size, can be a scalar for equal smoothing or Vector
for variable regression weigths;
- drv::Int=0 : order of the estimated derivative, 0 being the polynomial itself, 1 being the slope of the polynomial,
and so on;
- kernel::Symbol=:normal : type of kernel used, defaults to gaussian, currently unused;
- gridsize::Int=401 : number of equally-spaced grid points over which the function is to be estimated;
- wdisc::Int=25, : number of logarithmically-equally-spaced bandwidths on which bandwidth is discretised, to speed
up computation. More relevant when `bandwidth` is not a scalar;
- range_x::Vector{Float64}=Float64[] : vector containing the minimum and maximum values of x at which to compute
the estimate;
- binned::Bool=false : logical flag: if `true`, then x and y are taken to be grid counts rather than raw
data;
- truncate::Bool=true : logical flag: if `true`, data with x values outside the range specified by `range_x`
are ignored.

# Return
- ::T<:AbstractRange : the range (X values) of points estimated by regression;
- ::Vector{Float64} : the regression points (Y values) estimated.
"""
function locpoly(x::Vector{Float64}, y::Vector{Float64}, bandwidth::Union{Float64, Vector{Float64}};
    drv::Int = 0,
    degree::Int=drv+1,
    kernel::Symbol = :normal,
    gridsize::Int = 401,
    bwdisc::Int = 25,
    range_x::Vector{Float64} = Float64[],
    binned::Bool = false,
    truncate::Bool = true
)

    if any(bandwidth .<  0.0)
        error("'bandwidth' must be strictly positive")
    end

    if range_x == Float64[] && !binned
        maxx = maximum(x)
        minx = minimum(x)
        if y == Float64[]
            extra = 0.05 * (maxx - minx)
            range_x = [minx - extra, maxx + extra]
        else
            range_x = [minx, maxx]
        end
    end

    ## Rename common variables
    M = gridsize
    Q = bwdisc
    a = range_x[1]
    b = range_x[2]
    pp = degree + 1
    ppp = 2*degree + 1
    tau = 4

    ## Decide whether a density estimate or regression estimate is required.

    if y == Float64[]    # obtain density estimate
        n = length(x)
        gpoints = range(a, b, length=M)
        xcounts = linbin(x, gpoints, truncate)
        ycounts = (M-1) .* xcounts ./ (n*(b-a))
        xcounts = ones(M) #rep(1, M)
    else # obtain regression estimate
        ## Bin the data if not already binned
        if !binned
            gpoints = range(a, b, length=M)
            xcounts, ycounts = rlbin(x, y, gpoints, truncate)
        else
            xcounts = x
            ycounts = y
            M = length(xcounts)
            gpoints = range(a, b, length=M)
        end
    end

    ## Set the bin width
    delta = (b-a)/(M-1)

    (indic, Q, Lvec, hdisc) = discretise_the_bandwidths(bandwidth, M, Q, tau, delta)

    if minimum(Lvec) == 0
        error("Binning grid too coarse for current (small) bandwidth: consider increasing 'gridsize'")
    end

    ## Allocate space for the kernel vector and final estimate

    dimfkap = 2 * sum(Lvec) + Q
    fkap = zeros(convert(Int, dimfkap))
    curvest = zeros(M)
    midpts = zeros(Int, Q)
    ss = zeros(M, ppp)
    tt = zeros(M, pp)
    Smat = zeros(pp, pp)
    Tvec = zeros(pp)
    ipvt = zeros(Int, pp)


    mid = (Lvec[1] + 1) |> Int64
    for i in 1:(Q-1)
        midpts[i] = mid
        fkap[mid] = 1.0
        for j in 1:Lvec[i]
            fkap[mid+j] = exp(-(delta * j/hdisc[i]^(2/3)))
            fkap[mid-j] = fkap[mid+j]
        end
    end

    # TODO: Fix the weird converts everywhere
    midpts[Q] = mid
    fkap[mid] = 1.0
    for j in 1:Lvec[Q]
        fkap[convert(Int64, mid+j)] = exp(-(delta * j/hdisc[Q])) # technically we should a 2/2(?) exponent here
        fkap[convert(Int64, mid-j)] = fkap[convert(Int64, mid+j)]
    end

    for i in 1:M
        if xcounts[i] ≠ 0
            for j in 1:Q
                for k in max(1, i-Lvec[j]):min(M, i+Lvec[j])
                    if indic[convert(Int64, k)] == j
                        fac = 1.0
                        ss[convert(Int64, k), 1] += xcounts[convert(Int64, i)] * fkap[convert(Int64, i-k+midpts[j])]
                        tt[convert(Int64, k), 1] += ycounts[convert(Int64, i)] * fkap[convert(Int64, i-k+midpts[j])]
                        for ii in 2:ppp
                            fac = fac * delta * (i-k)
                            ss[convert(Int64, k), convert(Int64, ii)] += xcounts[convert(Int64, i)]*fkap[convert(Int64, i-k+midpts[j])]*fac
                            if ii <= pp 
                                tt[convert(Int64, k), convert(Int64, ii)] += ycounts[convert(Int64, i)]*fkap[convert(Int64, i-k+midpts[convert(Int64, j)])]*fac
                            end
                        end
                    end
                end
            end
        end
    end

    for i in 1:M
        for j in 1:pp
            for k in 1:pp
                indss = j + k - 1
                Smat[j, k] = ss[i, indss]
            end
            Tvec[j] = tt[i, j]
        end
        Smat, ipvt, _ = LinearAlgebra.LAPACK.getrf!(Smat)
        LinearAlgebra.LAPACK.getrs!('N', Smat, ipvt, Tvec)

        curvest[i] = Tvec[drv+1]
    end


    curvest = gamma(drv+1) .* curvest

    return (gpoints, curvest)
end

"""
    locpoly(x::Vector{Float64}, bandwidth::Union{Float64, Vector{Float64}}; args...)

Estimates a density pdf distribution of `x` using local polynomials.

For the density estimation the data is binned and the local fitting procedure is applied to the bin counts.
Bandwidth is an obligatory parameter, the method uses a fast binned implementation over an equally-spaced grid.    
If the passed bandwidth is a `Float64` the "width" of the bins for the regression is the same throughout all the data.  
If passed a `Vector{Float64}` it must first be of the same length as the optional keyword parameter `gridsize`.
The data is then equally split in "`gridsize`" parts and for each a corresponding value for the bandwidth is attributed.

# Args

The density estimation method can take the same keyword arguments as the regression method.

- x::Vector{Float64} : X values of the series;
- bandwidth::Union{Float64, Vector{Float64}} : bandwidth size, can be a scalar for equal smoothing or Vector
for variable regression weigths;
- drv::Int=0 : order of the estimated derivative, 0 being the polynomial itself, 1 being the slope of the polynomial,
and so on;
- kernel::Symbol=:normal : type of kernel use, defaults to gaussian, currently unused;
- gridsize::Int=401 : number of equally-spaced grid points over which the function is to be estimated;
- wdisc::Int=25, : number of logarithmically-equally-spaced bandwidths on which bandwidth is discretised, to speed
up computation. More relevant when `bandwidth` is not a scalar;
- range_x::Vector{Float64}=Float64[] : vector containing the minimum and maximum values of x at which to compute
the estimate;
- binned::Bool=false : logical flag: if `true`, then x and y are taken to be grid counts rather than raw
data;
- truncate::Bool=true : logical flag: if `true`, data with x values outside the range specified by `range_x`
are ignored.

# Return
- ::T<:AbstractRange : the range (X values) of values for the density estimation;
- ::Vector{Float64} : the Y values for the estimated density curve.
"""
locpoly(x::Vector{Float64}, bandwidth::Union{Float64, Vector{Float64}}; args...) = locpoly(x, Float64[], bandwidth, args...)
