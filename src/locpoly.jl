## For computing a binned local polynomial
## regression estimator of a univariate regression
## function or its derivative.
## The data are discretised on an equally
## spaced grid. The bandwidths are discretised on a
## logarithmically spaced grid.
using LinearAlgebra
using SpecialFunctions

function locpoly(x::Vector{Float64}, y::Vector{Float64}, bandwidth::Union{Float64, Vector{Float64}};
    drv::Int = 0,
    degree::Int=drv+1,
    kernel::Symbol = :normal,
    gridsize::Int = 401,
    bwdisc::Int = 25,
    range_x::Vector{Float64}=Float64[],
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

locpoly(x::Vector{Float64}, bandwidth::Union{Float64, Vector{Float64}}; args...) = locpoly(x, Float64[], bandwidth, args...)
