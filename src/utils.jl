# For application of linear binning to a univariate data set.
function linbin(X::Vector{Float64}, gpoints::AbstractRange, truncate::Bool = true)
    n = length(X)
    M = length(gpoints)

    a = gpoints[1]
    b = gpoints[M]

    trun = truncate ? 1 : 0

    res = zeros(Float64, M)

    delta = (b-a)/(M-1)

    for i in 1:n
        lxi = ((X[i]-a)/delta) + 1

        li = floor(lxi) |> Int64

        rem = lxi - li
        if li >= 1 & li < M
            res[li] += (1-rem)
            res[li+1] += rem
        end

        if li < 1 & trun == 0
            res[1] += 1
        end

        if li >= M & trun == 0
            res[m] += 1
        end
    end

    return res
end


# For application of linear binning to a regression data set.
function rlbin(X::Vector{Float64}, Y::Vector{Float64}, gpoints::AbstractRange, truncate::Bool = true)
    n = length(X)
    M = length(gpoints)
    a = gpoints[1]
    b = gpoints[M]

    trun = truncate ? 1 : 0

    xcnts = zeros(Float64, M)
    ycnts = zeros(Float64, M)

    delta = (b-a)/(M-1)

    for i in 1:n
        lxi = ((X[i]-a)/delta) + 1

        li = floor(lxi) |> Int64
        rem = lxi - li

        # Correction for right endpoint (not included if li == M)
        if X[i] == b
            li = M - 1
            rem = 1
        end

        if li >= 1 & li < M
            xcnts[li] += (1-rem)
            xcnts[li+1] += rem
            ycnts[li] += (1-rem)*Y[i]
            ycnts[li+1] += rem*Y[i]
        end

        if li < 1 & trun == 0
            xcnts[1] += 1
            ycnts[1] += Y[i]
        end

        if li >= M & trun == 0
            xcnts[M] += 1
            ycnts[M] += Y[i]
        end
    end

    return (xcnts, ycnts)
end


# Chooses the number of blocks for the preliminary
# step of a plug-in rule using Mallows' C_p.
function cpblock(X::Vector{Float64}, Y::Vector{Float64}, Nmax::Int, q::Int)
    n = length(X)

    # Sort the (X, Y) data with respect to the X's.
    sp = sortperm(X)
    X = X[sp] #don't sort in place so original x and y are not muted
    Y = Y[sp]

    # Set up arrays for FORTRAN subroutine "cp"
    qq = q + 1
    RSS = zeros(Float64, Nmax)
    Xj = zeros(n)
    Yj = zeros(n)
    coef = zeros(qq)
    Xmat = zeros(n, qq)
    Cpvals = zeros(Nmax)
    wk = zeros(n)
    qraux = zeros(qq)

    # For each number of partitions
    for Nval in 1:Nmax
        idiv = floor(Int, n/Nval) # Added the floor because idiv is defined in integer in Fortran

        # For each member of partitions
        for j in 1:Nval
            ilow = (j-1)*idiv + 1
            iupp = j*idiv

            j == Nval ? iupp = n : nothing
            
            nj = iupp - ilow + 1

            for k in 1:nj
                Xj[k] = X[ilow+k-1]
                Yj[k] = Y[ilow+k-1]
            end

            # Obtain the q'th degree fit over current member of partition
            for i in 1:nj
                Xmat[i, 1] = 1.0
                for k in 2:qq
                    Xmat[i, k] = Xj[i]^(k-1)
                end
            end
            #TODO: I have no idea if this is correct, need to check later
            coef = copy(Yj)
            LinearAlgebra.LAPACK.gels!('N', Xmat, coef)

            RSSj = 0.0

            for i in 1:nj
                fiti = coef[1]
                for k in 2:qq
                    fiti += coef[Int(k)] * Xj[Int(i)]^(k-1)
                end
                RSSj += (Yj[i] - fiti)^2
            end
            RSS[Nval] += RSSj
        end
    end

    # Now compute array of Mallow's C_p values.
    for i in 1:Nmax
        Cpvals[i] = ((n - qq*Nmax)*RSS[i]/RSS[Nmax]) + 2 * qq * i -n
    end

    # ccall((:cp_, libkernsmooth), Ptr{Void},
    #     (Ptr{Float64}, Ptr{Float64}, Ptr{Int}, Ptr{Int}, Ptr{Int}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
    #      Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
    #     X, Y, &n, &qq, &Nmax, RSS, Xj, Yj, coef, Xmat, wk, qraux, Cpvals)

    return sortperm(Cpvals)[1]
end


# For obtaining preliminary estimates of quantities required for the "direct plug-in"
# regression bandwidth selector based on blocked q'th degree polynomial fits.
function blkest(x::Vector{Float64}, y::Vector{Float64}, Nval::Int, q::Int)
    n = length(x)

    ## Sort the (x, y) data with respect to the x's.
    sp = sortperm(x)
    x = x[sp] #don't sort in place so original x and y are not muted
    y = y[sp]

    ## Set up arrays for FORTRAN program "blkest"
    qq = q + 1
    xj = zeros(n)
    yj = zeros(n)
    coef = zeros(qq)
    Xmat = zeros(n, qq)
    wk = zeros(n)
    qraux = zeros(qq)
    sigsqe = zeros(1)
    th22e = zeros(1)
    th24e = zeros(1)

    RSS = 0.0
    th22e = 0.0
    th24e = 0.0
    idiv = floor(n/Nval) |> Int # Added the floor because idiv is defined in integer in Fortran

    for j in 1:Nval
        ilow = (j-1)*idiv + 1
        iupp = j*idiv

        j == Nval ? iupp : nothing

        nj = iupp - ilow + 1

        for k in 1:nj
            xj[k] = x[ilow+k-1]
            yj[k] = y[ilow+k-1]
        end

        # Obtain a q'th degree fit over current member of partition
        # Set up "X" matrix
        for i in 1:nj
            Xmat[i,1] = 1.0
            for k in 2:qq
                Xmat[i, k] = xj[i]^(k-1)
            end
        end
        #TODO: I have no idea if this is correct, need to check later
        coef = copy(yj)
        LinearAlgebra.LAPACK.gels!('N', Xmat, coef)

        for i in 1:nj
            fiti = coef[1]
            ddm = 2*coef[3]
            ddddm = 24*coef[5]

            for k in 2:qq
                fiti += coef[k] * xj[i]^(k-1)
                if k <= (q-1)
                    ddm += k * (k+1) * coef[k+2] + xj[i]^(k-1)
                    if k <= (q-3)
                        ddddm += k * [k+1] * [k+2] * [k+3] * coef[k+4] * xj[i]^(k-1)
                    end
                end
            end
            th22e += ddm^2
            th24e += ddm * ddddm
            RSS += (yj[i] - fiti)^2
        end
    end

    sigsqe = RSS/(n-qq*Nval)
    th22e /= n
    th24e /= n

    # ccall((:blkest_, libkernsmooth), Ptr{Void},
    #       (Ptr{Float64}, Ptr{Float64}, Ptr{Int}, Ptr{Int}, Ptr{Int}, Ptr{Int}, Ptr{Float64},
    #        Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
    #        Ptr{Float64}, Ptr{Float64}),
    #       x, y, &n, &q, &qq, &Nval, xj, yj, coef, Xmat, wk, qraux, sigsqe, th22e, th24e)


    return (sigsqe[1], th22e[1], th24e[1])
end


function discretise_the_bandwidths(bandwidth, M, Q, tau, delta)
    if length(bandwidth) == M
        hlow = minimum(bandwidth)
        hupp = maximum(bandwidth)
        hdisc = [exp(h) for h in range(log(hlow), log(hupp), length=Q)]

        # Determine value of L for each member of "hdisc"
        Lvec = [floor(tau*h/delta) for h in hdisc]

        # Determine index of closest entry of "hdisc" to each member of "bandwidth"
        indic = if Q > 1
                    gap = (log(hdisc[Q])-log(hdisc[1]))/(Q-1)
                    if gap == 0
                        ones(Int, M)
                    else
                        lhlow = log(hlow)
                        [round((log(b) - lhlow)/gap + 1.0) for b in bandwidth]
                    end
                else
                    ones(Int, M)
                end
    elseif length(bandwidth) == 1
        indic = ones(Int, M)
        Q = 1
        Lvec = fill(floor(Int, tau*bandwidth/delta), Q)
        hdisc = fill(bandwidth, Q)
    else
        error("'bandwidth' must be a vector of length one or of length 'gridsize'")
    end

    return (indic, Q, Lvec, hdisc)
end


function sdiag(x::Vector{Float64}, bandwidth;
        drv::Int=0,
        degree::Int=1,
        kernel::Symbol = :normal,
        gridsize::Int = 401,
        bwdisc::Int = 25,
        range_x::Vector{Float64}=Float64[],
        binned::Bool = false,
        truncate::Bool = true
        )

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

    ## Bin the data if not already binned
    if !binned
        gpoints = range(a, b, length=M)
        xcounts = linbin(x, gpoints, truncate)
    else
        xcounts = x
        M = length(xcounts)
        gpoints = range(a, b, length=M)
    end

    ## Set the bin width
    delta = (b-a)/(M-1.0)

    (indic, Q, Lvec, hdisc) = discretise_the_bandwidths(bandwidth, M, Q, tau, delta)

    dimfkap = 2 * sum(Lvec) + Q
    fkap = zeros(dimfkap)
    midpts = zeros(Int, Q)
    ss = zeros(M, ppp)
    Smat = zeros(pp, pp)
    work = zeros(pp)
    det = zeros(2)
    ipvt = zeros(Int, pp)
    Sdg = zeros(M)

    # Obtain kernel weights
    mid = (Lvec[1] + 1) |> Int

    for i in 1:(Q-1)
        midpts[i] = mid
        fkap[mid] = 1.0
        for j in 1:Lvec[i]
            fkap[mid+j] = exp(-(delta*j/hdisc[i])) # The fortran exponentiate by 2/2(?)
            fkap[mid-j] = fkap[mid+j]
        end
        mid += Lvec[i] + Lvec[i+1] + 1
    end

    midpts[Q] = mid
    fkap[mid] = 1.0

    for j in 1:Lvec[Q]
        fkap[mid+j] = exp(-(delta*j/hdisc[Q])) # The fortran exponentiate by 2/2(?)
        fkap[mid-j] = fkap[mid+j]
    end

    # Combine kernel weights and grid counts
    for k in 1:M
        if xcounts[k] != 0
            for i in 1:Q
                for j in max(1, k-Lvec[i]):min(M, k+Lvec[i])
                    if indic[j] == i
                        fac = 1.0
                        ss[j, 1] += xcounts[k] * fkap[k-j+midpts[i]]
                        for ii in 2:ppp
                            fac *= delta * (k-j)
                            ss[j, ii] += xcounts[k]*fkap[k-j+midpts[i]]*fac
                        end
                    end
                end
            end
        end
    end

    for k in 1:M
        for i in 1:pp
            for j in 1:pp
                indss = i + j - 1
                Smat[i, j] = ss[k, indss]
            end
        end

        LinearAlgebra.LAPACK.getrf!(Smat, ipvt)
        LinearAlgebra.LAPACK.getri!(Smat, ipvt)

        Sdg[k] = Smat[1, 1]
    end

    # ccall((:sdiag_, libkernsmooth), Ptr{Void},
    #       (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int}, Ptr{Int}, Ptr{Int}, Ptr{Int},
    #        Ptr{Int}, Ptr{Float64}, Ptr{Int}, Ptr{Int}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
    #        Ptr{Float64}, Ptr{Int}, Ptr{Float64}),
    #       xcounts, &delta, hdisc, Lvec, indic, midpts, &M, &Q, fkap, &pp, &ppp, ss, Smat, work,
    #       det, ipvt, Sdg)

    return (gpoints, Sdg)
end


## For computing the binned diagonal entries of SS^T where S is a smoother matrix for
# local polynomial kernel regression.
function sstdiag(x::Vector{Float64}, bandwidth;
        drv::Int=0,
        degree::Int=1,
        kernel::Symbol = :normal,
        gridsize::Int = 401,
        bwdisc::Int = 25,
        range_x::Vector{Float64}=Float64[],
        binned::Bool = false,
        truncate::Bool = true
    )
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

    ## Bin the data if not already binned
    if !binned
        gpoints = range(a, b, length=M)
        xcounts = linbin(x, gpoints, truncate)
    else
        xcounts = x
        M = length(xcounts)
        gpoints = range(a, b, length=M)
    end

    ## Set the bin width

    delta = (b-a)/(M-1.0)

    ## Discretise the bandwidths
    (indic, Q, Lvec, hdisc) = discretise_the_bandwidths(bandwidth, M, Q, tau, delta)

    dimfkap = 2 * sum(Lvec) + Q
    fkap = zeros(dimfkap)
    midpts = zeros(Int, Q)
    ss = zeros(M, ppp)
    uu = zeros(M, ppp)
    Smat = zeros(pp, pp)
    Umat = zeros(pp, pp)
    work = zeros(pp)
    det = zeros(2)
    ipvt = zeros(Int, pp)
    SSTd = zeros(Float64, M)

    # Obtain kernel weights
    mid = Lvec[1] + 1

    for i in 1:(Q-1)
        midpts[i] = mid
        fkap[mid] = 1.0
        for j in 1:Lvec[i]
            fkap[mid+j] = exp(-(delta*j/hdisc[i])) # The fortran exponentiate by 2/2(?)
            fkap[mid-j] = fkap[mid+j]
        end
        mid += Lvec[i] + Lvec[i+1] + 1
    end

    midpts[Q] = mid
    fkap[mid] = 1.0

    for j in 1:Lvec[Q]
        fkap[mid+j] = exp(-(delta*j/hdisc[Q])) # The fortran exponentiate by 2/2(?)
        fkap[mid-j] = fkap[mid+j]
    end

    # Combine kernel weights and grid counts
    for k in 1:M
        if xcounts[k] != 0
            for i in 1:Q
                for j in max(1, k-Lvec[i]):min(M, k+Lvec[i])
                    if indic[j] == i
                        fac = 1.0
                        ss[j, 1] += xcounts[k] * fkap[k-j+midpts[i]]
                        uu[j, 1] += xcounts[k] * fkap[k-j+midpts[i]]^2
                        for ii in 2:ppp
                            fac *= delta*(k-j)
                            ss[j, ii] += xcounts[k] * fkap[k-j+midpts[i]] * fac
                            uu[j, ii] += xcounts[k] * fkap[k-j+midpts[i]]^2 * fac
                        end
                    end
                end
            end
        end
    end

    for k in 1:M
        for i in 1:pp
            for j in 1:pp
                indss = i + j - 1
                Smat[i, j] = ss[k, indss]
                Umat[i, j] = uu[k, indss]
            end
        end

        LinearAlgebra.LAPACK.getrf!(Smat, ipvt)
        LinearAlgebra.LAPACK.getri!(Smat, ipvt)

        for i in 1:pp
            for j in 1:pp
                SSTd[k] += Smat[1, i] * Umat[i, j] * Smat[j, 1]
            end
        end
    end


    # ccall((:sstdg_, libkernsmooth), Ptr{Void},
    #       (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int}, Ptr{Int}, Ptr{Int}, Ptr{Int},
    #        Ptr{Int}, Ptr{Float64}, Ptr{Int}, Ptr{Int}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
    #        Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int}, Ptr{Float64}),
    #       xcounts, &delta, hdisc, Lvec, indic, midpts, &M, &Q, fkap, &pp, &ppp, ss, uu, Smat,
    #       Umat, work, det, ipvt, SSTd)

    return (gpoints, SSTd)
end
