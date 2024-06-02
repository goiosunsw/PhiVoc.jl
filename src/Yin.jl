module Yin

# includes
using FFTW, DSP, PaddedViews

function AMDFtd(x, maxLag)
    n = length(x)
    maxLag = min(maxLag, n)
    diffs = zeros(maxLag)
    for i in (1:maxLag)
        diff = sum((x[1+i:(end-maxLag+i)] .- x[1:end-maxLag]).^2)
        diffs[i] = diff 
    end
    diffs
end

function AMDF(x, maxLag)
    n = length(x)
    maxLag = min(maxLag, n)
    xcs = [0; cumsum(x.^2)]
    nfft = nextfastfft(n+maxLag)
    xp = PaddedView(0., x, (nfft,))
    xf = rfft(xp)
    conv = irfft(xf .* conj(xf), nfft)[nfft:-1:nfft-maxLag+1]
    xcs[n:-1:(n-maxLag+1)] .+ xcs[n+1] .- xcs[1:maxLag] .- conv * 2
end

function NormCumAMDF(x, maxLag)
    lags = (1:maxLag)
    xd = AMDF(x, maxLag)
    cxd = cumsum(xd)
    xd .* (lags./cxd)
end 

function findFirstMinimum(x, threshold)
    local i
    for outer i in axes(x,1)
        if x[i] < threshold
            break
        end
    end

    while ((i < lastindex(x)) && (x[i+1] < x[i]))
        i+=1
    end 

    if i < lastindex(x)
        i
    else
        -1
    end
end

function parabolicInterp(v::Vector, idx::Integer)
    n = length(v)
    if idx == 1
        return (1, v[1])
    elseif idx == n
        return (n, v[n])
    else
        vl = v[idx-1]
        vr = v[idx+1]
        vc = v[idx]
        a = vl/2 + vr/2 - vc
        b = (vr - vl) / 2
        di = (vr - vl) / 2 / (2vc - vl - vr)
        vi = a * di^2 +b * di + vc
    end
    (idx+di, vi)
end

end # module Yin
