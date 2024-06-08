module Yin

# includes
using FFTW, DSP, PaddedViews

struct YinParams
    fs::Real
    fmin::Real
    fmax::Real
    nwind::Int
    nhop::Int
    threshold::Real
    fast::Bool
    minLag::Int
    maxLag::Int
end

function YinParams(fs::Real; fmin::Real=50, fmax::Real=4000, 
    threshold::Real=.1, nwind::Int=2048, nhop::Int=-1, fast::Bool=true)
    if nhop < 0
        nhop = nwind ÷ 2
    end
    YinParams(fs, fmin, fmax, nwind, nhop, threshold, fast, floor(fs/fmax), ceil(fs/fmin))
end

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

function yinFrame(frame::Vector{<: Number}, params::YinParams)
    camdf = NormCumAMDF(frame, params.maxLag)
    lag = findFirstMinimum(camdf, params.threshold)
    if lag <= 0
        return (-1, minimum(camdf))
    end
    (fineLag, val) = parabolicInterp(camdf, lag)
    (fineLag/params.fs, val)
end

function yin(x, params)
    startIdx = 1:params.nhop:(length(x)-params.nwind)
    nFrames = length(startIdx)
    freqs = Vector{Union{Missing,Float64}}(undef, nFrames)
    vals = zeros(nFrames)
    for (ii, stIdx) in enumerate(startIdx)
        frame = x[stIdx:stIdx+params.nwind]
        (lag, val) = yinFrame(frame, params)
        if lag > 0
            freqs[ii] = 1/lag
        end
        vals[ii] = val
    end
    (freqs, vals)

end

end # module Yin
