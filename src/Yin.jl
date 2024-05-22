module Yin

# includes
using FFTW, DSP

function AMDF(x, maxLag)
    diffs = zeros(maxLag)
    for i in (1:maxLag)
        diff = sum((x[1+i:(end-maxLag+i)] .- x[1:end-maxLag]).^2)
        diffs[i] = diff 
    end
    diffs
end

function NormCumAMDF(x, maxLag)
    lags = (1:maxLag)
    xd = AMDF(x, maxLag)
    cxd = cumsum(xd)
    xd .* (lags./cxd)
end 

function findMinima(x, threshold)
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

end # module Yin
