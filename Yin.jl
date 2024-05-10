module Yin

# includes
using FFTW, DSP

function AMDF(x)
    nx = length(x)
    xx = x-mean(x)
    # C0 = (((xx')*xx))
    Cx = conv(xx,reverse(xx))
    Cx[nx] .- 2 Cx 
end


end # module Yin
