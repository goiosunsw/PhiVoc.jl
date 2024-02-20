using Test, DSP, Revise
includet("PhiVoc.jl")
using .PhiVoc

x = ones(100)
xpeaks = [1,10,20,30]
x[xpeaks] .= [2,3,4,5]

# test peak finder
@test findall(findPeaks(x)) == [10,20,30]

# test n peak finder
@test highestNPeaks(x; n=2) == xpeaks[[end,end-1]]

# test n peak finder when n>actual nbr of peaks
@test highestNPeaks(x; n=8) == [30, 20, 10]

t = (0.0:1.0:2048.0)
y = sin.(2pi.*t./20)
@test isapprox(pvocSingleFrame(y; nfft=1024, lag=256, window=hanning, npk=10)[1][1],0.05,rtol=1e-6)

# test stft with defaults
@test size(PhiVoc.stft(y)) == (513, 5)

# test stft amplitude
@test isapprox(sum(abs.(PhiVoc.stft(y)[:,2]).^2), .5, rtol=0.1)