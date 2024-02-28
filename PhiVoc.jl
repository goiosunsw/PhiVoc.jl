
module PhiVoc
using FFTW, DSP, StatsBase

export findPeaks, highestNPeaks, pvocSignal, pvocSingleFrame

function findPeaks(x)
	pkmask = (x[1:end-2] .< x[2:end-1]) .& (x[3:end] .< x[2:end-1])
	pkmask = [false; pkmask; false]
end

function highestNPeaks(x; n=5)
	pmask = findPeaks(x)
	pkIdx = findall(pmask)
	n = min(n, length(pkIdx))
	sortIdx = sortperm(x[pmask])
	pkIdx[sortIdx[end:-1:end-n+1]]
end

function framePairsToPeakFreqs(framePair, lag::Integer, npk::Integer)
	nfft = (size(framePair,1)-1)*2
	# calculate average power spectrum for peak estimation
	frav = sum((abs.(framePair)).^2, dims=2)
	# select highest peaks
	pki = highestNPeaks(frav, n=npk)
	# phase differences between frames
	dphi = diff(angle.(framePair[pki,:]),dims=2)
	# peak bin centre frequencies
	cf = (pki.-1)./(nfft-1)
	# expected wrap factor for peak bins
	nwrap = round.(cf.*lag)
	# possible frequencies (different wrap factors possible in adjacent frames)
	fopt = (dphi/2pi .+ nwrap .+ (-1:1)' )./lag
	freqs = zeros(size(fopt,1))
	for ii in axes(fopt,1)
		ff = fopt[ii, argmin(abs.(cf[ii].-fopt[ii,:]))] 
		# FIXME: out of bounds are returned...
		if (ff > 0) && (ff < 0.5)
			freqs[ii] = ff
		end
	end
	ampls = frav[pki]
	(freqs, ampls, pki)
end
	

function pvocSingleFrame(x; nfft::Integer=1024, lag::Integer=256, window::Union{Function,Nothing}=nothing, npk::Integer=10)
	# Nyquist frame (non-redundant spectrum)
	inyq = nfft÷2+1
	# pre-allocate frames
	frames = zeros(ComplexF64,inyq,2)
	if !isnothing(window)
		win = window(nfft)
	end
	
	# calculate 2 frames of spectrum
	for (ii, ist) in enumerate([0,lag])
		xfr = x[1+ist:nfft+ist]
		if !isnothing(window)
			xfr .*= win
		end
		Xfr = fft(xfr)
		frames[:,ii] .= Xfr[1:inyq]
	end
	framePairsToPeakFreqs(frames, lag, npk)
end

function stft(x;  nfft::Integer=1024, lag::Integer=256, window::Union{Function,Nothing}=nothing)
	nx = length(x)
	if !isnothing(window)
		win = window(nfft)
	end
	
	fridx = 1:lag:(nx-nfft)
	nfr = length(fridx)
	nnyq = nfft ÷ 2 + 1
	sx = zeros(Complex{eltype(x)}, nnyq, nfr)
	for (ifr, ist) in enumerate(fridx)
		xfr = x[1+ist:nfft+ist]
		if !isnothing(window)
			xfr .*= win
		end
		Xfr = fft(xfr)
		sx[:,ifr] .= Xfr[1:nnyq]
	end
	if isnothing(window)
		wsum = sqrt(nfft)
	else
		wsum = sqrt(sum(win.^2))
	end
	wsum *= sqrt(nfft/2)
	sx./wsum
end

function pvocSignal(x; nfft::Integer=1024, lag::Integer=256, window::Union{Function,Nothing}=nothing, npk::Integer=10)
	sx = stft(x, nfft=nfft, lag=lag, window=window)
	nfr = size(sx,2)
	freqs = zeros(nfr-1, npk)
	pows = zeros(nfr-1, npk)
	nnyq = size(sx,1)
	for ii in 1:nfr-1
		framePair = sx[:,ii:ii+1]
		f, p, ib = framePairsToPeakFreqs(framePair, lag, npk)

		nf = length(f)
		freqs[ii, 1:nf] = f
		#pows[ii, 1:nf] = p./2
		# sum neighbouring bins to calculate partial power
		for (fidx, binidx) in enumerate(ib)
			minbin = max(binidx-3,1)
			maxbin = min(binidx+3,nnyq)
			pows[ii,fidx] = sum(abs.(framePair[minbin:maxbin,:]).^2)/2
		end
	end
	(freqs, pows)
end

function sinusoidTracks(freqs, pows; mindist=0.10)
	tracks = zeros(Int, size(freqs))

	# initialise track numbers for first time index
	tridx = 1
	for pidx in axes(freqs,2)
		if freqs[1,pidx] > 0
			tracks[1,tridx] = tridx
			tridx += 1
		end
	end

	# Run through the remaining partials and match to previous frame
	npart = lastindex(freqs,2)
	for tidx in 2:lastindex(freqs, 1)
		available = collect(1:npart)
		for pidx in 1:npart 
			dists = abs.(12 .* log2.(freqs[tidx-1,available] ./ freqs[tidx, pidx]))
			imin = argmin(dists)
			if (dists[imin] <= mindist) && (tracks[tidx-1,available[imin]] > 0)
				tracks[tidx,pidx] = tracks[tidx-1, available[imin]]
				filter!(x->(x!=imin), available)
			else
				if freqs[tidx-1, pidx] > 0.0
					tracks[tidx, pidx] = tridx
					tridx += 1
				end
			end
		end
	end

	tracks
end

struct SinusoidTrack
	freqs::Vector{Float64}
	powers::Vector{Float64}
	startIndex::Integer
end

function trackNumbersToTracks(freqs, pows, tracknbrs; minlength=10)
	ctr = countmap(tracknbrs)
	tracks = []
	for pno in keys(ctr)
		if (ctr[pno] > minlength) && (pno > 0)
			pmask = (tracknbrs.==pno)
			push!(tracks, SinusoidTrack(freqs[pmask],
			                            pows[pmask],
										findfirst(any(pmask, dims=2))[1]))
		end
	end
	tracks
end

end # module