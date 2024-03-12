using .PhiVoc
using GLMakie

function plot3SinusoidalTracks(tracks::Vector{SinusoidTrack})
    fig = Figure()
    ax = Axis3(fig[1,1])
    for t in tracks
        tt = t.startIndex - 1 .+ (1:length(t.freqs))
        lines!(ax, tt, t.freqs, 10 .* log10.(t.powers))
    end
    fig
end

function plotSinusoidalTracks(tracks::Vector{SinusoidTrack})
    fig = Figure()
    ax = Axis(fig[1,1])
    for t in tracks
        tt = t.startIndex - 1 .+ (1:length(t.freqs))
        GLMakie.scatter!(ax, tt, t.freqs, color=10 .* log10.(t.powers))
    end
    fig
end