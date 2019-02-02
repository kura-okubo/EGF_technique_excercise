module EGF_TOOLS

# cross-correlation module
export sac2seisjl, cumtrapz, vel2disp_withfiltering
using FFTW, Plots
using ..SeisJul

include("correlate_kurama.jl")
include("filter.jl")
using .Correlate
using .Filter

"""
    sac2seisjl(sacdata::struct)

store station data from SAC to SeisJul format.
"""


abstract type AbstructSacdata end

mutable struct sdata <: AbstructSacdata

    delta       ::Float64
    stationla   ::Float64
    stationlo   ::Float64
    stationel   ::Float64
    eventla     ::Float64
    eventlo     ::Float64
    eventdp     ::Float64
    dist        ::Float64
    nzyear      ::Int32
    nzday       ::Int32
    nzhour      ::Int32
    nzmin       ::Int32
    nzsec       ::Int32
    nzmsec      ::Int32
    npts        ::Int32
    kstnm       ::AbstractString
    kevnm       ::AbstractString
    kcmpnm      ::AbstractString
    x           ::AbstractArray

    sdata() = new()

end


function sac2seisjl(sacdata)

    seisjl_temp             = sdata()
    seisjl_temp.delta       = sacdata.delta
    seisjl_temp.stationla   = sacdata.stla
    seisjl_temp.stationlo   = sacdata.stlo
    seisjl_temp.stationel   = sacdata.stel
    seisjl_temp.eventla     = sacdata.evla
    seisjl_temp.eventlo     = sacdata.evlo
    seisjl_temp.eventdp     = sacdata.evdp
    seisjl_temp.dist        = sacdata.dist
    seisjl_temp.nzyear      = sacdata.nzyear
    seisjl_temp.nzday       = sacdata.nzjday
    seisjl_temp.nzhour      = sacdata.nzhour
    seisjl_temp.nzmin       = sacdata.nzmin
    seisjl_temp.nzsec       = sacdata.nzsec
    seisjl_temp.nzmsec      = sacdata.nzmsec
    seisjl_temp.npts        = sacdata.npts
    seisjl_temp.kstnm       = sacdata.kstnm
    seisjl_temp.kevnm       = sacdata.kevnm
    seisjl_temp.kcmpnm      = sacdata.kcmpnm
    seisjl_temp.x           = sacdata.t

    return seisjl_temp
end

"""
Cumulative Trapezoidal Integration
    cumtrapz(x, y, dim=1)
    :param x: vector describing time samples
    :param y: array describing response
    :param dim: dimension to integrate over
"""
function cumtrapz(x::AbstractArray, y::AbstractArray, dim::Integer=1)

    if ndims(y) == 1
        my = length(y);
    else
        error("error in cumtrapz because of the size of y.")
    end

    dt = diff(x);

    z = zeros(my,1)

    for i = 1:my-1

        z[i+1] = z[i] + 0.5 * dt[i] * (y[i] + y[i+1]);

    end

    return z
end

"""
Convert velocity to displacement with filtering.
    vel2disp_withfiltering(v, fs, freqmin, freqmax, corner=5)

"""
function vel2disp_withfiltering(v1::AbstractArray, t1::AbstractArray, fs::Float64, freqmin::Float64, freqmax::Float64; corner::Int=5, fft_cut_amp::Float64=5e1, IsPlot::Bool=false)

    v1_temp1 = detrend(v1)
    v1_temp2 = bandpass(v1_temp1, freqmin, freqmax, fs, corners=corner)
    d1_temp  = cumtrapz(t1, v1_temp2)

    fft_d1  = fft(d1_temp)
    L       = length(fft_d1)
    freq_d1 = fs .* collect(0:Int(L/2)) ./ L

    #cut low frequency noise
    fft_d1[ abs.(fft_d1/L) .> fft_cut_amp] .= 0.0
    cutID           = findall(x -> x < 0.1*freqmin, freq_d1)
    fft_d1[cutID]  .= 0.0

    if IsPlot
        Pu1_temp        =  abs.(fft_d1/L)
        Pu1_temp        = Pu1_temp[1:Int(L/2 + 1)]
        Pu1[2:end-1]    = 2 .* Pu1[2:end-1]
        
        plot(freq_d1, Pu1, line=(:red, 1, :solid),
        ylabel  = "|P1(f)|", 
        xlabel  = "Frequency [Hz]",
        title   = "Single-sided amplitude spectrum of u1",
        xlim    = (1e-4, 10.0),
        ylim    = (1e-6, 1e5),
        xscale  =:log10,
        yscale  =:log10
        )
        savefig("../fig/spectrum_after_filtering.png")
    end

    d1_temp2 = real(ifft(fft_d1))
    d1_temp3 = detrend(d1_temp2)
    d1       = bandpass(d1_temp3, freqmin, freqmax, fs)

    return d1
end

end