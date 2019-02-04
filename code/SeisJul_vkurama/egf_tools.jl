module EGF_TOOLS

# cross-correlation module
export sac2seisjl, cumtrapz, vel2disp_withfiltering, faultmodel, slipmodel, get_euijk, closest_index, synthesize_large_quake
using Dates, FFTW, Plots
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
abstract type AbstructFaultmodel end
abstract type AbstructFaultinfo end

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

mutable struct sfault <: AbstructFaultmodel

    ex       ::Float64
    ez       ::Float64
    eA       ::Float64
    es       ::Float64
    deltaL   ::Float64
    deltaW   ::Float64

    sfault() = new()

end

mutable struct sfaultinfo <: AbstructFaultinfo

    L               ::Float64
    W               ::Float64
    A               ::Float64
    AverageSlip     ::Float64
    NumofElement    ::Int

    sfaultinfo() = new()

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
        Pu1       =  abs.(fft_d1/L)
        Pu1       = Pu1[1:Int(L/2 + 1)]
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

"""
Make fault model for large event
    faultmodel(NL, NW, Faultaspectratio, m0, M0)

"""

function faultmodel(A::Float64, mu::Float64, NL::Int, NW::Int, Faultaspectratio::Float64, m0::Float64, M0::Float64; 
    slipmodel::String="homogeneous")

    if slipmodel == "homogeneous"
        #compute average slip on each element fault by M0 = mu*D*S
        #s_avrg = M0/(mu*A*NL*NW)
        s_avrg = M0/(mu*A)

        L = sqrt(Faultaspectratio*A)
        W = L/Faultaspectratio

        deltaL = L/NL
        deltaW = W/NW

        sf = Array{sfault}(undef, NL, NW)

        for i = 1:NL
            for j = 1:NW

                ex = (float(i)-0.5) * deltaL
                ez = (float(j)-0.5) * deltaW
                eA = deltaL * deltaW
                es = s_avrg

                sf[i,j] = sfault()
                sf[i,j].ex = ex
                sf[i,j].ez = ez
                sf[i,j].eA = eA
                sf[i,j].es = es
                sf[i,j].deltaL = deltaL
                sf[i,j].deltaW = deltaW
            end
        end

        sfinfo = sfaultinfo()
        sfinfo.L = L
        sfinfo.W = W
        sfinfo.A = A
        sfinfo.AverageSlip = s_avrg
        sfinfo.NumofElement = NL*NW

        return sf, sfinfo

    else
        error("error: fault model is not defined.")
    end
end

"""
Make slip model for large event
    slipmodel(Trise, maxslip, slipmodel)

"""

function slipmodel(t::Float64, Trise::Float64, maxslip::Float64; slipfunc::String="heaviside")

    if slipfunc == "heaviside"
        
        if t < Trise
            slip = 0
        else
            slip = maxslip
        end

        return slip

    else
        error("error: slip model is not defined.")
    end

end

function get_euijk(T_at_k::Float64, t::AbstractArray, disp::AbstractArray)

    if T_at_k < 0 #if T is negative, no signal contributes from small event
        euijk = 0.0
    else

        #improve searching time 
        tmid = (t[1] + t[end])/2

        if tmid > T_at_k
            tid =  findfirst(x -> x >= T_at_k, @views t)
            euijk = disp[tid]
        else
            tid =  length(t) + 2 - findfirst(x -> x < T_at_k, @views t[end:-1:1])
            euijk = disp[tid]
        end

    end
    
end

"""
Synthsize signal of large event
    synthesize_large_quake(t_syn, NL, NW, NRT, sfault, vr, slipfunc)

"""

function synthesize_large_quake(t_syn::Array, t_green::Array, disp_green::Array, NL::Int64, NW::Int64, NRT::Int64, 
    sfault::Array, vr::Float64, T_rise::Float64, mu::Float64, m0::Float64, dtslip::Float64;
    slipfunc::String="heaviside")


    #assign nucleation patch
    #Assume nucleatin from the middle of fault
    Nucleation_x = sfault[round(Int,NL/2), round(Int,NW/2)].ex
    Nucleation_z = sfault[round(Int,NL/2), round(Int,NW/2)].ez


    u_syn = zeros(Float64, length(t_syn),1)

    #For debug
    elapsedtime_i = Array{Float64,1}(undef,100000)
    elapsedtime_k = Array{Float64,1}(undef,100000)
    elapsedtime_f1 = Array{Float64,1}(undef,100000)
    elapsedtime_f2 = Array{Float64,1}(undef,100000)

    for l = 1:length(t_syn)
    #for l = 1:1500

        if mod(l,1000) == 0
            
            println("Last "*string(length(t_syn)-l))
            println(now())

        end

        t1 = t_syn[l]

        usum = 0

        elapsedtime_i[l] = @elapsed for i = 1:NL

            for j = 1:NW

                elapsedtime_k[l] = @elapsed for k = 1:NRT

                    #distance from the nucleation to the target element
                    rdist = sqrt((sfault[i,j].ex - Nucleation_x)^2 + (sfault[i,j].ez - Nucleation_z)^2)
                    #println(rdist)
                    tijr = rdist/vr #[s]
                    tijks = k*dtslip #[s] 

                    T_at_k = t1 - tijr - tijks

                    #get signal of small event at T_at_k

                    elapsedtime_f1[l] = @elapsed euijk = get_euijk(T_at_k, t_green, disp_green)

                    #get slip at ij and time k
                    elapsedtime_f2[l] = @elapsed slip_atk = slipmodel(tijks, T_rise, sfault[i,j].es, slipfunc=slipfunc) 

                    usum += (mu * sfault[i,j].eA * slip_atk / m0) * euijk

                end
            end
        end

        u_syn[l] = usum

    end

    fo = open("./elapsedtime.dat", "w")

    write(fo, "index, elapsed_i, elapsed_k, elapsed_f1, elapsed_f2 [ms]\n")

    for i1 = 1:length(t_syn)
        if mod(i1,1000) == 0
        
        write(fo, string([i1, elapsedtime_i[i1]*1e6, elapsedtime_k[i1]*1e6, elapsedtime_f1[i1]*1e6, elapsedtime_f2[i1]*1e6])*"\n")

        end
    end

    close(fo)

    return u_syn

end


end