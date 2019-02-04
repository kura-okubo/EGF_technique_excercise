"""
Empirical Green's function exercise
02/01/2019 Kurama Okubo
"""

include("./SeisJul_vkurama/SeisJul.jl")
include("./SeisJul_vkurama/correlate_kurama.jl")
include("./SeisJul_vkurama/tools.jl")
include("./SeisJul_vkurama/filter.jl")
include("./SeisJul_vkurama/egf_tools.jl")

using .EGF_TOOLS
using .Correlate
using .Filter

using Plots, FFTW, JLD

#-------------------------------------#
#A priori Parameters

IsFiltering = false #apply filtering on the synthetic data
IsOmega2model = true #apply filtering on the synthetic data

#-------------------------------------#


#Load data
t_small = load("../dataset/egfdata/syndata.jld", "t_small")
t_large = load("../dataset/egfdata/syndata.jld", "t_large")
u_small = load("../dataset/egfdata/syndata.jld", "u_small")
u_syn = load("../dataset/egfdata/syndata.jld", "u_syn")
u_obs = load("../dataset/egfdata/syndata.jld", "u_ovs")

if IsFiltering
    #Filtering u_syn
    println("Filtering synthetic result")

    freqmin = 0.1
    freqmax = 1.0
    fs = 1/(t_large[2] - t_large[1])
    corner = 12

    u_syn_temp = u_syn
    
    u1_temp1 = detrend(u_syn_temp)
    u1_temp2 = bandpass(u1_temp1, freqmin, freqmax, fs, corners=corner)

    u_syn = zeros(length(u_syn), 1)
    u_syn = u1_temp2
end

if IsOmega2model
    #Filtering u_syn
    println("Omega square constraint")

    u_syn_temp = u_syn
    
    fs = 1/(t_large[2]-t_large[1])

    fft_df  = fft(u_syn_temp)
    #shift fftdf

    L       = length(fft_df)
    fft_df_temp1 = vcat(fft_df[round(Int,L/2+1):end],fft_df[1:round(Int,L/2)])
    
    freq_df = fs .* collect(-round(Int, L/2):round(Int, L/2-1)) ./ L
    Puf        = abs.(fft_df_temp1/L)

    plot(freq_df, Puf, 
        ylabel  = "|P1(f)|", 
        xlabel  = "Frequency [Hz]",
        title   = "Amplitude spectrum of u1",
        yscale  =:log10
        )

    #falpha = 2.8e-1
    falpha = 8e-1

    fft_df_temp2 = Array{Complex{Float64}}(undef, length(freq_df))

    for i = 1:length(freq_df)

        if abs(freq_df[i]) > 1.0
            fft_df_temp2[i,1] = (falpha * real(fft_df_temp1[i]) * freq_df[i]^(-2)) + im * (falpha * imag(fft_df_temp1[i]) * freq_df[i]^(-2))
        else
            fft_df_temp2[i,1] = real(fft_df_temp1[i]) + im * imag(fft_df_temp1[i])
        end

    end


    #check result
    """
    pc = plot(freq_df[round(Int, L/2+2):end], real(abs.(fft_df_temp1[round(Int, L/2+2):end])/L),
        line=(:black, 1, :solid),
        ylabel  = "|P1(f)|", 
        xlabel  = "Frequency [Hz]",
        label="Positive freq",
        xscale  =:log10,
        yscale  =:log10
        )

    pc = plot!(abs.(freq_df[1:round(Int, L/2)]),  real(abs.(fft_df_temp1[1:round(Int, L/2)])/L),
        line=(:blue, 1, :dash),
        label   = "Negative freq"
        )
    """

    #put back shift

    fft_df_temp_reshift = vcat(fft_df_temp2[round(Int,L/2+1):end], fft_df_temp2[1:round(Int,L/2)])

    #u_syn1 = zeros(length(u_syn),1)
    u_syn = real(ifft(fft_df_temp_reshift))
    #u_syn1 = real(ifft(fft_df_temp2))

end

plot(u_syn)


p1 = plot(t_small, u_small*1e-7, xlabel = "Time [sec]",
                    ylabel = "cm",
                    linecolor = :black,
                    label = "Small event"
                    )
p2 = plot(t_large, u_syn*1e-7, xlabel = "Time [sec]",
                    ylabel = "cm",
                    linecolor = :red,
                    label="Synthetic",
                    ylim    = (-0.05, 0.05)
                    )

p2 = plot!(t_large, u_obs*1e-7, xlabel = "Time [sec]",
                    ylabel = "cm",
                    linecolor = :blue,
                    linestyle = :solid,
                    label="Observation",
                    ylim    = (-0.05, 0.05)
                    )

#FFT

fs = 1/(t_small[2]-t_small[1])

fft_d1  = fft(u_small)
L       = length(fft_d1)
freq_d1 = fs .* collect(0:Int(L/2)) ./ L
freq_d1[1] = 1e-6

Pu1        = abs.(fft_d1/L)
Pu1        = Pu1[1:Int(L/2 + 1)]
Pu1[2:end-1]    = 2 .* Pu1[2:end-1]

p3 = plot(freq_d1, Pu1, line=(:black, 1, :solid),
    ylabel  = "|P1(f)|", 
    xlabel  = "Frequency [Hz]",
    title   = "Single-sided amplitude spectrum of u1",
    label="small event",
    xlim    = (3e-1, 1e1),
    ylim    = (1e-3, 1e5),
    xscale  =:log10,
    yscale  =:log10,
    aspect_ratio = 1.5e-4
    )

#Observation
fft_d3  = fft(u_obs)
L       = length(fft_d3)
freq_d3 = fs .* collect(0:Int(L/2)) ./ L
freq_d3[1] = 1e-6

Pu3        =  abs.(fft_d3/L)
Pu3        = Pu3[1:Int(L/2 + 1)]
Pu3[2:end-1]    = 2 .* Pu3[2:end-1]


#Synthetic

p3 = plot!(freq_d3, Pu3, line=(:blue, 1, :solid),
    ylabel  = "|P1(f)|", 
    xlabel  = "Frequency [Hz]",
    label="Observation",
    xscale  =:log10,
    yscale  =:log10,
    legend=:topright
    )

fft_d2  = fft(u_syn)
L       = length(fft_d2)
freq_d2 = fs .* collect(0:Int(L/2)) ./ L
freq_d2[1] = 1e-6

Pu2        =  abs.(fft_d2/L)
Pu2        = Pu2[1:Int(L/2 + 1)]
Pu2[2:end-1]    = 2 .* Pu2[2:end-1]

fs = 1/(t_large[2]-t_large[1])

p3 = plot!(freq_d2, Pu2, line=(:red, 1, :solid),
    ylabel  = "|P1(f)|", 
    xlabel  = "Frequency [Hz]",
    label="Synthetic",
    xscale  =:log10,
    yscale  =:log10
    )

#save figure
savefig(p1, "../fig/smallevent.png")
savefig(p2, "../fig/largeevent.png")
savefig(p3, "../fig/comparisonspectra.png")

"""
p4 = plot(freq_d3, Pu2 ./ Pu1, line=(:black, 1, :solid),
    ylabel  = "|P1(f)/P2(f)|", 
    xlabel  = "Frequency [Hz]",
    label="Observation",
    xlim    = (5e-2, 1),
    ylim    = (-10, 1000),
    xscale  =:log10
    )


plot(p1, p3, p2, p4, layout = (2,2),
                 size = (1600, 1600))

#savefig("../fig/comparison.png")

"""
