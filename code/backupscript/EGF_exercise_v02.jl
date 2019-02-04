"""
Empirical Green's function exercise
02/01/2019 Kurama Okubo
"""

using SAC, Plots, FFTW, HDF5, JLD

include("./SeisJul_vkurama/SeisJul.jl")
include("./SeisJul_vkurama/correlate_kurama.jl")
include("./SeisJul_vkurama/tools.jl")
include("./SeisJul_vkurama/egf_tools.jl")

using .EGF_TOOLS

#-------------------------------------#
#A priori Parameters
name_dir_greenevent = "Mw4.3"
name_dir_largeevent = "Mw6.2"

comparison_component = "Z"

finame_greenevent   = "IU.COR..BH"*comparison_component*".M.1991.230.110726.SAC"
finame_largeevent   = "IU.COR..BH"*comparison_component*".M.1991.229.192945.SAC"
#-------------------------------------#

#Load dataset

greenevent_sacpath = "../dataset/"*name_dir_greenevent*"/"*finame_greenevent 
largeevent_sacpath = "../dataset/"*name_dir_largeevent*"/"*finame_largeevent 

function load_SAC_data(name::String)
    d = SAC.read(name)
end


d_green_sac = load_SAC_data(greenevent_sacpath)
d_large_sac = load_SAC_data(largeevent_sacpath)


d_green = sac2seisjl(d_green_sac)
d_large = sac2seisjl(d_large_sac)

t_green = d_green.delta .* collect(1:d_green.npts)
t_large = d_large.delta .* collect(1:d_large.npts)


#Filtering data
vel_raw_green = d_green.x
vel_raw_large = d_large.x

freqmin = 0.1    
freqmax = 5.0 

fs_green = 1/d_green.delta
fs_large = 1/d_large.delta

disp_filtered_green = vel2disp_withfiltering(vel_raw_green, t_green, fs_green, freqmin, freqmax; corner=5, fft_cut_amp=5e1, IsPlot=false)
disp_filtered_large = vel2disp_withfiltering(vel_raw_large, t_large, fs_large, freqmin, freqmax; corner=5, fft_cut_amp=5e5, IsPlot=false)

p1 = plot(t_green, disp_filtered_green*1e-7, xlabel = "Time [sec]",
                    ylabel = "cm",
                    linecolor = :black,
                    title = "Green's event"
                    )

p2 = plot(t_large, disp_filtered_large*1e-7, xlabel = "Time [sec]",
                    ylabel = "cm",
                    linecolor = :red,
                    title = "Large event"
                    )


plot(p1, p2, layout = (2,1),
                 size = (800, 1200))

savefig("../fig/displacement_comparison.png")

"""
NOT WORK WITH STRUCTURE
fid = h5open("../dataset/egfdata/egfdata.h5", "w")
fid["d_green"] = d_green
fid["d_large"] = d_large
fid["disp_green"] = disp_filtered_green
fid["disp_large"] = disp_filtered_large
close(fid)
"""
save("../dataset/egfdata/egfdata.jld", "d_green", d_green,
    "d_large", d_large, "disp_green", disp_filtered_green,
    "disp_large", disp_filtered_large)

