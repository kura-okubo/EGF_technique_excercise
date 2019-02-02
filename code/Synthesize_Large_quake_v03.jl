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

t_green = load("../dataset/egfdata/egfdata.jld", "t_green")
t_large = load("../dataset/egfdata/egfdata.jld", "t_large")
disp_green = load("../dataset/egfdata/egfdata.jld", "disp_green")
disp_large = load("../dataset/egfdata/egfdata.jld", "disp_large")


