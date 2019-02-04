"""
Empirical Green's function exercise
02/01/2019 Kurama Okubo
"""

using Plots, FFTW, JLD, Dates, BenchmarkTools

include("./SeisJul_vkurama/SeisJul.jl")
include("./SeisJul_vkurama/correlate_kurama.jl")
include("./SeisJul_vkurama/tools.jl")
include("./SeisJul_vkurama/egf_tools.jl")

using .EGF_TOOLS

#-------------------------------------#
#A priori Parameters

#Medium info
mu = 30e9 # Shear modulus [Pa]
cs = 3000 # shear wave velocity [m/s]

#Event info
mw0 = 4.3 # Moment magnitude of Small event
Mw0 = 6.2 # Moment magnitude of Large event

mdepth = 15.9 #[km]
mazimuth = 160
mdip = 53
mrake = 106

Mdepth = 12.0 #[km]
Mazimuth = 165
Mdip = 45
Mrake = 100

Faultaspectratio = 2.0 # L(length)/W(width)
NL = 11 # number of element area in the direction of L
NW = 6  # number of element area in the direction of W

totalfaultarea = 60e6 # a priori total slip area [m^2]

#Rupture info
vr = 0.6*cs # rupture velocity [m/s]
T_rise = 2.0  # rise time [s]
dtslip = 0.1 # time step of slipfunction [s]

slipfunc = "heaviside" # "heaviside"
#-------------------------------------#

#Load data
d_green = load("../dataset/egfdata/egfdata.jld", "d_green")
d_large = load("../dataset/egfdata/egfdata.jld", "d_large")
t_green = load("../dataset/egfdata/egfdata.jld", "t_green")
t_large = load("../dataset/egfdata/egfdata.jld", "t_large")
disp_green = load("../dataset/egfdata/egfdata.jld", "disp_green")
disp_large = load("../dataset/egfdata/egfdata.jld", "disp_large")

#compute synthetic large earthquake
#Assumptions:
#1. rise time is uniform on the all elements
#2. shear modulus is uniform on the all elements
#3. Homogeneous slip model (average slip is uniform on the all elements)
#4. Rupture is nucleated at the center of the fault

NRT = round(Int, T_rise / dtslip)

u_syn = zeros(length(disp_large), 1) 
t_syn = t_large

#Make fault model
#get the location and area of element fault 
m0 = 10^(1.5*(mw0+10.7)) * 1e-7 #[J]
M0 = 10^(1.5*(Mw0+10.7)) * 1e-7 #[J]

println(string([m0, M0]))

sfault, sfaultinfo = faultmodel(totalfaultarea, mu, NL, NW, Faultaspectratio, m0, M0, slipmodel = "homogeneous")

foinfo = open("../dataset/faultmodelinfo.dat", "w")
write(foinfo, "L [km], W [km], A [km^2], Slip Average [m], NumofElement, deltaL, deltaW\n")
str = string([sfaultinfo.L/1e3, sfaultinfo.W/1e3, sfaultinfo.A/1e6, sfaultinfo.AverageSlip,
	 sfaultinfo.NumofElement, sfault[1,1].deltaL, sfault[1,1].deltaW])
write(foinfo, str)
close(foinfo)

t_syn = t_large

u_syn = synthesize_large_quake(t_syn, t_green, disp_green, NL, NW, NRT, sfault, vr, T_rise, mu, m0, dtslip, slipfunc=slipfunc)

# save data
fid = jldopen("../dataset/egfdata/syndata.jld", "w")
fid["t_small"] = t_green
fid["u_small"] = disp_green
fid["t_large"] = t_large
fid["u_syn"] = u_syn
fid["u_ovs"] = disp_large
close(fid)


