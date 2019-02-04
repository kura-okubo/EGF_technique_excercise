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
NW = 5  # number of element area in the direction of W

totalfaultarea = 10e6 # a priori total slip area [m^2]

#Rupture info
vr = 0.6*cs # rupture velocity [m/s]
T_rise = 2.0  # rise time [s]
dtslip = 0.1 # time step of slipfunction [s]

slipfunc = "heaviside" # "heaviside", "brune"

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

sfault, sfaultinfo = faultmodel(totalfaultarea, mu, NL, NW, Faultaspectratio, m0, M0, slipmodel = "homogeneous")

foinfo = open("../dataset/faultmodelinfo.dat", "w")
write(foinfo, "L [km], W [km], A [km^2], Slip Average [m], NumofElement, deltaL, deltaW\n")
str = string([sfaultinfo.L/1e3, sfaultinfo.W/1e3, sfaultinfo.A/1e6, sfaultinfo.AverageSlip,
	 sfaultinfo.NumofElement, sfault[1,1].deltaL, sfault[1,1].deltaW])
write(foinfo, str)
close(foinfo)

#assign nucleation patch
Nucleation_x = sfault[round(Int,NL/2), round(Int,NW/2)].ex
Nucleation_z = sfault[round(Int,NL/2), round(Int,NW/2)].ez


t_syn = t_large
u_syn = zeros(Float64, length(t_syn),1)

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

fid = jldopen("../dataset/egfdata/syndata.jld", "w")
fid["t_small"] = t_green
fid["u_small"] = disp_green
fid["t_large"] = t_large
fid["u_syn"] = u_syn
fid["u_ovs"] = disp_large
close(fid)
#Plot

"""
p1 = plot(t_green, disp_green*1e-7, xlabel = "Time [sec]",
                    ylabel = "cm",
                    linecolor = :black,
                    title = "Green's event"
                    )

p2 = plot(t_large, disp_large*1e-7, xlabel = "Time [sec]",
                    ylabel = "cm",
                    linecolor = :red,
                    title = "Observation Large event"
                    )

p2 = plot!(t_syn, u_syn*1e-7, xlabel = "Time [sec]",
                    ylabel = "cm",
                    linecolor = :blue,
                    title = "Synthetic Large event"
                    )

plot(p1, p2, layout = (2,1),
                 size = (800, 1200))

"""

