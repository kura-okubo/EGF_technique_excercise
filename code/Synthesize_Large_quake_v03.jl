"""
Empirical Green's function exercise
02/01/2019 Kurama Okubo
"""

using Plots, FFTW, JLD, Dates

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

#Rupture info
vr = 0.6*cs # rupture velocity [m/s]
T_rise = 2.0  # rise time [s]
dtslip = 0.1 # time step of slipfunction [s]

totalfaultarea = 30e6 # apriori total slip area [m^2]
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

#assign nucleation patch
Nucleation_x = sfault[round(Int,NL/2), round(Int,NW/2)].ex
Nucleation_z = sfault[round(Int,NL/2), round(Int,NW/2)].ez


t_syn = t_large
u_syn = zeros(Float64, length(t_syn),1)


for l = 1:length(t_syn)

	if mod(l,100) == 0
		
		println("Last "*string(length(t_syn)-l))
		println(now())

	end

	t1 = t_syn[l]

	usum = 0

	for i = 1:NL
		for j = 1:NW
			for k = 1:NRT

				#distance from the nucleation to the target element
				rdist = sqrt((sfault[i,j].ex - Nucleation_x)^2 + (sfault[i,j].ez - Nucleation_z)^2)
				#println(rdist)
				tijr = rdist/vr #[s]
				tijks = k*dtslip #[s] 

				T_at_k = t1 - tijr - tijks

				#get signal of small event at T_at_k

				if T_at_k < 0 #if T is negative, no signal contributes from small event
					euijk = 0.0

				else 
					tid = findfirst(x -> x >= T_at_k, t_green)
					euijk = disp_green[tid]
				end

				#get slip at ij and time k
				slip_atk = slipmodel(tijks, T_rise, sfault[i,j].es, slipfunc=slipfunc)	

				usum += (mu * sfault[i,j].eA * slip_atk / m0) * euijk

			end
		end
	end

	u_syn[l] = usum

end

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

