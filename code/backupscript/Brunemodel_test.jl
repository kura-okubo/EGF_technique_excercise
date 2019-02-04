"""
Brune's model
02/03/2019 Kurama Okubo
"""

using Plots

#-------------------------------------#
#A priori Parameters

#Medium info
dt = 0.1

mu = 30e9


Risetime = 10
slipmax = 1

#-------------------------------------#


function get_u(t, slipmax, Risetime)

	for i = 1:length(t)
		u[i] = slipmax * (1-exp(-t[i]/(Risetime/4)))
	end

	return u
end

Tidmax = round(Int, Risetime*2/dt)
t = dt * collect(0:Tidmax)

u = zeros(Float64, length(t))

u = get_u(t, slipmax, Risetime)

p1 = plot(vcat(-dt*collect(1:10)[end:-1:1], t), vcat(zeros(10,1), u), line=(:black, 1, :solid),
    ylabel  = "displacement", 
    xlabel  = "Time [s]",
    label="Brune's Model"
    )

p1 = plot!([Risetime, Risetime], [0, 1.2*slipmax], line=(:brue, 1, :dash),
    label="Input Rise time"
    )
