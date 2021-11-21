
# ```julia
# using Coad
# pkg"add Coad, Plots"
# ```
push!(LOAD_PATH, pwd())

using Coad
using Plots
using Printf

# Set up grid / model simulation
model = Coad1D(
    n = 4 ,
    kernel = :long,
)

# Set up initial conditions
r̅=(10.0)*1e-6
dist = ExponentialDist((1.0)*1e-3, mass_from_r(r̅),)
set!(model, dist)

p = plot(
    model.rᵢ*1e6, model.gᵢ*1e3, label="t = 0 min", 
    xaxis=:log, xlabel="r (μm)",
    # xlim=(0.5, 5000), xticks=[1, 10, 100, 1000],
    # ylim=(0, 0.9), yticks=0:0.06:0.9,
    ylabel="g (g / m³)",
)
display(p)

tmax = 60*30 + 1 # seconds
Δt = 5.0 # s
Δt_plot = 10 # minutes
nt = ceil(Integer, tmax / Δt)
for i in 1:nt
    step!(model, Δt)
end

lmin = floor(Integer, nt*Δt/60)
display(plot!(p, model.rᵢ*1e6, model.gᵢ*1e3, 
              label = "t = $lmin min"))
println("End; press any key to close.")
xxx = readline()



# Set up integrator
# solver = Integrator1D(model, Δt=5.0, tmax=3601)

# Run
# run!(solver)

# Do something with the outputs?