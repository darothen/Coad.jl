
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
    n = 6,
    kernel = :long,
)

# Set up initial conditions
r̅=(10.0)*1e-6
dist = ExponentialDist(L̅=(1.0)*1e-3, x̅=mass_from_r(r̅))
set!(model, dist)

println("\n Plotting initial conditions...")
p = plot(
    model.rᵢ*1e6, model.gᵢ*1e3, label="t = 0 min", 
    xaxis=:log, xlabel="r (μm)",
    # xlim=(0.5, 5000), xticks=[1, 10, 100, 1000],
    # ylim=(0, 0.9), yticks=0:0.06:0.9,
    ylabel="g (g / m³)",
)
display(p)

t = Timestamp(0)
tmax = 60*60 + 1 # seconds
Δt = 5.0 # s
Δt_diag = 5 * 60
Δt_plot = 10 * 60
nt = ceil(Integer, tmax / Δt)

is_diag_step(t; Δt_diag=Δt_diag) = t % Δt_diag == Timestamp(0)
is_plot_step(t; Δt_plot=Δt_plot) = t % Δt_plot == Timestamp(0)

###
# using Profile
# using ProfileSVG
# Profile.clear()
# @profile step!(model, Δt)
# ProfileSVG.save(joinpath("assets", "prof.svg"))
###

println("\n Beginning simulation...")
println(" -----------------------")

g_diags = []
for i in 1:nt
    global t += Δt
    step!(model, Δt)

    mass_tot = sum(model.gᵢ)
    @printf "t = %s | mass %10.3e\n" t mass_tot

    if is_diag_step(t)
        push!(g_diags, model.gᵢ[:])
    end

    if is_plot_step(t)
        display(plot!(p, model.rᵢ*1e6, model.gᵢ*1e3, label = "t = $(t.minutes) min"))
    end
end

println("End; press any key to close.")
# using NPZ
# for (i, g) in enumerate(g_plots)
#     file_fn = @sprintf "output_%03d.npy" i
#     @printf "%03d %s %10.3e \n" i file_fn maximum(g)
#     file_pth = joinpath("output", file_fn)
#     npzwrite(file_pth, g)
# end
# npzwrite("output/r_grid.npy", model.rᵢ)
# xxx = readline()

# using NPZ
# for (i, g) in enumerate(g_diags)
#     file_fn = @sprintf "output_%03d.npy" i
#     @printf "%03d %s %10.3e \n" i file_fn maximum(g)
#     file_pth = joinpath("output", file_fn)
#     npzwrite(file_pth, g)
# end
# npzwrite("output/r_grid.npy", model.rᵢ)


# Set up integrator
# solver = Integrator1D(model, Δt=5.0, tmax=3601)

# Run
# run!(solver)

# Do something with the outputs?