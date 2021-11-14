
const ρw = 1000.0 # Density of water, kg/m³
const golovin_b = (1500.0)*1e-3 # Golovin collision coefficient, input cm³/g/s to -> m³/kg/s


mass_from_r(r; ρ=ρw) = (4*π/3)*ρ*(r^3)
r_from_mass(x; ρ=ρw) = (3*x/4/π/ρ)^(1/3)

# Initial cloud droplet distribution
nc(x; L=L, x̅=x̅) = (L / x̅^2) * exp(-x / x̅)

# Collision Kernels

#   inputs in kg
golovin_kernel(xᵢ, xⱼ) = golovin_b * (xᵢ + xⱼ)
