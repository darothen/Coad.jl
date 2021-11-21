
const r₁ = (0.1)*1e-6 # minimum droplet size bin, micron->m
const rm = (100_000.0)*1e-6 # maximum droplet size bin, micron->m
# NOTE: We use subscript "1" because Julia is a 1-indexed language. By convention
# .     we reserve the nought subscript for general parameters.
const gmin = 1e-60 # lower bound on permissible bin mass density

