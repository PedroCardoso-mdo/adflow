## Master file to store variables for optimization runs

# Inputs and outputs
output = "output"
gridFile = "vmesh_L2.cgns"
FFDFile = "fitted_ffdbox.xyz"

# Aero problem
CLStar = 0.8932
initialAoA = 1.893052821
mach = 0.08164
altitude = 304.8
areaRef = 2.1691
chordRef = 0.6091

# FFD deformations
xFraction = 0.25		# quarter-chord

# ADFlow CFD options
equationType = "RANS"
nkswitchtol = 1e-6
L2Convergence = 1e-6
nCycles = 10000

# ADFlow adjoint options
adjointMaxIter = 650
adjointL2Convergence = 1e-8

# Optimizer options
optAcc = 1e-8
optMaxIt = 500
