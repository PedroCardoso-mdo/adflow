# Import statements
import os
import numpy as np
from mpi4py import MPI
from baseclasses import AeroProblem
from adflow import ADFLOW
import vars_master as vars


# Allocate processor sets
comm = MPI.COMM_WORLD

# Create output folder
if not os.path.exists(vars.output):
    if comm.rank == 0:
        os.mkdir(vars.output)


# ------- CFD Solver -------

# ADFlow set up
aeroOptions = {
    # I/O Parameters
    "gridFile": vars.gridFile,
    "outputDirectory": vars.output,
    "monitorvariables": ["resrho", "cl", "cd"],
    "writeTecplotSurfaceSolution": True,
    # Physics Parameters
    "equationType": vars.equationType,
    # Solver Parameters
    "MGCycle": "sg",
    # ANK Solver Parameters
    "useANKSolver": True,
    "acousticScaleFactor": vars.mach,
    "ANKCharTimeStepType":"VLR",
    # NK Solver Parameters
    "useNKSolver": True,
    "NKSwitchTol": vars.nkswitchtol,
    # Termination Criteria
    "L2Convergence": vars.L2Convergence,
    "nCycles": vars.nCycles,
    # Lift index indication
    "liftindex": 3,
    # Variables to write to solution file
    "surfaceVariables": ['rho','p','vx','vy','vz','cp','mach','cf','cfx','cfy','cfz','yplus','blank'],  
    "volumeVariables": ['vort','vortx','vorty','vortz','cp','mach','macht','eddy','eddyratio','resrho','rhoe','resturb', 'blank'],
    #"restartFile" : "/home/ruben/Thesis/Optimization/analysis/0707_2017_Analysis_coarser_mesh/Data/wing_000_vol.cgns",
}

# Create solver
CFDSolver = ADFLOW(options=aeroOptions)

# Add features for plotting 
CFDSolver.addLiftDistribution(150, "y")
CFDSolver.addSlices("y", np.linspace(0.1, 3.6, 200))     # Alterar o linspace caso se pretendam slices diferentes

# Set up the aerodynamic problem
ap = AeroProblem(name="wing", alpha=1.8, mach=vars.mach, altitude=vars.altitude, areaRef=vars.areaRef, chordRef=vars.chordRef, evalFuncs=["cl", "cd"])


# ------- Evaluate and print functions -------

CFDSolver(ap)   # Solve
funcs = {}      # Evaluate and print
CFDSolver.evalFunctions(ap, funcs)
# Print the evaluated functions
if comm.rank == 0:
    print(funcs)
    
