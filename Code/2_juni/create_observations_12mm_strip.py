# Create observations, to be used in inverse_12mm_strip.py
#
# Joanneke E Jansen, May 2017

# Import the cbcbeat module
from cbcbeat import *
import numpy as np
import bisect
import argparse

# Set parameters
parser = argparse.ArgumentParser()
parser.add_argument("-gnaf", type=float, help="The factor of g_na")
parser.add_argument("-gk1f", type=float, help="The factor of g_k1")
parser.add_argument("-gkrf", type=float, help="The factor of g_kr")
parser.add_argument("-gcalf", type=float, help="The factor of g_cal")
args = parser.parse_args()
factor_gna = args.gnaf
factor_gk1 = args.gk1f
factor_gkr = args.gkrf
factor_gcal = args.gcalf

# If CreateInitial = True, we create initial conditions by letting the model run for 800 s
CreateInitial=False

# Turn on FFC/FEniCS optimizations
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["cpp_optimize"] = True
flags = ["-O3", "-ffast-math", "-march=native"]
parameters["form_compiler"]["cpp_optimize_flags"] = " ".join(flags)
parameters["form_compiler"]["quadrature_degree"] = 3

# Turn off adjoint functionality
parameters["adjoint"]["stop_annotating"] = True

# Define a [0,12]*[0,0.01] mm^2 rectangular domain
mesh = RectangleMesh(Point(0.0, 0.0), Point(12.0, 0.01), 200, 1)

# Define averaged nominal conductivities as found in Sepulveda, Roth, & Wikswo. (1989), Table 1
sigma_t = 20            # mu S / mm
sigma_l = 200           # mu S / mm

# Take total cell capacitance and volume as in the Paci2013 model, where we assume that the 
# total volume is the sum of the intracellular and SR volumes.
C_m = 9.87109e-5       # mu F
Volume = 9.3837e-6     # mm^3

# Scale conducitivites by Volume/C_m
M_t = sigma_t*Volume/C_m # mm^2 / s
M_l = sigma_l*Volume/C_m # mm^2 / s

# Define the conductivity (tensor) 
M = as_tensor(((M_t, 0.0), (0.0, M_l)))

# Define time
time = Constant(0.0)

# Define the external stimulus
stimulus = Expression('(x[0] < 0.5 && (t-floor(t)) < 0.005 ? v_amp : 0.0)', v_amp=5.5/0.987109, t=time, degree=1)

def initialmodel(model):
  # Collect this information into the CardiacModel class
  cardiac_model = CardiacModel(mesh, time, M, 'none', model, stimulus)

  # Customize and create a splitting solver
  ps = SplittingSolver.default_parameters()
  ps["theta"] = 1.0                       # First order Godunov splitting scheme
  #ps["theta"] = 0.5                       # Second order Strang splitting scheme
  ps["pde_solver"] = "monodomain"          # Use Monodomain model for the PDEs
  ps["CardiacODESolver"]["scheme"] = "GRL1" # 1st order Rush-Larsen for the ODEs
  #ps["MonodomainSolver"]["linear_solver_type"] = "direct"
  ps["MonodomainSolver"]["linear_solver_type"] = "iterative"
  ps["MonodomainSolver"]["algorithm"] = "cg"
  ps["MonodomainSolver"]["preconditioner"] = "petsc_amg"

  solver = SplittingSolver(cardiac_model, params=ps)

  # Extract the solution fields and set the initial conditions
  (vs_, vs, vur) = solver.solution_fields()
  vs_.assign(model.initial_conditions())

  # Set time stepping parameters
  h = 0.0001 # Time step size
  T = 800.0  # Final time
  interval = (0.0, T)

  # Solve
  for (timestep, fields) in solver.solve(interval, h):
      # Extract the components of the field (vs_ at previous timestep,
      # current vs, current vur)
      (vs_, vs, vur) = fields

  return vs

def main(model):
  # Collect the model information into the CardiacModel class
  cardiac_model = CardiacModel(mesh, time, M, 'none', model, stimulus)

  # Customize and create a splitting solver
  ps = SplittingSolver.default_parameters()
  ps["theta"] = 1.0                       # First order Godunov splitting scheme
  #ps["theta"] = 0.5                       # Second order Strang splitting scheme
  ps["pde_solver"] = "monodomain"          # Use Monodomain model for the PDEs
  ps["CardiacODESolver"]["scheme"] = "GRL1" # 1st order Rush-Larsen for the ODEs
  #ps["MonodomainSolver"]["linear_solver_type"] = "direct"
  ps["MonodomainSolver"]["linear_solver_type"] = "iterative"
  ps["MonodomainSolver"]["algorithm"] = "cg"
  ps["MonodomainSolver"]["preconditioner"] = "petsc_amg"

  solver = SplittingSolver(cardiac_model, params=ps)

  # Extract the solution fields and set the initial conditions
  (vs_, vs, vur) = solver.solution_fields()
  vs_.assign(model.initial_conditions())

  # Set time stepping parameters
  h = 0.0001 # Time step size
  T = 1.0    # Final time
  interval = (0.0, T)

  # Create files to save observed v and [CA]i
  hdf_v = HDF5File(mesh.mpi_comm(), "observed_v_{0}_{1}_{2}_{3}_12mm_strip.h5".format(factor_gna,factor_gk1,factor_gkr,factor_gcal), "w")
  hdf_cai = HDF5File(mesh.mpi_comm(), "observed_cai_{0}_{1}_{2}_{3}_12mm_strip.h5".format(factor_gna,factor_gk1,factor_gkr,factor_gcal), "w")
  times = []
  i=0

  # Solve
  for (timestep, fields) in solver.solve(interval, h):
      # Extract the components of the field (vs_ at previous timestep,
      # current vs, current vur)
      (vs_, vs, vur) = fields
      #Record vs at each ms
      if (abs(round(1000*timestep[1])-1000*timestep[1])) < (100*h):
        hdf_v.write(vs.split(deepcopy=True)[0],"v",i)
        hdf_cai.write(vs.split(deepcopy=True)[16],"cai",i)
        times.append(timestep[1])
        i=i+1
        print timestep[1]
        print "v(7.5,0.0) = ", vs((7.5,0.0))[0]
        #print i
  del hdf_v, hdf_cai
  
  # Save time steps
  np.savetxt("observed_times_{0}_{1}_{2}_{3}_12mm_strip.txt".format(factor_gna,factor_gk1,factor_gkr,factor_gcal), times)

# Record
if True:
    model = Paci2013_ventricular_adjusted()

    # Set parameters
    model.set_parameters(g_Na_factor=factor_gna)
    model.set_parameters(g_CaL_factor=factor_gcal)
    model.set_parameters(g_K1_factor=factor_gk1)
    model.set_parameters(g_Kr_factor=factor_gkr)

    if CreateInitial == True:
      # Create HDF5 file to save initial conditions
      hdf_vs = HDF5File(mesh.mpi_comm(), "initial_conditions_12mm_strip.h5", "w")
      vsinit = initialmodel(model)
      hdf_vs.write(vsinit,"vs",0)
      del hdf_vs 

    # Load initial conditions
    Q = VectorFunctionSpace(mesh, "CG", 1, 18)
    vsinit = Function(Q)
    hdf_vs = HDF5File(mesh.mpi_comm(), "initial_conditions_12mm_strip.h5", "r")
    hdf_vs.read(vsinit, "vs/vector_0")
    del hdf_vs

    # Set initial conditions
    model.set_initial_conditions(V=vsinit.split(deepcopy=True)[0], m=vsinit.split(deepcopy=True)[1], h=vsinit.split(deepcopy=True)[2], \
      j=vsinit.split(deepcopy=True)[3], d=vsinit.split(deepcopy=True)[4], f1=vsinit.split(deepcopy=True)[5], f2=vsinit.split(deepcopy=True)[6], \
      fCa=vsinit.split(deepcopy=True)[7], Xr1=vsinit.split(deepcopy=True)[8], Xr2=vsinit.split(deepcopy=True)[9], Xs=vsinit.split(deepcopy=True)[10], \
      Xf=vsinit.split(deepcopy=True)[11], q=vsinit.split(deepcopy=True)[12], r=vsinit.split(deepcopy=True)[13], Nai=vsinit.split(deepcopy=True)[14], \
      g=vsinit.split(deepcopy=True)[15], Cai=vsinit.split(deepcopy=True)[16], Ca_SR=vsinit.split(deepcopy=True)[17])
    main(model)

print "Success!"
