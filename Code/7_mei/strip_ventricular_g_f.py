# A basic test case
# ===========================================================
#
# This test case solves the monodomain equations on a [0,12]*[0,0.01] mm^2 rectangular domain,
# with the Paci2013 cell model, using the splittingsolver with the default 
# parameter values and default initial conditions. 
#
# Joanneke E Jansen, May 2017

# Import the cbcbeat module
from cbcbeat import *
import numpy as np
import bisect

# Turn on FFC/FEniCS optimizations
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["cpp_optimize"] = True
flags = ["-O3", "-ffast-math", "-march=native"]
parameters["form_compiler"]["cpp_optimize_flags"] = " ".join(flags)
parameters["form_compiler"]["quadrature_degree"] = 3

# Turn off adjoint functionality
parameters["adjoint"]["stop_annotating"] = True

# Define a [0,12]*[0,0.01] mm^2 rectangular domain
mesh = RectangleMesh(Point(0.0, 0.0), Point(12.0, 0.01), 100, 1)

# Define averaged nominal conductivities, surface to volume ratio and membrane
# capacitance, as found in Sepulveda, Roth, & Wikswo. (1989), Table 1
sigma_l = 150           # mu S / mm
sigma_t = 20            # mu S / mm
beta = 200.0            # mm^{-1}
C_m = 0.2               # mu F / mm^2 
#Note that the total cell capacitance is 9.87109e-2 mu F in the Paci2013 model, 
#so we assume a cell area of 9.87109e-2/0.2=0.49355 mm^2 TOO LARGE??

# Scale conducitivites by 1/(C_m * chi)
M_l = sigma_l/(C_m*beta) # mm^2 / s
M_t = sigma_t/(C_m*beta) # mm^2 / s

# Define the conductivity (tensor)
M = as_tensor(((M_l, 0.0), (0.0, M_t)))

# Define time
time = Constant(0.0)

# Define the external stimulus
stimulus = Expression('(x[0] < 3.0 && (t-floor(t)) < 0.005 ? v_amp : 0.0)', v_amp=5.5/0.987109, t=time, degree=1)

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
  T = 5.0   # Final time
  interval = (0.0, T)

  # Create files to save v and [CA]i
  times = []
  v_left = []
  v_right = []
  cai_left = []
  cai_right = []
  i=0

  # Solve
  for (timestep, fields) in solver.solve(interval, h):
      # Extract the components of the field (vs_ at previous timestep,
      # current vs, current vur)
      (vs_, vs, vur) = fields
      print "(t_0, t_1) = ", timestep
      print "v(5.0,0.0) = ", vs((5.0,0.0))[0]
      print "v(10.0,0.0) = ", vs((10.0,0.0))[0]
      
      #Record vs at each time step or at each ms
      #if (abs(round(timestep[1])-timestep[1])) < (0.1*h):
      if True:
        times.append(timestep[1])
        v_left.append(vs((5.0,0.0))[0])
        v_right.append(vs((10.0,0.0))[0])
        cai_left.append(vs((5.0,0.0))[16])
        cai_right.append(vs((10.0,0.0))[16])
        i=i+1
        #print i
  return times, v_left, v_right, cai_left, cai_right

for fact in np.array([0.5,1.5]):
    print fact
    model = Paci2013_ventricular()
    model.set_parameters(g_f=fact*params["g_f"])
    vsinit = initialmodel(model)
    model.set_initial_conditions(V=vsinit.split(deepcopy=True)[0], m=vsinit.split(deepcopy=True)[1], h=vsinit.split(deepcopy=True)[2], \
      j=vsinit.split(deepcopy=True)[3], d=vsinit.split(deepcopy=True)[4], f1=vsinit.split(deepcopy=True)[5], f2=vsinit.split(deepcopy=True)[6], \
      fCa=vsinit.split(deepcopy=True)[7], Xr1=vsinit.split(deepcopy=True)[8], Xr2=vsinit.split(deepcopy=True)[9], Xs=vsinit.split(deepcopy=True)[10], \
      Xf=vsinit.split(deepcopy=True)[11], q=vsinit.split(deepcopy=True)[12], r=vsinit.split(deepcopy=True)[13], Nai=vsinit.split(deepcopy=True)[14], \
      g=vsinit.split(deepcopy=True)[15], Cai=vsinit.split(deepcopy=True)[16], Ca_SR=vsinit.split(deepcopy=True)[17])
    (times, v_left, v_right, cai_left, cai_right) = main(model)
    np.savetxt("recorded_times_g_f_%s.txt" % fact, times)
    np.savetxt("v_left_g_f_%s.txt" % fact, v_left)
    np.savetxt("v_right_g_f_%s.txt" % fact, v_right)
    np.savetxt("cai_left_g_f_%s.txt" % fact, cai_left)
    np.savetxt("cai_right_g_f_%s.txt" % fact, cai_right)

print "Success!"
