# A basic test case
# ===========================================================
#
# This test case solves the monodomain equations on a [0,5]^2 square domain,
# with the Grandi cell model, using the splittingsolver with the default 
# parameter values and default initial conditions. 
#
# Joanneke E Jansen, April 2017

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

# Define a [0,5]^2 domain
mesh = RectangleMesh(Point(0.0, 0.0), Point(10.0, 0.01), 100, 1)

# Define averaged nominal conductivities, surface to volume ratio and membrane
# capacitance, as found in Sepulveda, Roth, & Wikswo. (1989), Table 1
sigma_l = 0.15          # mS / mm
sigma_t = 0.02          # mS / mm
beta = 200.0            # mm^{-1}
C_m = 0.2               # mu F / mm^2

# Scale conducitivites by 1/(C_m * chi)
M_l = sigma_l/(C_m*beta) # mm^2 / ms
M_t = sigma_t/(C_m*beta) # mm^2 / ms

# Define the conductivity (tensor)
M = as_tensor(((M_l, 0.0), (0.0, M_t)))

# Define time
time = Constant(0.0)

# Define the external stimulus
stimulus = Expression('(x[0] < 0.5 && t < 10.0 ? p : 0)', p=10.0, t=time, degree=1)

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
  vs_.assign(model.initial_conditions(), solver.VS)
  vs.assign(model.initial_conditions())

  # Set time stepping parameters
  h = 0.01 # Time step size
  T = 500.00   # Final time
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
      #print "(t_0, t_1) = ", timestep
      #print "v(10.0,0.0) = ", vs((10.0,0.0))[0]
      #print "v(0.0,0.0) = ", vs((0.0,0.0))[0]
      
      #Record vs at each time step or at each ms
      #if (abs(round(timestep[1])-timestep[1])) < (0.1*h):
      if True:
        times.append(timestep[1])
        v_left.append(vs((2.5,0.0))[0])
        v_right.append(vs((7.5,0.0))[0])
        cai_left.append(vs((2.5,0.0))[36])
        cai_right.append(vs((7.5,0.0))[36])
        i=i+1
        #print i
  return times, v_left, v_right, cai_left, cai_right

for fact in np.linspace(0.1,2,20):
    print fact
    model = Grandi_pasqualini_bers_2010_adjusted()
    model.set_parameters(g_K1_factor=fact)
    (times, v_left, v_right, cai_left, cai_right) = main(model)
    np.savetxt("recorded_times.txt" % fact, times)
    np.savetxt("v_left_g_K1_%s.txt" % fact, v_left)
    np.savetxt("v_right_g_K1_%s.txt" % fact, v_right)
    np.savetxt("cai_left_g_K1_%s.txt" % fact, cai_left)
    np.savetxt("cai_right_g_K1_%s.txt" % fact, cai_right)

print "Success!"
