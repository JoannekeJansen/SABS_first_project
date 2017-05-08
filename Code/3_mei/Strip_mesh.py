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

# Set to 'True' to create files that can be opened in Paraview
Postprocessing=False

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

# Use the Grandi cell model
cell_model = Grandi_pasqualini_bers_2010_adjusted()
cell_model.set_parameters(g_CaL_factor=1.0)

# Define the external stimulus
stimulus = Expression('(x[0] < 0.5 && t < 10.0 ? p : 0)', p=10.0, t=time, degree=1)

# Collect this information into the CardiacModel class
cardiac_model = CardiacModel(mesh, time, M, 'none', cell_model, stimulus)

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
vs_.assign(cell_model.initial_conditions(), solver.VS)
vs.assign(cell_model.initial_conditions())

# Set time stepping parameters
h = 0.01 # Time step size
T = 500.00   # Final time
interval = (0.0, T)

# Create HDF5 file to save vs
#hdf_v = HDF5File(mesh.mpi_comm(), "Basic_test_case_monodomain_synthetic_observations_v_adjusted_mini.h5", "w")
#hdf_cai = HDF5File(mesh.mpi_comm(), "Basic_test_case_monodomain_synthetic_observations_cai_adjusted_mini.h5", "w")
times = []
v_left = []
v_right = []
cai_left = []
cai_right = []
i=0

#timer = Timer("xxx:Solve")

# Solve
for (timestep, fields) in solver.solve(interval, h):
    # Extract the components of the field (vs_ at previous timestep,
    # current vs, current vur)
    (vs_, vs, vur) = fields
    #print "(t_0, t_1) = ", timestep
    #print "v(10.0,0.0) = ", vs((10.0,0.0))[0]
    #print "v(0.0,0.0) = ", vs((0.0,0.0))[0]
    
    #Record vs at each time step or at each ms
    if (abs(round(timestep[1])-timestep[1])) < (0.1*h):
      #hdf_v.write(vs.split(deepcopy=True)[0],"v",i)
      #hdf_cai.write(vs.split(deepcopy=True)[36],"cai",i)
      times.append(timestep[1])
      v_left.append(vs((0.0,0.0))[0])
      v_right.append(vs((10.0,0.0))[0])
      cai_left.append(vs((0.0,0.0))[36])
      cai_right.append(vs((10.0,0.0))[36])
      i=i+1
      print i
#del hdf_v, hdf_cai

# Save time steps
np.savetxt("recorded_times.txt", times)

# Save v and [Ca]_i at point in centre and point in corner
np.savetxt("v_left_mesh_100.txt", v_left)
np.savetxt("v_right_mesh_100.txt", v_right)
np.savetxt("cai_left_mesh_100.txt", cai_left)
np.savetxt("cai_right_mesh_100.txt", cai_right)

#timer.stop()
#list_timings(TimingClear_keep, [TimingType_wall])

# Postprocess recorded_times to be able to plot certain times in Paraview
if Postprocessing == True:
    times_to_plot=[0.05, 0.1] # Pick times to plot
    hdf_vs = HDF5File(mesh.mpi_comm(), "Results/Basic_test_case_monodomain_synthetic_observations_vs.h5", "r")
    xf_v = XDMFFile(mesh.mpi_comm(), "Results/Basic_test_case_monodomain_synthetic_observations_v_plot.xdmf")
    xf_cai = XDMFFile(mesh.mpi_comm(), "Results/Basic_test_case_monodomain_synthetic_observations_cai_plot.xdmf")
    attr_vs = hdf_vs.attributes("vs")
    N = attr_vs['count']
    Q = vs.function_space()
    vs_obs = Function(Q)
    for i in range(len(times_to_plot)):
      # Find closest value of times_to_plot[i] to one of the recorded times
      pos=bisect.bisect_left(times, times_to_plot[i])
      before = times[pos - 1]
      after = times[pos]
      if after - times_to_plot[i] < times_to_plot[i] - before:
        idx = pos
      else:
        idx = pos-1
      # Retrieve respective values of v and [Ca]_i at that time step and save
      dataset_vs = "vs/vector_%d"%idx
      hdf_vs.read(vs_obs, dataset_vs)
      xf_v.write(vs_obs.split(deepcopy=True)[0], times[idx])
      xf_cai.write(vs_obs.split(deepcopy=True)[36], times[idx])
    del hdf_vs

print "Success!"
