# A basic test case
# ===========================================================
#
# This test case solves the monodomain equations on a [0,5]^2 square domain,
# with the Grandi cell model, using the splittingsolver with the default 
# parameter values and default initial conditions. 

# Import the cbcbeat module
from cbcbeat import *
from mshr import *
import numpy as np

# Set to 'True' to create synthetic obervations, 
# to be used in Basic_test_case_monodomain_inverse_problem.py
Synthetic_observations=True

# Turn on FFC/FEniCS optimizations
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["cpp_optimize"] = True
flags = ["-O3", "-ffast-math", "-march=native"]
parameters["form_compiler"]["cpp_optimize_flags"] = " ".join(flags)
parameters["form_compiler"]["quadrature_degree"] = 3

# Turn off adjoint functionality
parameters["adjoint"]["stop_annotating"] = True

# Define a [0,5]^2 domain
domain = Rectangle(Point(0.0, 0.0), Point(5.0, 5.0))
mesh = generate_mesh(domain, 100)

# Define averaged nominal conductivities, surface to volume ratio and membrane
# capacitance, as found in Sepulveda, Roth, & Wikswo. (1989), Table 1
sigma_l = 0.15			# mS / mm
sigma_t = 0.02			# mS / mm
beta = 200.0			# mm^{-1}
C_m = 0.2				# mu F / mm^2

# Scale conducitivites by 1/(C_m * chi)
M_l = sigma_l/(C_m*beta) # mm^2 / ms
M_t = sigma_t/(C_m*beta) # mm^2 / ms

# Define the conductivity (tensor)
M = as_tensor(((M_l, 0.0), (0.0, M_t)))

# Define time
time = Constant(0.0)

# Use the Grandi cell model
cell_model = Grandi_pasqualini_bers_2010()

# Define the external stimulus
p = Expression('10', degree=1, domain=mesh)
stimulus = Expression('(x[0] > 2.0 and x[0] < 3.0 and x[1] > 2.0 \
	and x[1] < 3.0 and t < 3.0 ? p : 0)', p=p, t=time, degree=1)

# Collect this information into the CardiacModel class
cardiac_model = CardiacModel(mesh, time, M, 'none', cell_model, stimulus)

# Customize and create a splitting solver
ps = SplittingSolver.default_parameters()
ps["theta"] = 0.5                        # Second order splitting scheme
ps["pde_solver"] = "monodomain"          # Use Monodomain model for the PDEs
ps["CardiacODESolver"]["scheme"] = "GRL1" # 1st order Rush-Larsen for the ODEs
ps["MonodomainSolver"]["linear_solver_type"] = "iterative"
ps["MonodomainSolver"]["algorithm"] = "cg"
ps["MonodomainSolver"]["preconditioner"] = "petsc_amg"

solver = SplittingSolver(cardiac_model, params=ps)

# Extract the solution fields and set the initial conditions
(vs_, vs, vur) = solver.solution_fields()
vs_.assign(cell_model.initial_conditions(), solver.VS)
vs.assign(cell_model.initial_conditions())

# Time stepping parameters
h = 1 # Time step size
T = 10.0 # Final time
interval = (0.0, T)

vtkfile_v = File('Results/Basic_test_case_monodomain_v.pvd')
vtkfile_Ca_sl = File('Results/Basic_test_case_monodomain_Ca_sl.pvd')

if Synthetic_observations == True:
    # Create HDF5 file to save vs
    hdf_vs = HDF5File(mesh.mpi_comm(), "Results/Basic_test_case_monodomain_synthetic_observations_vs.h5", "w")
    times = []
    i=0

# Solve
for (timestep, fields) in solver.solve(interval, h):
    print "(t_0, t_1) = ", timestep

    # Extract the components of the field (vs_ at previous timestep,
    # current vs, current vur)
    (vs_, vs, vur) = fields
    vtkfile_v << vs.split(deepcopy=True)[0]
    vtkfile_Ca_sl << vs.split(deepcopy=True)[38]
    # If Synthetic_observations == True, record vs each ms
    if Synthetic_observations == True and (round(timestep[1])-timestep[1]) < (0.1*h):
        hdf_vs.write(vs,"vs",i)
        times.append(timestep[1])
        i=i+1

if Synthetic_observations == True:
    np.savetxt("Results/recorded_times.txt", times)
    del hdf_vs

print "Success!"