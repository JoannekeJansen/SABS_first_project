# Generate synthetic observations
# ===========================================================
#
# This test case solves the monodomain equations on a [0,10]^2 square,
# with the Grandi cell model, using the splittingsolver with the default 
# parameter values and default initial conditions and saves the data in hdf5 format.

# Import the cbcbeat module
from cbcbeat import *
from mshr import *

# Turn on FFC/FEniCS optimizations
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["cpp_optimize"] = True
flags = ["-O3", "-ffast-math", "-march=native"]
parameters["form_compiler"]["cpp_optimize_flags"] = " ".join(flags)
parameters["form_compiler"]["quadrature_degree"] = 3

# Turn off adjoint functionality
parameters["adjoint"]["stop_annotating"] = True

# Define a [0,10]^2 domain
from mshr import *
domain = Rectangle(Point(0.0, 0.0), Point(10.0, 10.0))
mesh = generate_mesh(domain, 100)
time = Constant(0.0)

# Define the conductivity (tensor)
M_i = as_matrix(((0.255, 0.0), (0.0, 0.0775)))

# We use the Grandi cell model
#cell_model = Grandi_pasqualini_bers_2010()
cell_model = FitzHughNagumoManual()

# Define the external stimulus
p = Expression('25', degree=1, domain=mesh)
stimulus = Expression('(x[0] > 4.5 and x[0] < 5.5 and x[1] > 4.5 \
	and x[1] < 5.5 and t > 0 and t < 100 ? p : 0)', p=p, t=time, degree=1)

# Collect this information into the CardiacModel class
cardiac_model = CardiacModel(mesh, time, M_i, 'none', cell_model, stimulus)

# Customize and create a splitting solver
ps = SplittingSolver.default_parameters()
ps["theta"] = 0.5                        # Second order splitting scheme
ps["pde_solver"] = "monodomain"          # Use Monodomain model for the PDEs
ps["CardiacODESolver"]["scheme"] = "GRL1" #  1st order Rush-Larsen for the ODEs
ps["MonodomainSolver"]["linear_solver_type"] = "iterative"
ps["MonodomainSolver"]["algorithm"] = "cg"
ps["MonodomainSolver"]["preconditioner"] = "petsc_amg"

solver = SplittingSolver(cardiac_model, params=ps)

# Extract the solution fields and set the initial conditions
(vs_, vs, vur) = solver.solution_fields()
vs_.assign(cell_model.initial_conditions())

# Time stepping parameters
dt = 0.5
T = 5.5
interval = (0.0, T)

# Create HDF5 files to save v and c at N time steps
hdf_v = HDF5File(mesh.mpi_comm(), "Results/Basic_test_case_monodomain_synthetic_observations_v.h5", "w")
hdf_c = HDF5File(mesh.mpi_comm(), "Results/Basic_test_case_monodomain_synthetic_observations_c.h5", "w")
i = 1

vtkfile = File('Results/Basic_test_case_monodomain.pvd')

# Solve
for (timestep, fields) in solver.solve(interval, dt):
    # Extract the components of the field (vs_ at previous timestep,
    # current vs, current vur) and save to HDF5 files.
    (vs_, vs, vur) = fields

    # We save v and c every ms 
    if  (round(timestep[0])-timestep[0]) < (0.1*dt):
    	hdf_v.write(vs.split(deepcopy=True)[0], "v", i)
    	hdf_c.write(vs.split(deepcopy=True)[1], "c", i)
    	i = i+1
    	print "Saved u and c at t=", timestep[0]
    vtkfile << vur


# Delete HDF5File object, closing file
del hdf_v, hdf_c

print "Generated N =", i-1, "synthetic observations of u and c!"
