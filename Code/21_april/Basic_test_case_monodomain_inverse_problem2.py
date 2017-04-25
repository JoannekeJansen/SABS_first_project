# Solver for the inverse problem of Basic_test_case_monodomain.py
# ===========================================================
#

# Import the cbcbeat module
from cbcbeat import *
from mshr import *
import numpy.random

# Turn on FFC/FEniCS optimizations
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["cpp_optimize"] = True
flags = ["-O3", "-ffast-math", "-march=native"]
parameters["form_compiler"]["cpp_optimize_flags"] = " ".join(flags)
parameters["form_compiler"]["quadrature_degree"] = 3

# Define a [0,10]^2 domain
from mshr import *
domain = Rectangle(Point(0.0, 0.0), Point(10.0, 10.0))
mesh = generate_mesh(domain, 100)

# Define the conductivity (tensor)
M=0.255

def forward(GNa):
	# Define time
	time = Constant(0.0)
	
	# We use the Grandi cell model
	cell_model = Grandi_pasqualini_bers_2010()
	cell_model.set_parameters(GNa=GNa)

	# Define the external stimulus
	p = Expression('0', degree=1, domain=mesh)
	stimulus = Expression('(x[0] > 4.5 and x[0] < 5.5 and x[1] > 4.5 \
		and x[1] < 5.5 and t > 0 and t < 100 ? p : 0)', p=p, t=time, degree=1)

	# Collect this information into the CardiacModel class
	cardiac_model = CardiacModel(mesh, time, M, 'none', cell_model)
	#cardiac_model = CardiacModel(mesh, time, M, 'none', cell_model, stimulus)

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

	# Solve forward problem
	for (timestep, fields) in solver.solve(interval, dt):
	    print "(t_0, t_1) = ", timestep

	    # Extract the components of the field (vs_ at previous timestep,
	    # current vs, current vur)
	    (vs_, vs, vur) = fields
	#    v.assign(vur)

	return vur

# # Load synthetic data
# Q = FunctionSpace(mesh, "CG", 1)
# h = Function(Q)
# g = Function(Q)
# j = Function(Q)
# hdf_v = HDF5File(mesh.mpi_comm(), "Results/Basic_test_case_monodomain_synthetic_observations_v.h5", "r")
# hdf_c = HDF5File(mesh.mpi_comm(), "Results/Basic_test_case_monodomain_synthetic_observations_c.h5", "r")
# attr_v = hdf_v.attributes("v")
# attr_c = hdf_c.attributes("c")
# N = attr_v['count']
# for timestep in range(N):
#     dataset_v = "v/vector_%d"%timestep
#     dataset_c = "c/vector_%d"%timestep
#     print 'Retrieving time step t=', timestep
#     hdf_v.read(h, dataset_v)
#     g.assign(h)
#     hdf_c.read(h, dataset_c)
#     j.assign(h)
# del hdf_v, hdf_c

if __name__ == "__main__":
    GNa = Constant(23)	      	  # initial guess, default = 23.
    v = forward(GNa)              # solve the forward problem once. 

    # Define functional of interest
    J = Functional(inner(v, v)*dx*dt[FINISH_TIME])

    # Indicate the control parameter of interest
    dJdnu = compute_gradient(J, Control(GNa))
    Jnu = assemble(inner(v, v)*dx) # current value

    parameters["adjoint"]["test_derivative"] = True
    parameters["adjoint"]["stop_annotating"] = True # stop registering equations

    def Jhat(GNa): # the functional as a pure function of GNa
         v = forward(GNa)
         return assemble(inner(v, v)*dx)

    conv_rate = taylor_test(Jhat, Control(GNa), Jnu, dJdnu, seed=(0.0001))

# Define the reduced functional and solve the optimisation problem:
#rf = ReducedFunctional(J, GNa)
#opt_ctrls = minimize(rf, options={"maxiter": 50})

print "Success!"