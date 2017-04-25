# Solver for the inverse problem of Basic_test_case_monodomain.py
# ===========================================================
#

# Import the cbcbeat module
from cbcbeat import *
from mshr import *

# Turn on FFC/FEniCS optimizations
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["cpp_optimize"] = True
flags = ["-O3", "-ffast-math", "-march=native"]
parameters["form_compiler"]["cpp_optimize_flags"] = " ".join(flags)
parameters["form_compiler"]["quadrature_degree"] = 3

# Define a [0,10]^2 domain
from mshr import *
domain = Rectangle(Point(0.0, 0.0), Point(5.0, 5.0))
mesh = generate_mesh(domain, 100)

# Averaged nominal conductivities, surface to volume ratio and membrane
# capacitance, as found in Sepulveda, Roth, & Wikswo. (1989), Table 1
#sigma_l = 0.15			# mS / mm
sigma_t = 0.02			# mS / mm
beta = 200.0			# mm^{-1}
C_m = 0.2				# mu F / mm^2

def forward(sigma_l):
	# Define time
	time = Constant(0.0)

	# We use the Grandi cell model
	cell_model = Grandi_pasqualini_bers_2010()

	# Define the external stimulus
	p = Expression('10', degree=1, domain=mesh)
	stimulus = Expression('(x[0] > 2.0 and x[0] < 3.0 and x[1] > 2.0 \
	and x[1] < 3.0 and t < 3.0 ? p : 0)', p=p, t=time, degree=1)

	# Scale conducitivites by 1/(C_m * chi)
	M_l = sigma_l/(C_m*beta) # mm^2 / ms
	M_t = sigma_t/(C_m*beta) # mm^2 / ms

	# Define the conductivity (tensor)
	M = as_tensor(((M_l, 0.0), (0.0, M_t)))

	# Collect this information into the CardiacModel class
	cardiac_model = CardiacModel(mesh, time, M, 'none', cell_model, stimulus)

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
	dt = 1.0
	T = 10.0
	interval = (0.0, T)

	# Solve forward problem
	for (timestep, fields) in solver.solve(interval, dt):
	    print "(t_0, t_1) = ", timestep

	    # Extract the components of the field (vs_ at previous timestep,
	    (vs_, vs, vur) = fields

	return vur

# # Load synthetic data
# Q = FunctionSpace(mesh, "CG", 1)
# h = Function(Q)
# g = Function(Q)
# # j = Function(Q)
# hdf_v = HDF5File(mesh.mpi_comm(), "Results/Basic_test_case_monodomain_synthetic_observations_v.h5", "r")
# hdf_v.read(h, "v")
# # # hdf_c = HDF5File(mesh.mpi_comm(), "Results/Basic_test_case_monodomain_synthetic_observations_c.h5", "r")
# # attr_v = hdf_v.attributes("v")
# # # attr_c = hdf_c.attributes("c")
# # N = attr_v['count']
# # for timestep in range(N):
# #      dataset_v = "v/vector_%d"%timestep
# # #     dataset_c = "c/vector_%d"%timestep
# #      print 'Retrieving time step t=', timestep
# #      hdf_v.read(h, dataset_v)
# #      g.assign(h)
# # #     hdf_c.read(h, dataset_c)
# # #     j.assign(h)
# del hdf_v#, hdf_c

if __name__ == "__main__":
    sigma_l = Constant(0.15)     	  # initial guess.
    v = forward(sigma_l)              # solve the forward problem once. 

    #success = replay_dolfin(tol=0.0, stop=True)
    # adj_html("forward.html", "forward")
    # adj_html("adjoint.html", "adjoint")
    
    # Define functional of interest
    J = Functional(inner(v, v)*dx*dt[FINISH_TIME])

    # Indicate the control parameter of interest
    dJdnu = compute_gradient(J, Control(sigma_l))
    Jnu = assemble(inner(v, v)*dx) # current value

    parameters["adjoint"]["test_derivative"] = True
    parameters["adjoint"]["stop_annotating"] = True # stop registering equations

    def Jhat(sigma_l): # the functional as a pure function of M
         v = forward(sigma_l)
         return assemble(inner(v, v)*dx)

    conv_rate = taylor_test(Jhat, Control(sigma_l), Jnu, dJdnu, seed=(0.01))
    #conv_rate = taylor_test(Jhat, Control(M), Jnu, dJdnu)


# Define the reduced functional and solve the optimisation problem:
#rf = ReducedFunctional(J, M)
#opt_ctrls = minimize(rf, options={"maxiter": 50})

print "Success!"