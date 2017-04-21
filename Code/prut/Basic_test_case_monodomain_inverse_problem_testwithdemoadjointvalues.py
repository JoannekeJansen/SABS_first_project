# Solver for the inverse problem of Basic_test_case_monodomain.py
# ===========================================================
#

# Import the cbcbeat module
from cbcbeat import *
import numpy.random

# Turn on FFC/FEniCS optimizations
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["cpp_optimize"] = True
flags = ["-O3", "-ffast-math", "-march=native"]
parameters["form_compiler"]["cpp_optimize_flags"] = " ".join(flags)
parameters["form_compiler"]["quadrature_degree"] = 3

# Define the computational domain
mesh = UnitSquareMesh(100, 100)

def forward(M_i):
	time = Constant(0.0)
	# Pick a cell model (see supported_cell_models for tested ones)
	#cell_model = Tentusscher_panfilov_2006_epi_cell()
	cell_model = FitzHughNagumoManual()

	# Define some external stimulus
	#stimulus = Expression("(x[0] > 0.9 && t <= 1.0) ? 30.0 : 0.0",
     #                 t=time, degree=0)
	p = Expression('25', degree=1, domain=mesh)
	stimulus = Expression('(x[0] > 4.5 and x[0] < 5.5 and x[1] > 4.5 \
		and x[1] < 5.5 and t > 0 and t < 100 ? p : 0)', p=p, t=time, degree=1)

	# Collect this information into the CardiacModel class
	cardiac_model = CardiacModel(mesh, time, M_i, None, cell_model, stimulus)

	# Customize and create a splitting solver
	ps = SplittingSolver.default_parameters()
	ps["theta"] = 0.5                        # Second order splitting scheme
	ps["pde_solver"] = "monodomain"          # Use Monodomain model for the PDEs
	ps["CardiacODESolver"]["scheme"] = "RL1" # 1st order Rush-Larsen for the ODEs
	ps["MonodomainSolver"]["linear_solver_type"] = "iterative"
	ps["MonodomainSolver"]["algorithm"] = "cg"
	ps["MonodomainSolver"]["preconditioner"] = "petsc_amg"

	solver = SplittingSolver(cardiac_model, params=ps)

	# Extract the solution fields and set the initial conditions
	(vs_, vs, vur) = solver.solution_fields()
	vs_.assign(cell_model.initial_conditions(), solver.VS)

	# Time stepping parameters
	k = 0.01
	T = 0.1
	interval = (0.0, T)

	# Solve forward problem
	for (timestep, fields) in solver.solve(interval, k):
	    print "(t_0, t_1) = (%g, %g)" % timestep
	    (vs_, vs, vur) = fields

	return vs

if __name__ == "__main__":
    M_i = Constant(2.0)	       	  # initial guess.
    vs = forward(M_i)             # solve the forward problem once.

    # Define functional of interest
    j = inner(vs, vs)*dx*dt[FINISH_TIME]
    J = Functional(j)

    # Indicate the control parameter of interest
    m = Control(M_i) 

    # Indicate the control parameter of interest
    dJdnu = compute_gradient(J, Control(M_i))
    Jnu = assemble(inner(vs, vs)*dx) # current value

    parameters["adjoint"]["stop_annotating"] = True # stop registering equations

    def Jhat(M_i): # the functional as a pure function of Gna
         vs = forward(M_i)
         return assemble(inner(vs, vs)*dx)

    conv_rate = taylor_test(Jhat, Control(M_i), Jnu, dJdnu, seed=(1))


print "Success!"