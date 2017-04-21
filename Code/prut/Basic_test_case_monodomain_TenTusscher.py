# A basic test case, based on demo_monodomain.py
# ===========================================================
#
# This test case solves the monodomain equations on a [0,10]^2 square,
# with the Ten Tusscher Panfilov cell model, using the splittingsolver with the default 
# parameter values and default initial conditions. 

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

# Define the conductivity (tensors)
M_i = 2.0
M_e = 1.0

# We use the Tentusscher Panfilov cell model
cell_model = Tentusscher_panfilov_2006_epi_cell()

# Define the external stimulus
p = Expression('500', degree=1, domain=mesh)
stimulus = Expression('(x[0] > 4.5 and x[0] < 5.5 and x[1] > 4.5 \
	and x[1] < 5.5 and t > 5.0 and t < 15.0 ? p : 0)', p=p, t=time, degree=1)

# Collect this information into the CardiacModel class
cardiac_model = CardiacModel(mesh, time, M_i, M_e, cell_model, stimulus)

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
vs_.assign(cell_model.initial_conditions())

# Time stepping parameters
dt = 0.01
T = 25.0
interval = (0.0, T)

vtkfile = File('Results/Basic_test_case_monodomain_TenTusscher.pvd')

# Solve
for (timestep, fields) in solver.solve(interval, dt):
    print "(t_0, t_1) = ", timestep

    # Extract the components of the field (vs_ at previous timestep,
    # current vs, current vur)
    (vs_, vs, vur) = fields
    vtkfile << vur

print "Success!"
