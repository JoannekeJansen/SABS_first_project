# Solver for the inverse problem of Basic_test_case_monodomain.py
# ===========================================================
#

# Import the cbcbeat module
from cbcbeat import *
from mshr import *
import numpy as np

# Turn on FFC/FEniCS optimizations
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["cpp_optimize"] = True
flags = ["-O3", "-ffast-math", "-march=native"]
parameters["form_compiler"]["cpp_optimize_flags"] = " ".join(flags)
parameters["form_compiler"]["quadrature_degree"] = 3

# Define a [0,5]^2 domain
from mshr import *
domain = Rectangle(Point(0.0, 0.0), Point(5.0, 5.0))
mesh = generate_mesh(domain, 100)

# Averaged nominal conductivities, surface to volume ratio and membrane
# capacitance, as found in Sepulveda, Roth, & Wikswo. (1989), Table 1
sigma_l = 0.15         # mS / mm
sigma_t = 0.02          # mS / mm
beta = 200.0            # mm^{-1}
C_m = 0.2               # mu F / mm^2

#def forward(sigma_t):
def forward(a):
    # Define time
    time = Constant(0.0)

    # We use the Grandi cell model
    # cell_model = Grandi_pasqualini_bers_2010()
    cell_model = FitzHughNagumoManual() # For fast tests
    cell_model.set_parameters(a=a)

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
    vs_.assign(cell_model.initial_conditions(), solver.VS)
    vs.assign(cell_model.initial_conditions())

    # Time stepping parameters
    h = 0.5  # Time step size
    T = 5.0 # Final time
    interval = (0.0, T)

    # Solve forward problem
    for (timestep, fields) in solver.solve(interval, h):
        print "(t_0, t_1) = ", timestep
        # Extract the components of the field (vs_ at previous timestep,
        # current vs, current vur)
        (vs_, vs, vur) = fields

    return vs_, vs, vur

if __name__ == "__main__":
    a= Constant(0.15)                # initial guess, a=0.13 in the test problem.
    #b=Constant(0.013)                # initial guess, b=0.013 in the test problem.

    ctrl1 = a
    (vs_, vs, vur) = forward(a)   # solve the forward problem once. 

    # Load recorded data and define functional of interest
    times = np.loadtxt("Results/recorded_times_FitzH.txt")
    hdf_vs = HDF5File(mesh.mpi_comm(), "Results/Basic_test_case_monodomain_synthetic_observations_FitzH.h5", "r")
    attr_vs = hdf_vs.attributes("vs")
    N = attr_vs['count']
    Q = vs.function_space()
    vs_obs = {}
    I = 0
    for i in range(N):
      vs_obs[i] = Function(Q, annotate=False)
      dataset_vs = "vs/vector_%d"%i
      hdf_vs.read(vs_obs[i], dataset_vs)
      I = I + inner(split(vs)[0] - split(vs_obs[i])[0], split(vs)[0] - split(vs_obs[i])[0])*dx*dt[times[i]]
      #I = I + inner(split(vs)[38] - split(vs_obs[i])[38], split(vs)[38] - split(vs_obs[i])[38])*dx*dt[times[i]]
    del hdf_vs
    J=Functional(I/N)

    # Define the reduced functional and solve the optimisation problem:
    rf = ReducedFunctional(J, [Control(ctrl1)])
    assert rf.taylor_test(ctrl1, seed=1e-2) > 1.5
    #opt_ctrls = minimize(rf, tol=1e-02, options={"maxiter": 10})
    opt_ctrls = minimize(rf)
    print("ctrl1 = %f" %float(opt_ctrls))

# Taylor remainder without gradient information: [0.06534623028985469, 0.029487833810073666, 0.013959639519622466, 0.006783277021912969, 0.0033526482204846164]
# Convergence orders for Taylor remainder without gradient information (should all be 1): [1.147984166612108, 1.0788581579945282, 1.0412073699745041, 1.0166812971174115]
# Taylor remainder with gradient information: [0.010559247516349297, 0.002094342423320971, 0.00026289382624611855, 6.509582477520499e-05, 7.153820285947044e-05]
# Convergence orders for Taylor remainder with gradient information (should all be 2): [2.333937781424893, 2.9939451735844, 2.0138433448634547, -0.13614886388181136]

    

  