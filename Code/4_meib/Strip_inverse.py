# Solver for the inverse problem of Basic_test_case_monodomain.py
# ===========================================================
#
# Joanneke E Jansen, April 2017

# Import the cbcbeat module
from cbcbeat import *
import numpy as np

# Set log level
#set_log_level(WARNING)

#timer = Timer("xxx:Solve")

#adj_checkpointing(strategy='multistage', steps=int(100*2)+1,
#                  snaps_on_disk=2, snaps_in_ram=2, verbose=True)

# Turn on FFC/FEniCS optimizations
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["cpp_optimize"] = True
flags = ["-O3", "-ffast-math", "-march=native"]
parameters["form_compiler"]["cpp_optimize_flags"] = " ".join(flags)
parameters["form_compiler"]["quadrature_degree"] = 3

# Define a [0,5]x [0,0.01] domain
mesh = RectangleMesh(Point(0.0, 0.0), Point(5.0, 0.01), 100, 1)

# Define averaged nominal conductivities, surface to volume ratio and membrane
# capacitance, as found in Sepulveda, Roth, & Wikswo. (1989), Table 1
sigma_l = 0.15         # mS / mm
sigma_t = 0.02          # mS / mm
beta = 200.0            # mm^{-1}
C_m = 0.2               # mu F / mm^2

def forward(GNa):
    # Define time
    time = Constant(0.0)

    # Use the Grandi cell model
    cell_model = Grandi_pasqualini_bers_2010()
    cell_model.set_parameters(GNa=GNa)

    # Define the external stimulus
    stimulus = Expression('(x[0] < 0.5 && t < 10.0 ? p : 0)', p=10.0, t=time, degree=1)

    # Scale conducitivites by 1/(C_m * chi)
    M_l = sigma_l/(C_m*beta) # mm^2 / ms
    M_t = sigma_t/(C_m*beta) # mm^2 / ms

    # Define the conductivity (tensor)
    M = as_tensor(((M_l, 0.0), (0.0, M_t)))

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

    # Time stepping parameters
    h = 0.01 # Time step size
    T = 500.0   # Final time
    interval = (0.0, T)

    # Solve forward problem
    for (timestep, fields) in solver.solve(interval, h):
        # Extract the components of the field (vs_ at previous timestep,
        # current vs, current vur)
        (vs_, vs, vur) = fields

    return vs, vur

if __name__ == "__main__":
    Optimize = True
    Evaluate_J = False
    GNaf = Constant(0.50)          # initial guess
    ctrl1 = GNaf
    Gna = Constant(ctrl1*23)
    (vs, vur) = forward(Gna)    # solves the forward problem once. 

    # Load recorded data and define functional of interest
    times = np.loadtxt("recorded_times.txt")
    hdf_v = HDF5File(mesh.mpi_comm(), "Strip_synthetic_observations_v.h5", "r")
    hdf_cai = HDF5File(mesh.mpi_comm(), "Strip_synthetic_observations_cai.h5", "r")
    attr_v = hdf_v.attributes("v")
    attr_cai = hdf_cai.attributes("cai")
    N = attr_v['count']
    print "N=", N
    Q = vur.function_space()
    I = 0
    for i in range(N):
        print i
        v_obs = Function(Q)
        cai_obs = Function(Q)
        dataset_v = "v/vector_%d"%i
        dataset_cai = "cai/vector_%d"%i
        hdf_v.read(v_obs, dataset_v)
        hdf_cai.read(cai_obs, dataset_cai)
        I = I + (inner(split(vs)[0] - v_obs, split(vs)[0] - v_obs)*dx*dt[float(times[i])]) \
            /(assemble(inner(v_obs, v_obs)*dx))
        I = I + (inner(split(vs)[36] - cai_obs, split(vs)[36] - cai_obs)*dx*dt[float(times[i])]) \
            /(assemble(inner(cai_obs, cai_obs)*dx))
    del hdf_v, hdf_cai
    J = Functional(I/N)

    def hessian_cb(j, m, mdot, h):
        """
        output after each function evaluation
        """
        print "<<<< eval"
        print 'j  = ', j
        print 'm  = ', [float(va) for va in m]
        print 'mdot  = ', [float(va) for va in mdot]
        print 'h  = ', [float(va) for va in h]
        print ">>>>"

    # Define the reduced functional
    def derivative_cb_post(j, dj, m):
        """
        output after each derivative evaluation
        """
        print "<<<< grad"
        print 'j  = ', j
        print 'dj = ', [float(va) for va in dj]
        print 'm  = ', [float(va) for va in m]
        print ">>>>"

    def eval_cb_post(j, m):
        """
        output after each function evaluation
        """
        print "<<<< eval"
        print 'j  = ', j
        print 'm  = ', [float(va) for va in m]
        print ">>>>"
 
    rf = ReducedFunctional(J, [Control(ctrl1)] , eval_cb_post = eval_cb_post, derivative_cb_post = derivative_cb_post, hessian_cb = hessian_cb)

       # Solve the optimisation problem
    if Optimize == True:
        # assert rf.taylor_test(ctrl1, seed=1e-2) > 1.5
        # rf.taylor_test(ctrl1, seed=1e-1)
        opt_ctrls = minimize(rf, tol=1e-03, options={"maxiter": 10, "gtol": 1e-15})
    
    # Compute value of J for different values of ctrl1 and save to file
    if Evaluate_J == True:
        M = 19 # Number of functional evaluations
        J_values = np.zeros(M)
        orig_val = 1.0    # the value of the control variable that was used to generate the recorded data
        control_values=np.array([0.5*orig_val,0.6*orig_val,0.7*orig_val,0.8*orig_val,0.9*orig_val, 1.0*orig_val,1.1*orig_val,1.2*orig_val,1.3*orig_val,1.4*orig_val,1.5*orig_val, \
            0.4*orig_val,0.3*orig_val,0.2*orig_val,0.1*orig_val,1.6*orig_val,1.7*orig_val,1.8*orig_val,1.9*orig_val])
        for i in range(np.size(control_values)):
            J_values[i]=rf(Constant(control_values[i]))
            print "Control values and J value", control_values[i], J_values[i]
        np.savetxt("J_values.txt", np.column_stack((J_values, control_values)))

    #timer.stop()
    #list_timings(TimingClear_keep, [TimingType_wall])
    
    print "Success!"


