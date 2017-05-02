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

# Define a [0,5]^2 domain
mesh = RectangleMesh(Point(0.0, 0.0), Point(0.25, 0.25), 5, 5)

# Define averaged nominal conductivities, surface to volume ratio and membrane
# capacitance, as found in Sepulveda, Roth, & Wikswo. (1989), Table 1
sigma_l = 0.15         # mS / mm
sigma_t = 0.02          # mS / mm
beta = 200.0            # mm^{-1}
C_m = 0.2               # mu F / mm^2

def forward(g_CaL_factor, g_K1_factor, g_Kr_factor):
    # Define time
    time = Constant(0.0)

    # Use the Grandi cell model
    cell_model = Grandi_pasqualini_bers_2010_adjusted()
    cell_model.set_parameters(g_CaL_factor=g_CaL_factor)
    cell_model.set_parameters(g_K1_factor=g_K1_factor)
    cell_model.set_parameters(g_Kr_factor=g_Kr_factor)
    #cell_model = FitzHughNagumoManual() # For fast tests

    # Define the external stimulus
    stimulus = Expression('(x[0] < 0.1 && t < 10.0 ? p : 0)', p=10.0, t=time, degree=1)

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
    T = 10.0   # Final time
    interval = (0.0, T)

    # Solve forward problem
    for (timestep, fields) in solver.solve(interval, h):
        # Extract the components of the field (vs_ at previous timestep,
        # current vs, current vur)
        (vs_, vs, vur) = fields
        #print "(t_0, t_1) = ", timestep
        #print "v(2.5,2.5) = ", vs((2.5,2.5))[0]

    return vs, vur

if __name__ == "__main__":
    Optimize = True
    Evaluate_J = False
    g_CaL_factor = Constant(0.50)          # initial guess
    g_K1_factor = Constant(0.50)          # initial guess
    g_Kr_factor = Constant(0.50)          # initial guess
    ctrl1 = g_CaL_factor
    ctrl2 = g_K1_factor
    ctrl3 = g_Kr_factor
    (vs, vur) = forward(ctrl1, ctrl2, ctrl3)               # solves the forward problem once. 

    # Load recorded data and define functional of interest
    times = np.loadtxt("recorded_times_adjusted_Ca100_mini.txt")
    hdf_v = HDF5File(mesh.mpi_comm(), "Basic_test_case_monodomain_synthetic_observations_v_adjusted_mini.h5", "r")
    hdf_cai = HDF5File(mesh.mpi_comm(), "Basic_test_case_monodomain_synthetic_observations_cai_adjusted_mini.h5", "r")
    attr_v = hdf_v.attributes("v")
    attr_cai = hdf_cai.attributes("cai")
    N = attr_v['count']
    print "N=", N
    Q = vur.function_space()
    I = 0
    for i in range(N-490):
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
 
    # def hessian_cb(j, m, mdot, h):
    #     print "j = %f, dj = %f, m = %f,  mdot = %f,  h = %f." % (j, dj, float(m), float(mdot), float(h))

    rf = ReducedFunctional(J, [Control(ctrl1), Control(ctrl2)] , eval_cb_post = eval_cb_post, derivative_cb_post = derivative_cb_post)
    #rf = ReducedFunctional(J, [Control(ctrl3)] , eval_cb_post = eval_cb_post, derivative_cb_post = derivative_cb_post)

       # Solve the optimisation problem
    if Optimize == True:
        # assert rf.taylor_test(ctrl1, seed=1e-2) > 1.5
        # rf.taylor_test(ctrl1, seed=1e-1)
        opt_ctrls = minimize(rf, tol=1e-03, options={"maxiter": 10, "gtol": 1e-15})
        #print("ctrl1 = %f" %float(opt_ctrls))
    
    # Compute value of J for different values of ctrl1 and save to file
    if Evaluate_J == True:
        M = 11 # Number of functional evaluations
        J_values = np.zeros(M)
        orig_val = 1.0    # the value of the control variable that was used to generate the recorded data
        control_values=np.array([0.5*orig_val,0.6*orig_val,0.7*orig_val,0.8*orig_val,0.9*orig_val, 1.0*orig_val,1.1*orig_val,1.2*orig_val,1.3*orig_val,1.4*orig_val,1.5*orig_val])
        for i in range(np.size(control_values)):
            J_values[i]=rf(Constant(control_values[i]))
            print "Control values and J value", control_values[i], J_values[i]
        np.savetxt("J_values.txt", np.column_stack((J_values, control_values)))

    #timer.stop()
    #list_timings(TimingClear_keep, [TimingType_wall])
    print "Success!"

# g_CaL
# 0.5 1.93527937638e-05
# 0.6 1.18803170559e-05
# 0.7 6.71780805334e-06
# 0.8 3.55040157179e-06
# 0.9 8.14333575251e-07
# 1.0 5.00047004073e-07
# 1.1 1.33222176686e-06
# 1.2 2.59156339615e-06
# 1.3 5.61789511944e-06
# 1.4 9.96937546274e-06
# 1.5 1.46064834398e-05

# g_K1
# 5.370963280342857735e-03 5.999999999999999778e-01
# 4.828581202845943360e-03 6.999999999999999556e-01
# 3.432023925061003251e-03 8.000000000000000444e-01
# 7.070887115231877301e-04 9.000000000000000222e-01
# 1.962715577320342202e-05 1.000000000000000000e+00
# 7.352872828873287378e-04 1.100000000000000089e+00
# 2.563403554148015358e-03 1.199999999999999956e+00
# 6.825743546286359466e-03 1.300000000000000044e+00
# 1.367088387259096822e-02 1.399999999999999911e+00
# 1.935968460402388500e-02 1.500000000000000000e+00

# g_Kr
# 5.790092186946272854e-03 5.000000000000000000e-01
# 5.776799527707745553e-03 5.999999999999999778e-01
# 5.780594989351353456e-03 6.999999999999999556e-01
# 5.785016890513867971e-03 8.000000000000000444e-01
# 5.771505782650183068e-03 9.000000000000000222e-01
# 5.782347917867306647e-03 1.000000000000000000e+00
# 5.776781740422678828e-03 1.100000000000000089e+00
# 5.764588662044315380e-03 1.199999999999999956e+00
# 5.774912376675597932e-03 1.300000000000000044e+00
# 5.774468913767038213e-03 1.399999999999999911e+00
# 5.761168727772681200e-03 1.500000000000000000e+00


# j  =  1.88921402486e-05
# dj =  [-8.086349720218391e-05, -0.0003479982830444868]
# m  =  [0.5112754539579218, 0.9970791725958807]
# >>>>

# At iterate    2    f=  1.88921D-05    |proj g|=  3.47998D-04

#            * * *



