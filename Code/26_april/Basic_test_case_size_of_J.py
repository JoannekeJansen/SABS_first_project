# Solver for the inverse problem of Basic_test_case_monodomain.py
# ===========================================================
#
# Joanneke E Jansen, April 2017

# Import the cbcbeat module
from cbcbeat import *
import numpy as np

timer = Timer("xxx:Solve")

# Turn off adjoint functionality
parameters["adjoint"]["stop_annotating"] = True

# Turn on FFC/FEniCS optimizations
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["cpp_optimize"] = True
flags = ["-O3", "-ffast-math", "-march=native"]
parameters["form_compiler"]["cpp_optimize_flags"] = " ".join(flags)
parameters["form_compiler"]["quadrature_degree"] = 3

# Define a [0,5]^2 domain
mesh = RectangleMesh(Point(0.0, 0.0), Point(5.0, 5.0), 50, 50)

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
    #cell_model = FitzHughNagumoManual() # For fast tests

    # Define the external stimulus
    stimulus = Expression('(x[0] > 2.25 && x[0] < 2.75 && x[1] > 2.25 \
    && x[1] < 2.75 && t < 3.0 ? p : 0)', p=10.0, t=time, degree=1)

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
    h = 0.05 # Time step size
    T = 2.0   # Final time
    interval = (0.0, T)

    # Create HDF5 file to save vs
    hdf_vs_model = HDF5File(mesh.mpi_comm(), "Results/Basic_test_case_monodomain_model_observations_vs.h5", "w")
    i=0

    # Solve forward problem
    for (timestep, fields) in solver.solve(interval, h):
        # Extract the components of the field (vs_ at previous timestep,
        # current vs, current vur)
        (vs_, vs, vur) = fields
        #print "(t_0, t_1) = ", timestep
        #Record vs at each time step or at each ms
        if (abs(round(timestep[1])-timestep[1])) < (0.1*h):
          hdf_vs_model.write(vs,"vs",i)
          i=i+1
    del hdf_vs_model
    return vs

def size_of_cost_function(ctrl1):
    vs = forward(ctrl1)
    times = np.loadtxt("Results/recorded_times.txt")
    hdf_vs = HDF5File(mesh.mpi_comm(), "Results/Basic_test_case_monodomain_synthetic_observations_vs.h5", "r")
    hdf_vs_model = HDF5File(mesh.mpi_comm(), "Results/Basic_test_case_monodomain_model_observations_vs.h5", "r")
    attr_vs = hdf_vs.attributes("vs")
    N = attr_vs['count']
    print "N=", N
    Q = vs.function_space()
    vs_obs = {}
    I = 0
    for i in range(N):
        vs_obs = Function(Q, annotate=False)
        vs_model = Function(Q, annotate=False)
        dataset_vs = "vs/vector_%d"%i
        hdf_vs.read(vs_obs, dataset_vs)
        hdf_vs_model.read(vs_model, dataset_vs)
        I = I + (assemble(inner(split(vs_model)[0] - split(vs_obs)[0], split(vs_model)[0] - split(vs_obs)[0])*dx) \
            /assemble(inner(split(vs_obs)[0], split(vs_obs)[0])*dx))
        I = I + (assemble(inner(split(vs_model)[36] - split(vs_obs)[36], split(vs_model)[36] - split(vs_obs)[36])*dx) \
            /assemble(inner(split(vs_obs)[36], split(vs_obs)[36])*dx))
    del hdf_vs, hdf_vs_model
    return (I/N)

if __name__ == "__main__":
    #GNa= Constant(23)                # initial guess, GNa=23 in the test problem.
    M = 2 # Number of functional evaluations
    J_values = np.zeros(M)
    orig_val = 23    # the value of the control variable that was used to generate the recorded data
    control_values=np.linspace(0.5*orig_val, 2*orig_val, num=M)
    for i in range(np.size(control_values)):
        print "Control_value =", control_values[i]
        J_values[i]=size_of_cost_function(control_values[i])
        print "J_value = ", J_values[i]
    np.savetxt("Results/J_values.txt", np.column_stack((J_values, control_values)))

    print "Success!"
    timer.stop()
    list_timings(TimingClear_keep, [TimingType_wall])
