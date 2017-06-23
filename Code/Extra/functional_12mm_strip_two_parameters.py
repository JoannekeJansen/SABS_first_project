# Solver for the inverse problem of Basic_test_case_monodomain.py
# ===========================================================
#
# Joanneke E Jansen, April 2017

# Import the cbcbeat module
from cbcbeat import *
import numpy as np
import time as time1
import argparse

# Set log level
set_log_level(WARNING)

# Turn off adjoint functionality
parameters["adjoint"]["stop_annotating"] = True

# Retrieve parameters
parser = argparse.ArgumentParser()
parser.add_argument("-gnaf", type=float, help="The factor of g_na")
parser.add_argument("-gk1f", type=float, help="The factor of g_k1")
parser.add_argument("-gkrf", type=float, help="The factor of g_kr")
parser.add_argument("-gcalf", type=float, help="The factor of g_cal")
parser.add_argument("-gna", type=int, help="The added percentage of noise to g_na")
parser.add_argument("-gk1", type=int, help="The added percentage of noise to g_k1")
parser.add_argument("-gkr", type=int, help="The added percentage of noise to g_kr")
parser.add_argument("-gcal", type=int, help="The added percentage of noise to g_cal")
parser.add_argument("-v", type=int, help="The added percentage of noise to the observed v")
parser.add_argument("-cai", type=int, help="The added percentage of noise to the observed cai")
parser.add_argument("-control_var", type=str, help="The first of the two control variables, options are gna, gcal, gk1 or gkr")
parser.add_argument("-control_var2", type=str, help="The second of the two control variables, options are gna, gcal, gk1 or gkr")
args = parser.parse_args()

# Set parameter factors
global gnaf
global gk1f
global gkrf
global gcalf
gnaf = args.gnaf
gk1f = args.gk1f
gkrf = args.gkrf
gcalf = args.gcalf

# Set noise parameters
global noise_percentage_gna
global noise_percentage_gk1
global noise_percentage_gkr
global noise_percentage_gcal
global noise_percentage_v
global noise_percentage_cai

if args.gna:
  noise_percentage_gna = args.gna
else:
  noise_percentage_gna = 0
if args.gk1:
  noise_percentage_gk1 = args.gk1
else:
  noise_percentage_gk1 = 0
if args.gkr:
  noise_percentage_gkr = args.gkr
else:
  noise_percentage_gkr = 0
if args.gcal:
  noise_percentage_gcal = args.gcal
else:
  noise_percentage_gcal = 0
if args.v:
  noise_percentage_v = args.v
else:
  noise_percentage_v = 0
if args.cai:
  noise_percentage_cai = args.cai
else:
  noise_percentage_cai = 0
if args.gna or args.gk1 or args.gkr or args.gcal or args.v or args.cai:
   observed_times ='observed_times_{0}_{1}_{2}_{3}_12mm_strip_{4}_{5}_{6}_{7}_{8}_{9}.txt'.format(gnaf,gk1f,gkrf,gcalf, \
     noise_percentage_gna, noise_percentage_gk1, noise_percentage_gkr, noise_percentage_gcal, noise_percentage_v, noise_percentage_cai)
   observed_v ='observed_v_{0}_{1}_{2}_{3}_12mm_strip_{4}_{5}_{6}_{7}_{8}_{9}.h5'.format(gnaf,gk1f,gkrf,gcalf, \
     noise_percentage_gna, noise_percentage_gk1, noise_percentage_gkr, noise_percentage_gcal, noise_percentage_v, noise_percentage_cai)
   observed_cai ='observed_cai_{0}_{1}_{2}_{3}_12mm_strip_{4}_{5}_{6}_{7}_{8}_{9}.h5'.format(gnaf,gk1f,gkrf,gcalf, \
     noise_percentage_gna, noise_percentage_gk1, noise_percentage_gkr, noise_percentage_gcal, noise_percentage_v, noise_percentage_cai)
else:
  observed_times = "observed_times_{0}_{1}_{2}_{3}_12mm_strip.txt".format(gnaf,gk1f,gkrf,gcalf)
  observed_v = "observed_v_{0}_{1}_{2}_{3}_12mm_strip.h5".format(gnaf,gk1f,gkrf,gcalf)
  observed_cai = "observed_cai_{0}_{1}_{2}_{3}_12mm_strip.h5".format(gnaf,gk1f,gkrf,gcalf)

# Turn on FFC/FEniCS optimizations
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["cpp_optimize"] = True
flags = ["-O3", "-ffast-math", "-march=native"]
parameters["form_compiler"]["cpp_optimize_flags"] = " ".join(flags)
parameters["form_compiler"]["quadrature_degree"] = 3

# Define a [0,12]*[0,0.01] mm^2 rectangular domain
mesh = RectangleMesh(Point(0.0, 0.0), Point(12.00, 0.01), 200, 1)

# Load initial conditions
Q = VectorFunctionSpace(mesh, "CG", 1, 18)
vsinit = Function(Q)
hdf_vs = HDF5File(mesh.mpi_comm(), "initial_conditions_12mm_strip.h5", "r")
hdf_vs.read(vsinit, "vs/vector_0")
del hdf_vs

# Load recorded data
times = np.loadtxt(observed_times)
hdf_v = HDF5File(mesh.mpi_comm(), observed_v, "r")
hdf_cai = HDF5File(mesh.mpi_comm(), observed_cai, "r")
attr_v = hdf_v.attributes("v")
attr_cai = hdf_cai.attributes("cai")
N = attr_v['count']
Q = FunctionSpace(mesh, "CG", 1, 1)

def forward(factor_gna, factor_gcal, factor_gk1, factor_gkr):
    # Define time
    time = Constant(0.0)

    # Define the external stimulus
    stimulus = Expression('(x[0] < 0.5 && (t-floor(t)) < 0.005 ? v_amp : 0.0)', v_amp=5.5/0.987109, t=time, degree=1)

    # Define averaged nominal conductivities as found in Sepulveda, Roth, & Wikswo. (1989), Table 1
    sigma_t = 20.0                        # mu S / mm
    sigma_l = 200.0                       # mu S / mm

    # Take total cell capacitance and volume as in the Paci2013 model, where we assume that the 
    # total volume is the sum of the intracellular and SR volumes.
    C_m = 9.87109e-5                      # mu F
    Volume = 9.3837e-6                    # mm^3

    # Scale conducitivites by Volume/C_m
    M_t = sigma_t*Volume/C_m              # mm^2 / s
    M_l = sigma_l*Volume/C_m              # mm^2 / s

    # Define the conductivity (tensor) 
    M = as_tensor(((M_t, 0.0), (0.0, M_l)))

    # Use the Paci2013 ventricular cell model
    model = Paci2013_ventricular_adjusted()

    # Set initial conditions
    model.set_initial_conditions(V=vsinit.split(deepcopy=True)[0], m=vsinit.split(deepcopy=True)[1], h=vsinit.split(deepcopy=True)[2], \
      j=vsinit.split(deepcopy=True)[3], d=vsinit.split(deepcopy=True)[4], f1=vsinit.split(deepcopy=True)[5], f2=vsinit.split(deepcopy=True)[6], \
      fCa=vsinit.split(deepcopy=True)[7], Xr1=vsinit.split(deepcopy=True)[8], Xr2=vsinit.split(deepcopy=True)[9], Xs=vsinit.split(deepcopy=True)[10], \
      Xf=vsinit.split(deepcopy=True)[11], q=vsinit.split(deepcopy=True)[12], r=vsinit.split(deepcopy=True)[13], Nai=vsinit.split(deepcopy=True)[14], \
      g=vsinit.split(deepcopy=True)[15], Cai=vsinit.split(deepcopy=True)[16], Ca_SR=vsinit.split(deepcopy=True)[17])
    
    # Set new parameters
    model.set_parameters(g_Na_factor=factor_gna)
    model.set_parameters(g_CaL_factor=factor_gcal)
    model.set_parameters(g_K1_factor=factor_gk1)
    model.set_parameters(g_Kr_factor=factor_gkr)

    # Collect this information into the CardiacModel class
    cardiac_model = CardiacModel(mesh, time, M, 'none', model, stimulus)

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
    vs_.assign(model.initial_conditions(), solver.VS)
    vs.assign(model.initial_conditions())

    # Set time stepping parameters
    h = 0.0001 # Time step size
    T = 1.00    # Final time
    interval = (0.0, T)

    # 
    I = 0.0
    I1= 0.0
    I2=0.0
    i=0

    # Solve forward problem
    for (timestep, fields) in solver.solve(interval, h):
        # Extract the components of the field (vs_ at previous timestep,
        # current vs, current vur)
        (vs_, vs, vur) = fields
        if (abs(round(1000*timestep[1])-1000*timestep[1]) < 100*h) and (i < 1000):
            v_obs = Function(Q)
            cai_obs = Function(Q)
            dataset_v = "v/vector_%d"%i
            dataset_cai = "cai/vector_%d"%i
            hdf_v.read(v_obs, dataset_v)
            hdf_cai.read(cai_obs, dataset_cai)
            I = I + (assemble(inner(split(vs)[0] - v_obs, split(vs)[0] - v_obs)*dx)) \
              /(assemble(inner(v_obs, v_obs)*dx))
            I1 = I1 + (assemble(inner(split(vs)[0] - v_obs, split(vs)[0] - v_obs)*dx)) \
              /(assemble(inner(v_obs, v_obs)*dx))
            I = I + (assemble(inner(split(vs)[16] - cai_obs, split(vs)[16] - cai_obs)*dx)) \
              /(assemble(inner(cai_obs, cai_obs)*dx))
            I2 = I2 + (assemble(inner(split(vs)[16] - cai_obs, split(vs)[16] - cai_obs)*dx)) \
              /(assemble(inner(cai_obs, cai_obs)*dx))
            i=i+1

    J = I/i
    J1 = I1/i
    J2 = I2/i

    return J, J1, J2

if __name__ == "__main__":

    controls = {'gna': 1.0, 'gcal': 1.0, 'gk1': 1.0, 'gkr': 1.0} # initial guess

    # Compute value of J for different values of the control_var and save to file
    control_values=np.arange(0.75,1.26,0.01)
    control_values2=np.arange(0.75,1.26,0.01)
    M = np.size(control_values)
    M2 = np.size(control_values2)
    J_values = np.zeros(M*M2)
    J1_values = np.zeros(M*M2)
    J2_values = np.zeros(M*M2)
    C = np.zeros(M*M2)
    C2 = np.zeros(M*M2)
    k=0
    for i in range(np.size(control_values)):
      for j in range(np.size(control_values2)):
        controls[args.control_var]=control_values[i]
        controls[args.control_var2]=control_values2[j]
        C[k]=controls[args.control_var]
        C2[k]=controls[args.control_var2]
        (J_values[k], J1_values[k], J2_values[k]) = forward(controls['gna'], controls['gcal'], controls['gk1'], controls['gkr']) 
        print "Control values and J value", C[k], C2[k], J_values[k]
        np.savetxt('J_values_{0}_{1}_{2}_{3}_12mm_strip_{4}_{5}.txt'.format(gnaf,gk1f,gkrf,gcalf,args.control_var,args.control_var2), np.column_stack((J_values, J1_values, J2_values, C, C2)))
        k=k+1

    del hdf_v, hdf_cai
        
    print "Success!"