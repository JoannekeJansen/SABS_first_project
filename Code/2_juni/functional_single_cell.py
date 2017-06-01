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
parser.add_argument("--gna", type=int, help="The added percentage of noise to g_na")
parser.add_argument("--gk1", type=int, help="The added percentage of noise to g_k1")
parser.add_argument("--gkr", type=int, help="The added percentage of noise to g_kr")
parser.add_argument("--gcal", type=int, help="The added percentage of noise to g_cal")
parser.add_argument("--v", type=int, help="The added percentage of noise to the observed v")
parser.add_argument("--cai", type=int, help="The added percentage of noise to the observed cai")
parser.add_argument("--control_var", type=str, help="The control variable, options are gna, gcal, gk1 or gkr")
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
   observed_v ='observed_v_{0}_{1}_{2}_{3}_single_cell_{4}_{5}_{6}_{7}_{8}_{9}.txt'.format(gnaf,gk1f,gkrf,gcalf, \
     noise_percentage_gna, noise_percentage_gk1, noise_percentage_gkr, noise_percentage_gcal, noise_percentage_v, noise_percentage_cai)
   observed_cai ='observed_cai_{0}_{1}_{2}_{3}_single_cell_{4}_{5}_{6}_{7}_{8}_{9}.txt'.format(gnaf,gk1f,gkrf,gcalf, \
     noise_percentage_gna, noise_percentage_gk1, noise_percentage_gkr, noise_percentage_gcal, noise_percentage_v, noise_percentage_cai)
else:
  observed_times = "observed_times_{0}_{1}_{2}_{3}_single_cell.txt".format(gnaf,gk1f,gkrf,gcalf)
  observed_v = "observed_v_{0}_{1}_{2}_{3}_single_cell.txt".format(gnaf,gk1f,gkrf,gcalf)
  observed_cai = "observed_cai_{0}_{1}_{2}_{3}_single_cell.txt".format(gnaf,gk1f,gkrf,gcalf)

# For computing faster
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["cpp_optimize"] = True
flags = "-O3 -ffast-math -march=native"
parameters["form_compiler"]["cpp_optimize_flags"] = flags

# Define the external stimulus
class Stimulus(Expression):
    "Stimulus as in Paci2013."
    def __init__(self, **kwargs):
        self.t = kwargs["time"]
    def eval(self, value, x):
        if (float(self.t)-float(int(float(self.t)))) <= 0.005:
            v_amp = 5.5/0.987109 # A/F
            value[0] = v_amp
        else:
            value[0] = 0.0

# Load initial conditions
mesh =UnitIntervalMesh(1)
Q = VectorFunctionSpace(mesh, "CG", 1, 18)
vsinit = Function(Q)
hdf_vs = HDF5File(mesh.mpi_comm(), "initial_conditions_single_cell.h5", "r")
hdf_vs.read(vsinit, "vs/vector_0")
del hdf_vs

# Load recorded data
times = np.loadtxt(observed_times)
v_obs = np.loadtxt(observed_v)
cai_obs = np.loadtxt(observed_cai)

def forward(factor_gna, factor_gcal, factor_gk1, factor_gkr):
    model = Paci2013_ventricular_adjusted()

    # Set initial conditions
    model.set_initial_conditions(V=vsinit((0.0))[0], m=vsinit((0.0))[1], h=vsinit((0.0))[2], \
      j=vsinit((0.0))[3], d=vsinit((0.0))[4], f1=vsinit((0.0))[5], f2=vsinit((0.0))[6], \
      fCa=vsinit((0.0))[7], Xr1=vsinit((0.0))[8], Xr2=vsinit((0.0))[9], Xs=vsinit((0.0))[10], \
      Xf=vsinit((0.0))[11], q=vsinit((0.0))[12], r=vsinit((0.0))[13], Nai=vsinit((0.0))[14], \
      g=vsinit((0.0))[15], Cai=vsinit((0.0))[16], Ca_SR=vsinit((0.0))[17])

    # Set new parameters
    model.set_parameters(g_Na_factor=factor_gna)
    model.set_parameters(g_CaL_factor=factor_gcal)
    model.set_parameters(g_K1_factor=factor_gk1)
    model.set_parameters(g_Kr_factor=factor_gkr)        
    
    time = Constant(0.0)
    model.stimulus = Stimulus(time=time, degree=0)

    #model.set_initial_conditions(V=3.6153e-04, Cai=2.0586e-05)

    # Initialize solver
    params = SingleCellSolver.default_parameters()
    params["scheme"] = "GRL1"
    solver = SingleCellSolver(model, time, params)

    # Assign initial conditions
    (vs_, vs) = solver.solution_fields()
    vs_.assign(model.initial_conditions())

    # Solve and extract values
    h = 0.0001
    T = 1.0
    interval = (0.0, T)

    I = 0.0
    I1= 0.0
    I2=0.0
    i=0

    # Solve forward problem
    for (timestep, fields) in solver.solve(interval, h):
        # Extract the components of the field (vs_ at previous timestep,
        # current vs, current vur)
        if (abs(round(1000*timestep[1])-1000*timestep[1]) < 100*h) and (i < 1000):
            #print vs.vector().array()[0], v_obs[i]
            I = I + inner(vs.vector().array()[0] - v_obs[i], vs.vector().array()[0] - v_obs[i]) \
              / inner(v_obs[i], v_obs[i])
            I1 = I1 + inner(vs.vector().array()[0] - v_obs[i], vs.vector().array()[0]- v_obs[i]) \
              / inner(v_obs[i], v_obs[i])
            I = I +  inner(vs.vector().array()[16] - cai_obs[i], vs.vector().array()[16] - cai_obs[i]) \
              / inner(cai_obs[i], cai_obs[i])
            I2 = I2 + inner(vs.vector().array()[16] - cai_obs[i], vs.vector().array()[16] - cai_obs[i]) \
              / inner(cai_obs[i], cai_obs[i])
            i=i+1

    J = I/i
    J1 = I1/i
    J2 = I2/i

    return J, J1, J2

if __name__ == "__main__":

    controls = {'gna': 1.0, 'gcal': 1.0, 'gk1': 1.0, 'gkr': 1.0} # initial guess

    # Compute value of J for different values of the control_var and save to file
    control_values=np.arange(0.75,1.26,0.01)
    M = np.size(control_values)
    J_values = np.zeros(M)
    J1_values = np.zeros(M)
    J2_values = np.zeros(M)
    for i in range(np.size(control_values)):
      controls[args.control_var]=control_values[i]
      (J_values[i], J1_values[i], J2_values[i]) = forward(controls['gna'], controls['gcal'], controls['gk1'], controls['gkr']) 
      print "Control values and J value", control_values[i], J_values[i]
      np.savetxt('J_values_{0}_{1}_{2}_{3}_single_cell_{4}.txt'.format(gnaf,gk1f,gkrf,gcalf,args.control_var), np.column_stack((J_values, J1_values, J2_values, control_values)))

    del hdf_v, hdf_cai
        
    print "Success!"