# Create observations, to be used in inverse_12mm_strip.py
#
# Joanneke E Jansen, May 2017

# Import the cbcbeat module
from cbcbeat import *
import numpy as np
import bisect
import argparse

# Set parameters
parser = argparse.ArgumentParser()
parser.add_argument("-gnaf", type=float, help="The factor of g_na")
parser.add_argument("-gk1f", type=float, help="The factor of g_k1")
parser.add_argument("-gkrf", type=float, help="The factor of g_kr")
parser.add_argument("-gcalf", type=float, help="The factor of g_cal")
args = parser.parse_args()
factor_gna = args.gnaf
factor_gk1 = args.gk1f
factor_gkr = args.gkrf
factor_gcal = args.gcalf

# Disable adjointing
parameters["adjoint"]["stop_annotating"] = True

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

def main(factor_gna, factor_gcal, factor_gk1, factor_gkr):
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

    # Initialize solver
    params = SingleCellSolver.default_parameters()
    params["scheme"] = "GRL1"
    solver = SingleCellSolver(model, time, params)

    # Assign initial conditions
    (vs_, vs) = solver.solution_fields()
    vs_.assign(model.initial_conditions())

    # Solve and extract values
    dt = 0.0001
    T = 1.0
    interval = (0.0, T)

    # Create files to save observed v and [CA]i
    v_obs = np.zeros(1000)
    cai_obs = np.zeros(1000)
    times = np.zeros(1000)
    i=0

    for ((t0, t1), fields) in solver.solve(interval, dt):
      print vs.vector().array()[0]
      #Record vs at each ms
      if (abs(round(1000*t1)-1000*t1)) < (100*dt) and (i < 1000):
        times[i] = t1
        v_obs[i] = vs.vector().array()[0]
        cai_obs[i] = vs.vector().array()[16]
        i=i+1
        #print i

    # Save
    np.savetxt("observed_times_{0}_{1}_{2}_{3}_single_cell.txt".format(factor_gna,factor_gk1,factor_gkr,factor_gcal), times)
    np.savetxt("observed_v_{0}_{1}_{2}_{3}_single_cell.txt".format(factor_gna,factor_gk1,factor_gkr,factor_gcal), v_obs)
    np.savetxt("observed_cai_{0}_{1}_{2}_{3}_single_cell.txt".format(factor_gna,factor_gk1,factor_gkr,factor_gcal), cai_obs)

# Record
main(factor_gna, factor_gcal, factor_gk1, factor_gkr)

print "Success!"
