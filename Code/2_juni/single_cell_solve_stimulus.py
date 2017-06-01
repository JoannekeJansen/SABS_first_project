# A single cell solve
# ===========================================================
#
# This test case solves the Paci2013 model (either atrial or ventricular) in a single cell, with a 1Hz stimulus.
#
# Joanneke E Jansen, May 2017

from cbcbeat import *
import numpy as np

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

def main(factor_gna, factor_gcal, factor_gk1, factor_gkr, factor_gto, factor_gf, factor_gks):
    #model = Paci2013_atrial()
    model = Paci2013_ventricular_adjusted()

    # Set initial conditions
    model.set_initial_conditions(V=vsinit((0.0))[0], m=vsinit((0.0))[1], h=vsinit((0.0))[2], \
    j=vsinit((0.0))[3], d=vsinit((0.0))[4], f1=vsinit((0.0))[5], f2=vsinit((0.0))[6], \
    fCa=vsinit((0.0))[7], Xr1=vsinit((0.0))[8], Xr2=vsinit((0.0))[9], Xs=vsinit((0.0))[10], \
    Xf=vsinit((0.0))[11], q=vsinit((0.0))[12], r=vsinit((0.0))[13], Nai=vsinit((0.0))[14], \
    g=vsinit((0.0))[15], Cai=vsinit((0.0))[16], Ca_SR=vsinit((0.0))[17])

    # Extract default parameters
    params = Paci2013_ventricular.default_parameters()
    model.set_parameters(g_to=factor_gto*params["g_to"])
    model.set_parameters(g_f=factor_gf*params["g_f"])
    model.set_parameters(g_Ks=factor_gks*params["g_Ks"])

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
    dt = 0.0001
    T = 1.0
    interval = (0.0, T)

    # Create arrays to save time, v and cai
    times = np.zeros(int(T/dt))
    values_v = np.zeros(int(T/dt))
    values_c = np.zeros(int(T/dt))
    i=0

    for ((t0, t1), fields) in solver.solve(interval, dt):
        if True:
        #if (abs(round(t1)-t1) < 0.1*dt) & (i < int(500)):
            #print t1
            print vs.vector().array()[0]
            #print vs.vector().array()[16]
            times[i] = t1
            values_v[i] = vs.vector().array()[0]
            values_c[i] = vs.vector().array()[16]
            i=i+1

    # Save time steps
    np.savetxt("recorded_times_single_cell__{0}_{1}_{2}_{3}_{4}_{5}_{6}.txt".format(factor_gna, factor_gcal, factor_gk1, factor_gkr, factor_gto, factor_gf, factor_gks), times)
    np.savetxt("recorded_v_single_cell__{0}_{1}_{2}_{3}_{4}_{5}_{6}.txt".format(factor_gna, factor_gcal, factor_gk1, factor_gkr, factor_gto, factor_gf, factor_gks), values_v)
    np.savetxt("recorded_cai_single_cell__{0}_{1}_{2}_{3}_{4}_{5}_{6}.txt".format(factor_gna, factor_gcal, factor_gk1, factor_gkr, factor_gto, factor_gf, factor_gks), values_c)
    return vs

if __name__ == "__main__":

    i = main(1.0,1.0,1.0,1.0,1.0,1.0,1.0)

    i = main(0.5,1.0,1.0,1.0,1.0,1.0,1.0)
    i = main(1.5,1.0,1.0,1.0,1.0,1.0,1.0)
 
    i = main(1.0,0.5,1.0,1.0,1.0,1.0,1.0)
    i = main(1.0,1.5,1.0,1.0,1.0,1.0,1.0)

    i = main(1.0,1.0,0.5,1.0,1.0,1.0,1.0)
    i = main(1.0,1.0,1.5,1.0,1.0,1.0,1.0)

    i = main(1.0,1.0,1.0,0.5,1.0,1.0,1.0)
    i = main(1.0,1.0,1.0,1.5,1.0,1.0,1.0)

    i = main(1.0,1.0,1.0,1.0,0.5,1.0,1.0)
    i = main(1.0,1.0,1.0,1.0,1.5,1.0,1.0)

    i = main(1.0,1.0,1.0,1.0,1.0,0.5,1.0)
    i = main(1.0,1.0,1.0,1.0,1.0,1.5,1.0)

    i = main(1.0,1.0,1.0,1.0,1.0,1.0,0.5)
    i = main(1.0,1.0,1.0,1.0,1.0,1.0,1.5)


