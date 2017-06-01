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

def main():
    #model = Paci2013_atrial()
    model = Paci2013_ventricular_adjusted()
    model.set_parameters()
        
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
    T = 800.0
    interval = (0.0, T)

    for ((t0, t1), fields) in solver.solve(interval, dt):
        print vs.vector().array()[0]

    return vs

if __name__ == "__main__":

  # Create HDF5 file to save initial conditions
  mesh = UnitIntervalMesh(1)
  hdf_vs = HDF5File(mesh.mpi_comm(), "initial_conditions_single_cell.h5", "w")
  vsinit = main()
  hdf_vs.write(vsinit,"vs",0)
  del hdf_vs 