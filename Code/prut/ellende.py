import math
from cbcbeat import *

parameters["adjoint"]["test_derivative"] = True

# For computing faster
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["cpp_optimize"] = True
flags = "-O3 -ffast-math -march=native"
parameters["form_compiler"]["cpp_optimize_flags"] = flags

class Stimulus(Expression):
    "Some self-defined stimulus."
    def __init__(self, **kwargs):
        self.t = kwargs["time"]
    def eval(self, value, x):
        if float(self.t) >= 2 and float(self.t) <= 11:
            v_amp = 125
            value[0] = 0.05*v_amp
        else:
            value[0] = 0.0


def forward(Gna):
    "Solve a single cell model on some time frame."

        # Choose your favorite cell model
    params = Tentusscher_2004_mcell.default_parameters()
    model = Tentusscher_2004_mcell(params={"g_Na":Gna})
    
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
    dt = 0.05
    interval = (0.0, 5.0)

    solutions = solver.solve(interval, dt)

    for ((t0, t1), fields) in solutions:
        vs = fields

    return split(vs)[0]

if __name__ == "__main__":

    Gna = Constant(14.838)              # initial guess.
    v = forward(Gna)

    #split(vs)[0]

    # Define functional of interest
    J = Functional(inner(v, v)*dx*dt[FINISH_TIME])

    # Indicate the control parameter of interest
    dJdnu = compute_gradient(J, Control(Gna))
    Jnu = assemble(inner(v, v)*dx) # current value

    parameters["adjoint"]["stop_annotating"] = True # stop registering equations

    def Jhat(Gna): # the functional as a pure function of Gna
         v = forward(Gna)
         #v=split(vs)[0]
         return assemble(inner(v, v)*dx)

    conv_rate = taylor_test(Jhat, Control(Gna), Jnu, dJdnu, seed=(0.001))

