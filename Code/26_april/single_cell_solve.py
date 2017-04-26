from cbcbeat import *

# Disable adjointing
parameters["adjoint"]["stop_annotating"] = True

# For computing faster
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["cpp_optimize"] = True
flags = "-O3 -ffast-math -march=native"
parameters["form_compiler"]["cpp_optimize_flags"] = flags

# # Define the external stimulus
# class Stimulus(Expression):
#     "Some self-defined stimulus."
#     def __init__(self, **kwargs):
#         self.t = kwargs["time"]
#     def eval(self, value, x):
#         if float(self.t) <= 3:
#             value[0] = 10.0
#         else:
#             value[0] = 0.0

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

def main():
    model = Grandi_pasqualini_bers_2010()
        
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
    interval = (0.0, 600.0)
    
    solutions = solver.solve(interval, dt)
    times = []
    values = []
    
    for ((t0, t1), vs) in solutions:
        #times.append(t1)
        #values.append(vs.vector().array())
        print "Current time: %g" % t1
        print "v = ", vs.vector().array()[0]

    return times, values

if __name__ == "__main__":

    (times, values) = main()