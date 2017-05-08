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
    #model = Grandi_pasqualini_bers_2010_adjusted()
    #model = Paci_hyttinen_aaltosetala_severi_ventricularVersion()
    #model = Paci_hyttinen_aaltosetala_severi_atrialVersion()
    model = Paci2013_atrial()
    #model = Paci2013_ventricular()
    model.set_parameters()
        
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
    T = 5.0
    interval = (0.0, T)

    # Create arrays to save time, v and cai
    times = np.zeros(int(T/dt))
    values_v = np.zeros(int(T/dt))
    values_c = np.zeros(int(T/dt))
    i=0

    for ((t0, t1), fields) in solver.solve(interval, dt):
        if True:
        #if (abs(round(t1)-t1) < 0.1*dt) & (i < int(500)):
            print t1
            print vs.vector().array()[0]
            times[i] = t1
            values_v[i] = vs.vector().array()[0]
            values_c[i] = vs.vector().array()[16]
            i=i+1

    # Save time steps
    np.savetxt("single_cell_times_atrial_stim.txt", times)
    np.savetxt("single_cell_v_atrial_stim.txt", values_v)
    np.savetxt("single_cell_cai_atrial_stim.txt", values_c)
    # np.savetxt("single_cell_times_ventricular_stim.txt", times)
    # np.savetxt("single_cell_v_ventricular_stim.txt", values_v)
    # np.savetxt("single_cell_cai_ventricular_stim.txt", values_c)

    return vs

if __name__ == "__main__":

    vs = main()
