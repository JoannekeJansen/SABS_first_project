# Main Ionic currents
# ===========================================================
#
# Save the main ionic currents during the action potential.
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

def main():
    #model = Paci2013_atrial()
    model = Paci2013_ventricular_adjusted()

    # Set initial conditions
    model.set_initial_conditions(V=vsinit((0.0))[0], m=vsinit((0.0))[1], h=vsinit((0.0))[2], \
    j=vsinit((0.0))[3], d=vsinit((0.0))[4], f1=vsinit((0.0))[5], f2=vsinit((0.0))[6], \
    fCa=vsinit((0.0))[7], Xr1=vsinit((0.0))[8], Xr2=vsinit((0.0))[9], Xs=vsinit((0.0))[10], \
    Xf=vsinit((0.0))[11], q=vsinit((0.0))[12], r=vsinit((0.0))[13], Nai=vsinit((0.0))[14], \
    g=vsinit((0.0))[15], Cai=vsinit((0.0))[16], Ca_SR=vsinit((0.0))[17])
    
    time = Constant(0.0)
    model.stimulus = Stimulus(time=time, degree=0)

    # Initialize solver
    params = SingleCellSolver.default_parameters()
    #params["scheme"] = "GRL1"
    solver = SingleCellSolver(model, time, params)

    # Assign initial conditions
    (vs_, vs) = solver.solution_fields()
    vs_.assign(model.initial_conditions())

    # Solve and extract values
    dt = 0.0001
    T = 1.0
    interval = (0.0, T)

    # Create arrays to save time, v and cai
    times = np.zeros(int(T/dt)+1)
    Ina = np.zeros(int(T/dt)+1)
    Ical = np.zeros(int(T/dt)+1)
    Ikr = np.zeros(int(T/dt)+1)
    Ik1 = np.zeros(int(T/dt)+1)
    Iks = np.zeros(int(T/dt)+1)
    If = np.zeros(int(T/dt)+1)
    Ito = np.zeros(int(T/dt)+1)
    v_values = np.zeros(int(T/dt)+1)
    i=0

    # Extract default parameters
    param = Paci2013_ventricular_adjusted.default_parameters()

    # Assign parameters
    g_Na_factor = param["g_Na_factor"]
    g_CaL_factor = param["g_CaL_factor"]
    g_Kr_factor = param["g_Kr_factor"]
    g_K1_factor = param["g_K1_factor"]
    PkNa = param["PkNa"]
    g_Na = param["g_Na"]
    g_CaL = param["g_CaL"]
    g_Kr = param["g_Kr"]
    g_Ks = param["g_Ks"]
    g_K1 = param["g_K1"]
    E_f = param["E_f"]
    g_f = param["g_f"]
    g_b_Na = param["g_b_Na"]
    g_b_Ca = param["g_b_Ca"]
    Km_K = param["Km_K"]
    Km_Na = param["Km_Na"]
    PNaK = param["PNaK"]
    KmCa = param["KmCa"]
    KmNai = param["KmNai"]
    Ksat = param["Ksat"]
    alpha = param["alpha"]
    gamma = param["gamma"]
    kNaCa = param["kNaCa"]
    KPCa = param["KPCa"]
    g_PCa = param["g_PCa"]
    g_to = param["g_to"]
    Cao = param["Cao"]
    F = param["F"]
    Ki = param["Ki"]
    Ko = param["Ko"]
    Nao = param["Nao"]
    R = param["R"]
    T = param["T"]

    for ((t0, t1), fields) in solver.solve(interval, dt):
        print vs.vector().array()[0]
        V=vs.vector().array()[0]
        m=vs.vector().array()[1] 
        h=vs.vector().array()[2]
        j=vs.vector().array()[3]
        d=vs.vector().array()[4] 
        f1=vs.vector().array()[5] 
        f2=vs.vector().array()[6] 
        fCa=vs.vector().array()[7] 
        Xr1=vs.vector().array()[8] 
        Xr2=vs.vector().array()[9] 
        Xs=vs.vector().array()[10] 
        Xf=vs.vector().array()[11] 
        q=vs.vector().array()[12]
        r=vs.vector().array()[13] 
        Nai=vs.vector().array()[14]
        g=vs.vector().array()[15]
        Cai=vs.vector().array()[16] 
        Ca_SR=vs.vector().array()[17]

        # Expressions for the Electric potentials component
        E_Na = R*T*ufl.ln(Nao/Nai)/F
        E_K = R*T*ufl.ln(Ko/Ki)/F
        E_Ks = R*T*ufl.ln((Ko + Nao*PkNa)/(Ki + PkNa*Nai))/F
        E_Ca = 0.5*R*T*ufl.ln(Cao/Cai)/F

        # Expressions for the i_Na component
        i_Na = g_Na_factor*g_Na*(m*m*m)*(-E_Na + V)*h*j

        # Expressions for the i_CaL component
        i_CaL = 4*g_CaL_factor*g_CaL*(F*F)*(-0.341*Cao +\
            Cai*ufl.exp(2*F*V/(R*T)))*V*d*f1*f2*fCa/(R*T*(-1 +\
            ufl.exp(2*F*V/(R*T))))

        # Expressions for the i_Kr component
        i_Kr = 0.430331482912*g_Kr_factor*g_Kr*ufl.sqrt(Ko)*(-E_K + V)*Xr1*Xr2

        # Expressions for the i_Ks component
        i_Ks = g_Ks*(Xs*Xs)*(1 + 0.6/(1 +\
            6.48182102606e-07*ufl.elem_pow(1.0/Cai, 1.4)))*(-E_Ks +\
            V)

        # Expressions for the i_K1 component
        alpha_K1 = 3.91/(1 + 2.44592857399e-52*ufl.exp(594.2*V - 594.2*E_K))
        beta_K1 = (0.00277806676906*ufl.exp(588.6*V - 588.6*E_K) -\
            1.5394838221*ufl.exp(0.2*V - 0.2*E_K))/(1 + ufl.exp(454.7*V -\
            454.7*E_K))
        XK1_inf = alpha_K1/(alpha_K1 + beta_K1)
        i_K1 = 0.430331482912*g_K1_factor*g_K1*ufl.sqrt(Ko)*(-E_K + V)*XK1_inf

        # Expressions for the i_f component
        i_f = g_f*(-E_f + V)*Xf

        # Expressions for the i_b Na component
        i_b_Na = g_b_Na*(-E_Na + V)

        # Expressions for the i_b Ca component
        i_b_Ca = g_b_Ca*(-E_Ca + V)

        # Expressions for the i_NaK component
        i_NaK = Ko*PNaK*Nai/((Km_K + Ko)*(Km_Na + Nai)*(1 +\
            0.0353*ufl.exp(-F*V/(R*T)) + 0.1245*ufl.exp(-0.1*F*V/(R*T))))

        # Expressions for the i_NaCa component
        i_NaCa = kNaCa*(Cao*(Nai*Nai*Nai)*ufl.exp(F*gamma*V/(R*T)) -\
            alpha*(Nao*Nao*Nao)*Cai*ufl.exp(F*(-1 + gamma)*V/(R*T)))/((1 +\
            Ksat*ufl.exp(F*(-1 + gamma)*V/(R*T)))*(Cao +\
            KmCa)*((KmNai*KmNai*KmNai) + (Nao*Nao*Nao)))

        # Expressions for the i_PCa component
        i_PCa = g_PCa*Cai/(KPCa + Cai)

        # Expressions for the i_to component
        i_to = g_to*(-E_K + V)*q*r

        times[i] = t1
        v_values[i] = vs.vector().array()[0]
        Ina[i] = i_Na
        Ical[i] = i_CaL
        Ikr[i] = i_Kr
        Ik1[i] = i_K1
        Iks[i] = i_Ks
        If[i] = i_f
        Ito[i] = i_to
        i=i+1

    # Save time steps
    np.savetxt("recorded_single_cell_times.txt", times)
    np.savetxt("recorded_single_cell_v_values.txt", v_values)
    np.savetxt("recorded_single_cell_Ina.txt", Ina)
    np.savetxt("recorded_single_cell_Ical.txt", Ical)
    np.savetxt("recorded_single_cell_Ikr.txt", Ikr)
    np.savetxt("recorded_single_cell_Ik1.txt", Ik1)
    np.savetxt("recorded_single_cell_Iks.txt", Iks)
    np.savetxt("recorded_single_cell_If.txt", If)
    np.savetxt("recorded_single_cell_Ito.txt", Ito)
    return vs

if __name__ == "__main__":

    i = main()

