from dolfin import *
from dolfin_adjoint import *
from cbcbeat import *

def main(GNa):

    parameters["form_compiler"]["quadrature_degree"] = 2
    parameters["form_compiler"]["cpp_optimize"] = True
    #flags = ["-O3", "-ffast-math", "-march=native"]
    #parameters["form_compiler"]["cpp_optimize_flags"] = " ".join(flags)

    # Create cell model
    cell = FitzHughNagumoManual()
    cell.set_parameters(c_1=GNa)
    num_states = cell.num_states()

    # Create function spaces
    mesh = UnitSquareMesh(10, 10)
    V = FiniteElement("CG", mesh.ufl_cell(), 1)
    S  = FiniteElement("CG", mesh.ufl_cell(), 1)
    VS = FunctionSpace(mesh, V * S)
    #S = BasicSplittingSolver.state_space(mesh, num_states)

    # Create solution function and set its initial value
    vs = Function(VS, name="vs")
    vs.assign(project(cell.initial_conditions(), VS, annotate=False))
    (v, s) = split(vs)

    # Define the right-hand-side of the system of ODEs: Dt(u) = rhs(u)
    # Note that sign of the ionic current
    (w, r) = TestFunctions(VS)
    rhs = inner(cell.F(v, s), r) + inner(- cell.I(v, s), w)
    form = rhs*dP

    # In the beginning...
    time = Constant(0.0)

    # Create scheme
    scheme = BackwardEuler(form, vs, time)
    scheme.t().assign(float(time))  # FIXME: Should this be scheme or
                                    # solver and why is this needed?
    info(scheme)

    # Create solver
    solver = PointIntegralSolver(scheme)

    # Time step
    k = 0.1
    solver.step(k)

    return vs

if __name__ == "__main__":

    # Run forward
    GNa = Constant(0.26)              # initial guess.
    v = main(GNa)

    # Replay
    # print "Replaying"
    # success = replay_dolfin(tol=0.0, stop=True)
    # assert (success == True), "Replay failed." # !? How does this actually work

    # Define functional of interest
    J = Functional(inner(v, v)*dx*dt[FINISH_TIME])

    # Indicate the control parameter of interest
    dJdnu = compute_gradient(J, Control(GNa))
    Jnu = assemble(inner(v, v)*dx) # current value

    parameters["adjoint"]["test_derivative"] = True
    parameters["adjoint"]["stop_annotating"] = True # stop registering equations

    def Jhat(GNa): # the functional as a pure function of GNa
         v = main(GNa)
         return assemble(inner(v, v)*dx)

    conv_rate = taylor_test(Jhat, Control(GNa), Jnu, dJdnu, seed=(0.01))
