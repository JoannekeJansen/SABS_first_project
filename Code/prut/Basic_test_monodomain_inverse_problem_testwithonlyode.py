from cbcbeat import *

# For computing faster
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["cpp_optimize"] = True
flags = "-O3 -ffast-math -march=native"
parameters["form_compiler"]["cpp_optimize_flags"] = flags
parameters["form_compiler"]["quadrature_degree"] = 4

def forward(Gna):
    # Set-up domain in space and time
    N = 10
    mesh = UnitSquareMesh(N, N)
    time = Constant(0.0)
  
    # Choose your favorite cell model
    params = Tentusscher_2004_mcell.default_parameters()
    model = Tentusscher_2004_mcell(params={"g_Na":Gna})

    # You can set spatially varying cell model parameters e.g. as:
    model.set_parameters(K_mNa=Expression("40*sin(pi*x[0])", degree=4))

    # Add some stimulus
    stimulus = Expression("100*t", t=time, degree=0)

    Solver = CardiacODESolver
    params = Solver.default_parameters()
    solver = Solver(mesh, time, model, I_s=stimulus, params=params)

    # Set-up initial conditions
    (vs_, vs) = solver.solution_fields()
    vs_.assign(model.initial_conditions())

    # Set-up other solution parameters
    dt = 0.2
    interval = (0.0, 1.0)

    # Generator for solutions
    solutions = solver.solve(interval, dt)

    for ((t0, t1), fields) in solutions:
        vs = fields
    
    return vs


if __name__ == "__main__":
    Gna = Constant(14.838)              # initial guess.
    vs = forward(Gna)                  # solve the forward problem once. 
    v=split(vs)[0]

    #success = replay_dolfin(tol=0.0, stop=True)
    # adj_html("forward.html", "forward")
    # adj_html("adjoint.html", "adjoint")

    
    # Define functional of interest
    J = Functional(inner(v, v)*dx*dt[FINISH_TIME])

    # Indicate the control parameter of interest
    dJdnu = compute_gradient(J, Control(Gna))
    Jnu = assemble(inner(v, v)*dx) # current value

    parameters["adjoint"]["test_derivative"] = True
    parameters["adjoint"]["stop_annotating"] = True # stop registering equations

    def Jhat(Gna): # the functional as a pure function of Gna
         vs = forward(Gna)
         v=split(vs)[0]
         return assemble(inner(v, v)*dx)

    conv_rate = taylor_test(Jhat, Control(Gna), Jnu, dJdnu, seed=(0.00001))

