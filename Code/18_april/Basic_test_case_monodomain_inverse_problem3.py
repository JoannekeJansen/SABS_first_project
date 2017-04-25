# Solver for the inverse problem of Basic_test_case_monodomain.py
# ===========================================================
#

# Import the cbcbeat module
from cbcbeat import *
from mshr import *
import numpy.random

set_log_level(20)


# Turn on FFC/FEniCS optimizations
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["cpp_optimize"] = True
flags = ["-O3", "-ffast-math", "-march=native"]
parameters["form_compiler"]["cpp_optimize_flags"] = " ".join(flags)
parameters["form_compiler"]["quadrature_degree"] = 3

# Define a [0,10]^2 domain
from mshr import *
domain = Rectangle(Point(0.0, 0.0), Point(10.0, 10.0))
mesh = generate_mesh(domain, 100)

# Define the conductivity (tensor)
#M = as_matrix(((0.255, 0.0), (0.0, 0.0775)))
#M=0.255

Q = FunctionSpace(mesh, "CG", 1)

def forward(M):
	# Define time
	time = Constant(0.0)
	
	#v = Function(Q)

	# We use the Grandi cell model
	cell_model = Grandi_pasqualini_bers_2010()
	cell_model.set_parameters(GNa=GNa)

	#cell_model = FitzHughNagumoManual()
	#cell_model.set_parameters(v_rest=GNa)

	# Define the external stimulus
	p = Expression('25', degree=1, domain=mesh)
	stimulus = Expression('(x[0] > 4.5 and x[0] < 5.5 and x[1] > 4.5 \
		and x[1] < 5.5 and t > 0 and t < 100 ? p : 0)', p=p, t=time, degree=1)

	# Collect this information into the CardiacModel class
	cardiac_model = CardiacModel(mesh, time, M, 'none', cell_model, stimulus)

	# Customize and create a splitting solver
	ps = SplittingSolver.default_parameters()
	ps["theta"] = 0.5                        # Second order splitting scheme
	ps["pde_solver"] = "monodomain"          # Use Monodomain model for the PDEs
	ps["CardiacODESolver"]["scheme"] = "GRL1" #  1st order Rush-Larsen for the ODEs
	ps["MonodomainSolver"]["linear_solver_type"] = "iterative"
	ps["MonodomainSolver"]["algorithm"] = "cg"
	ps["MonodomainSolver"]["preconditioner"] = "petsc_amg"

	solver = SplittingSolver(cardiac_model, params=ps)

	# Extract the solution fields and set the initial conditions
	(vs_, vs, vur) = solver.solution_fields()
	vs_.assign(cell_model.initial_conditions())

	# Time stepping parameters
	dt = 0.5
	T = 5.5
	interval = (0.0, T)

	# Solve forward problem
	for (timestep, fields) in solver.solve(interval, dt):
	    #print "(t_0, t_1) = ", timestep

	    # Extract the components of the field (vs_ at previous timestep,
	    # current vs, current vur)
	    (vs_, vs, vur) = fields
	#    v.assign(vur)

	return vur

# # Load synthetic data
# Q = FunctionSpace(mesh, "CG", 1)
# h = Function(Q)
# g = Function(Q)
# j = Function(Q)
# hdf_v = HDF5File(mesh.mpi_comm(), "Results/Basic_test_case_monodomain_synthetic_observations_v.h5", "r")
# hdf_c = HDF5File(mesh.mpi_comm(), "Results/Basic_test_case_monodomain_synthetic_observations_c.h5", "r")
# attr_v = hdf_v.attributes("v")
# attr_c = hdf_c.attributes("c")
# N = attr_v['count']
# for timestep in range(N):
#     dataset_v = "v/vector_%d"%timestep
#     dataset_c = "c/vector_%d"%timestep
#     print 'Retrieving time step t=', timestep
#     hdf_v.read(h, dataset_v)
#     g.assign(h)
#     hdf_c.read(h, dataset_c)
#     j.assign(h)
# del hdf_v, hdf_c

if __name__ == "__main__":
    GNa = Constant(23)	      	  # initial guess.
    v = forward(GNa)                  # solve the forward problem once. 

    #success = replay_dolfin(tol=0.0, stop=True)
    # adj_html("forward.html", "forward")
    # adj_html("adjoint.html", "adjoint")
    
    # Define functional of interest
    J = Functional(inner(v, v)*dx*dt[FINISH_TIME])

    # Indicate the control parameter of interest
    dJdnu = compute_gradient(J, Control(GNa))
    Jnu = assemble(inner(v, v)*dx) # current value

    parameters["adjoint"]["test_derivative"] = True
    parameters["adjoint"]["stop_annotating"] = True # stop registering equations

    def Jhat(GNa): # the functional as a pure function of GNa
         v = forward(GNa)
         return assemble(inner(v, v)*dx)


    conv_rate = taylor_test(Jhat, Control(GNa), Jnu, dJdnu, seed=(1000))
    print conv_rate
    print "when seed was ", 1000
    conv_rate = taylor_test(Jhat, Control(GNa), Jnu, dJdnu, seed=(100))
    print "when seed was ", 100
    conv_rate = taylor_test(Jhat, Control(GNa), Jnu, dJdnu, seed=(10))
    print "when seed was ", 10
    conv_rate = taylor_test(Jhat, Control(GNa), Jnu, dJdnu, seed=(1))
    print "when seed was ", 1
    conv_rate = taylor_test(Jhat, Control(GNa), Jnu, dJdnu, seed=(0.1))
    print "when seed was ", 0.1
    conv_rate = taylor_test(Jhat, Control(GNa), Jnu, dJdnu, seed=(0.01))
    print "when seed was ", 0.01
    conv_rate = taylor_test(Jhat, Control(GNa), Jnu, dJdnu, seed=(0.001))
    print "when seed was ", 0.001
    conv_rate = taylor_test(Jhat, Control(GNa), Jnu, dJdnu, seed=(0.0001))
    print "when seed was ", 0.0001
    conv_rate = taylor_test(Jhat, Control(GNa), Jnu, dJdnu, seed=(0.0001))
    print "when seed was ", 0.00001

print "Success!"

# Output:
# Taylor remainder without gradient information: [48.20048773754388, 44.47637924703304, 61.54938862065319, 48.90004040184431, 36.26104794675484]
# Convergence orders for Taylor remainder without gradient information (should all be 1): [0.11600839978583023, -0.46870518099492703, 0.33190886912252904, 0.4314150374817749]
# Taylor remainder with gradient information: [939.7235459657641, 449.48563760462093, 185.4316198051738, 74.59046381106918, 25.484204159701903]
# Convergence orders for Taylor remainder with gradient information (should all be 2): [1.0639613796179586, 1.27738774476256, 1.3138241706146465, 1.5493878972552297]
# 1.06396137962
# when seed was  1000
# Taylor remainder without gradient information: [44.36995293048676, 31.45715211587958, 20.084189843619242, 15.937913833651692, 13.053287267917767]
# Convergence orders for Taylor remainder without gradient information (should all be 1): [0.49619495632895305, 0.6473277992576171, 0.33359746454851263, 0.288049629633938]
# Taylor remainder with gradient information: [54.42245043984403, 17.939049569285814, 4.613910998963455, 3.588863412360343, 6.878762057272093]
# Convergence orders for Taylor remainder with gradient information (should all be 2): [1.6010984605229326, 1.9590413764776715, 0.3624631580752703, -0.9386219353147917]
# when seed was  100
# Taylor remainder without gradient information: [14.866653769044206, 12.386162636103109, 3.543618398020044, 3.2777575738728046, 3.128309537889436]
# Convergence orders for Taylor remainder without gradient information (should all be 1): [0.2633506632568921, 1.8054341368239208, 0.1125140973563783, 0.06732588472297753]
# Taylor remainder with gradient information: [4.987413432011126, 7.446542467586569, 1.073808313761774, 2.0428525317436694, 2.5108570168248683]
# Convergence orders for Taylor remainder with gradient information (should all be 2): [-0.5782789175526704, 2.7938342379733836, -0.9278485835420254, -0.2975948118924035]
# when seed was  10
# Taylor remainder without gradient information: [3.2193250907585025, 3.0970652197720483, 0.01943372213281691, 0.00995469547342509, 0.0050368859665468335]
# Convergence orders for Taylor remainder without gradient information (should all be 1): [0.055856504570332074, 7.31619570802833, 0.9651131582875645, 0.9828451135552808]
# Taylor remainder with gradient information: [2.2314010570551943, 2.603103202920394, 0.22754728629301008, 0.1135358087394884, 0.05670836613990991]
# Convergence orders for Taylor remainder with gradient information (should all be 2): [-0.22228267290995335, 3.5159942240452224, 1.0030189914402576, 1.0015138938990744]
# when seed was  1
# Taylor remainder without gradient information: [0.008001836715266109, 0.0040390705689787865, 0.002029058989137411, 0.0010169490706175566, 0.0005090739578008652]
# Convergence orders for Taylor remainder without gradient information (should all be 1): [0.9863078370481568, 0.9932125444629328, 0.9965633779482871, 0.9983002599628034]
# Taylor remainder with gradient information: [0.09079056665506469, 0.04535713111618661, 0.022669041853445288, 0.011332101350673793, 0.0056654512528448095]
# Convergence orders for Taylor remainder with gradient information (should all be 2): [1.0012130147652258, 1.0006059769171016, 1.0003080044602386, 1.0001526335807942]
# when seed was  0.1
# Taylor remainder without gradient information: [0.0008139272686094046, 0.0004073419841006398, 0.00020374578889459372, 0.00010188866872340441, 5.102483555674553e-05]
# Convergence orders for Taylor remainder without gradient information (should all be 1): [0.9986593647366251, 0.9994702769126494, 0.99977662646865, 0.9977220820941154]
# Taylor remainder with gradient information: [0.009065313068423675, 0.0045322781844159, 0.002266064295363676, 0.0011330163734057305, 0.0005664276855078219]
# Convergence orders for Taylor remainder with gradient information (should all be 2): [1.00012042963167, 1.000047618796903, 1.000020085598082, 1.0002050225046604]
# when seed was  0.01
# Taylor remainder without gradient information: [8.151587098836899e-05, 4.078890196979046e-05, 2.0420411601662636e-05, 1.0199262760579586e-05, 5.095149390399456e-06]
# Convergence orders for Taylor remainder without gradient information (should all be 1): [0.9989043061796805, 0.9981647245269958, 1.001547073248752, 1.0012685206578027]
# Taylor remainder with gradient information: [0.000906408162714939, 0.00045317311488186354, 0.00022656059682416436, 0.00011329124145233391, 5.665010271605729e-05]
# Convergence orders for Taylor remainder with gradient information (should all be 2): [1.0000985796663757, 1.000165302828263, 0.9998606403087947, 0.9998858536863525]
# when seed was  0.001
# Taylor remainder without gradient information: [8.138129487633705e-06, 4.121684469282627e-06, 2.025510184466839e-06, 1.0497169569134712e-06, 4.790490493178368e-07]
# Convergence orders for Taylor remainder without gradient information (should all be 1): [0.981463170365284, 1.0249487262218149, 0.9482849634558447, 1.13175509096786]
# Taylor remainder with gradient information: [9.06542738826971e-05, 4.5274517215882774e-05, 2.267259065811586e-05, 1.1299333464377879e-05, 5.695476161327838e-06]
# Convergence orders for Taylor remainder with gradient information (should all be 2): [1.0016757801217653, 0.997750008674555, 1.0047115763169046, 0.9883493057713338]
# when seed was  0.0001
# Taylor remainder without gradient information: [8.138129487633705e-06, 4.121684469282627e-06, 2.025510184466839e-06, 1.0497169569134712e-06, 4.790490493178368e-07]
# Convergence orders for Taylor remainder without gradient information (should all be 1): [0.981463170365284, 1.0249487262218149, 0.9482849634558447, 1.13175509096786]
# Taylor remainder with gradient information: [9.06542738826971e-05, 4.5274517215882774e-05, 2.267259065811586e-05, 1.1299333464377879e-05, 5.695476161327838e-06]
# Convergence orders for Taylor remainder with gradient information (should all be 2): [1.0016757801217653, 0.997750008674555, 1.0047115763169046, 0.9883493057713338]
# when seed was  1e-05


