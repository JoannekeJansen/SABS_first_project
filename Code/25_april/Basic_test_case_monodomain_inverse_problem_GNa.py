# Solver for the inverse problem of Basic_test_case_monodomain.py
# ===========================================================
#

# Import the cbcbeat module
from cbcbeat import *
from mshr import *
import numpy as np

# Turn on FFC/FEniCS optimizations
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["cpp_optimize"] = True
flags = ["-O3", "-ffast-math", "-march=native"]
parameters["form_compiler"]["cpp_optimize_flags"] = " ".join(flags)
parameters["form_compiler"]["quadrature_degree"] = 3

# Define a [0,5]^2 domain
from mshr import *
domain = Rectangle(Point(0.0, 0.0), Point(5.0, 5.0))
mesh = generate_mesh(domain, 100)

# Averaged nominal conductivities, surface to volume ratio and membrane
# capacitance, as found in Sepulveda, Roth, & Wikswo. (1989), Table 1
sigma_l = 0.15         # mS / mm
sigma_t = 0.02          # mS / mm
beta = 200.0            # mm^{-1}
C_m = 0.2               # mu F / mm^2

#def forward(sigma_t):
def forward(GNa):
    # Define time
    time = Constant(0.0)

    # We use the Grandi cell model
    cell_model = Grandi_pasqualini_bers_2010()
    cell_model.set_parameters(GNa=GNa)
    #cell_model = FitzHughNagumoManual() # For fast tests

    # Define the external stimulus
    p = Expression('10', degree=1, domain=mesh)
    stimulus = Expression('(x[0] > 2.0 and x[0] < 3.0 and x[1] > 2.0 \
    and x[1] < 3.0 and t < 3.0 ? p : 0)', p=p, t=time, degree=1)

    # Scale conducitivites by 1/(C_m * chi)
    M_l = sigma_l/(C_m*beta) # mm^2 / ms
    M_t = sigma_t/(C_m*beta) # mm^2 / ms

    # Define the conductivity (tensor)
    M = as_tensor(((M_l, 0.0), (0.0, M_t)))

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
    vs_.assign(cell_model.initial_conditions(), solver.VS)
    vs.assign(cell_model.initial_conditions())

    # Time stepping parameters
    h = 5  # Time step size
    T = 50.0 # Final time
    interval = (0.0, T)

    # Solve forward problem
    for (timestep, fields) in solver.solve(interval, h):
        print "(t_0, t_1) = ", timestep
        # Extract the components of the field (vs_ at previous timestep,
        # current vs, current vur)
        (vs_, vs, vur) = fields

    return vs_, vs, vur

if __name__ == "__main__":
    GNa= Constant(20)                # initial guess, GNa=23 in the test problem.
    ctrl1 = GNa
    (vs_, vs, vur) = forward(GNa)   # solve the forward problem once. 

    # Load recorded data and define functional of interest
    times = np.loadtxt("Results/recorded_times.txt")
    hdf_vs = HDF5File(mesh.mpi_comm(), "Results/Basic_test_case_monodomain_synthetic_observations_vs.h5", "r")
    attr_vs = hdf_vs.attributes("vs")
    N = attr_vs['count']
    Q = vs.function_space()
    vs_obs = {}
    I = 0
    for i in range(N):
      vs_obs[i] = Function(Q, annotate=False)
      dataset_vs = "vs/vector_%d"%i
      hdf_vs.read(vs_obs[i], dataset_vs)
      I = I + inner(split(vs)[0] - split(vs_obs[i])[0], split(vs)[0] - split(vs_obs[i])[0])*dx*dt[times[i]]
      I = I + inner(split(vs)[38] - split(vs_obs[i])[38], split(vs)[38] - split(vs_obs[i])[38])*dx*dt[times[i]]
    del hdf_vs
    J=Functional(I/N)

    # Define the reduced functional and solve the optimisation problem:
    rf = ReducedFunctional(J, [Control(ctrl1)])
    rf.taylor_test(ctrl1, seed=0.1)
    #assert rf.taylor_test(ctrl1, seed=1) > 1.5
    opt_ctrls = minimize(rf, tol=1e-02, options={"maxiter": 10})
    #opt_ctrls = minimize(rf)
    #print("ctrl1 = %f" %float(opt_ctrls))
    # print "seed=0.000001:"
    # rf.taylor_test(ctrl1, seed=0.000001)
    # print "seed=0.00001:"
    # rf.taylor_test(ctrl1, seed=0.00001)
    # print "seed=0.0001:"
    # rf.taylor_test(ctrl1, seed=0.0001)
    # print "seed=0.001:"
    # rf.taylor_test(ctrl1, seed=0.001)
    # print "seed=0.01:"
    # rf.taylor_test(ctrl1, seed=0.01)
    # print "seed=0.1:"
    # rf.taylor_test(ctrl1, seed=0.1)
    # print "seed=1:"
    # rf.taylor_test(ctrl1, seed=1)
    # print "seed=10:"
    # rf.taylor_test(ctrl1, seed=10)
    # print "seed=100:"
    # rf.taylor_test(ctrl1, seed=100)
    # print "seed=1000:"
    # rf.taylor_test(ctrl1, seed=1000)

# seed=0.000001
# Taylor remainder without gradient information: [0.0001020319372218026, 0.0005180796911079533, 0.0014478646307964915, 0.001404236378665491, 0.00031983884997366374]
# Convergence orders for Taylor remainder without gradient information (should all be 1): [-2.3441532268206755, -1.4826807871623733, 0.044140914263333976, 2.134368713693494]
# Taylor remainder with gradient information: [8.864614525506544e-05, 0.0005113867951245847, 0.0014512110787881758, 0.001405909602661333, 0.0003206754619715848]
# Convergence orders for Taylor remainder with gradient information (should all be 2): [-2.5282851056860154, -1.5047705614590123, 0.04575353997692158, 2.1323179664677183]

# seed=0.00001:
# Taylor remainder without gradient information: [0.001731591581509928, 0.00044477077502591555, 0.0005608998002788468, 0.00044808342148483504, 9.650497845692385e-05]
# Convergence orders for Taylor remainder without gradient information (should all be 1): [1.9609647919596074, -0.33468107466751906, 0.3239757197141785, 2.215092075062728]
# Taylor remainder with gradient information: [0.0018654559778182693, 0.0005117029731800862, 0.0005943658993559321, 0.0004648164710233777, 8.813845368765253e-05]
# Convergence orders for Taylor remainder with gradient information (should all be 2): [1.8661497938812242, -0.2160447303370932, 0.35469015396314035, 2.3988177012053367]

# seed=0.0001:
# Taylor remainder without gradient information: [0.0010400387025795155, 0.002633435596272271, 0.0021765954646699726, 0.0013026772191730984, 0.0003728982584041418]
# Convergence orders for Taylor remainder without gradient information (should all be 1): [-1.34030896129184, 0.2748728795565626, 0.7405936442115062, 1.8046256891402748]
# Taylor remainder with gradient information: [0.002378841425384936, 0.0033028369576749814, 0.0025112961453713275, 0.0014700275595237759, 0.000289223088228803]
# Convergence orders for Taylor remainder with gradient information (should all be 2): [-0.4734466494280939, 0.39527358481228414, 0.7725889664793514, 2.345588572321628]

# seed=0.001:
# Taylor remainder without gradient information: [0.0043253525770756696, 2.3412955097334986e-05, 0.0005043056678246671, 0.0013182049632693982, 0.0003871592982198635]
# Convergence orders for Taylor remainder without gradient information (should all be 1): [7.529366889636296, -4.428919496912437, -1.386204362843609, 1.7675755114020895]
# Taylor remainder with gradient information: [0.009059474499799873, 0.006715826493535106, 0.0028419011013942185, 0.00035489842134004464, 0.0004493923940848579]
# Convergence orders for Taylor remainder with gradient information (should all be 2): [0.4318624102488313, 1.2407086086271004, 3.0013782873132433, -0.3405695502612113]

# seed=0.01:
# Taylor remainder without gradient information: [0.04612080523088391, 0.02388384603664173, 0.010077083381336838, 0.006645434201342226, 0.0028217039304796465]
# Convergence orders for Taylor remainder without gradient information (should all be 1): [0.9493825277088591, 1.2449570349332635, 0.6006427678941215, 1.2357968458595647]
# Taylor remainder with gradient information: [0.08773075999967764, 0.04304193657863904, 0.02338580792630355, 0.010086011452477967, 0.00554401889643045]
# Convergence orders for Taylor remainder with gradient information (should all be 2): [1.0273397761912788, 0.8801097188513926, 1.2132775028775173, 0.863351688633085]

# seed=0.1:
# Taylor remainder without gradient information: [0.4940282456910077, 0.240197290894391, 0.11672765846302102, 0.059327796266956945, 0.02697608599717194]
# Convergence orders for Taylor remainder without gradient information (should all be 1): [1.0403736497087575, 1.0410734332386131, 0.976366346269002, 1.137027154262772]
# Taylor remainder with gradient information: [0.8456594114460905, 0.42964653767415806, 0.2181942558212535, 0.10813316087518032, 0.056754392573896695]
# Convergence orders for Taylor remainder with gradient information (should all be 2): [0.9769264672030992, 0.9775371472178503, 1.0128041038659288, 0.9300050569866347]

# seed=1:
# Taylor remainder without gradient information: [6.908552880808297, 2.9266809348689424, 1.3312069065984753, 0.6304268131175528, 0.30539720123691305]
# Convergence orders for Taylor remainder without gradient information (should all be 1): [1.2391180716342562, 1.136530648871989, 1.0783340221652953, 1.0456420555061368]
# Taylor remainder with gradient information: [6.563141726314522, 3.809166368692467, 2.0367167451822294, 1.0535350127727996, 0.5365837117082631]
# Convergence orders for Taylor remainder with gradient information (should all be 2): [0.7849112871443796, 0.9032299477032376, 0.951007091631543, 0.9733630943381869]

# seed=10:
# Taylor remainder without gradient information: [218.2566111803762, 73.32426002583018, 25.318792263852238, 9.607060385910902, 3.996477124168002]
# Convergence orders for Taylor remainder without gradient information (should all be 1): [1.5736628432189868, 1.5340820184994324, 1.3980416272225062, 1.2653662245499109]
# Taylor remainder with gradient information: [75.60149844710844, 1.9967036591962994, 10.344985919464705, 8.22482870574757, 4.919467421661233]
# Convergence orders for Taylor remainder with gradient information (should all be 2): [5.242722693254494, -2.3732395453747026, 0.33087414354039363, 0.7414834941047641]

# seed=100:
# Taylor remainder without gradient information: [2498.2473162997057, 1450.6865486267739, 778.4878462402919, 323.954327100433, 117.6960208149046]
# Convergence orders for Taylor remainder without gradient information (should all be 1): [0.7841804765058672, 0.8979894072134498, 1.2648840879101013, 1.4607248830810171]
# Taylor remainder with gradient information: [303.299990776467, 353.2128858651545, 229.7510148594822, 49.58591141002813, 19.488187030297823]
# Convergence orders for Taylor remainder with gradient information (should all be 2): [-0.2197925220646829, 0.6204667419319163, 2.212069055797404, 1.3473303891835522]

# seed=1000:
# Taylor remainder without gradient information: [6787.629038620802, 5323.514212585675, 4012.3385089545504, 2887.7818657654393, 1770.980432865043]
# Convergence orders for Taylor remainder without gradient information (should all be 1): [0.35052879324522707, 0.4079355998621837, 0.4744815573767478, 0.7054134973686024]
# Taylor remainder with gradient information: [74343.179747738, 35241.89018059373, 16270.36368763515, 7253.569232529411, 3299.695116282382]
# Convergence orders for Taylor remainder with gradient information (should all be 2): [1.0769090907940868, 1.115044804941725, 1.1654835239125743, 1.1363583414169185]



    

#   