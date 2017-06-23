# Average the recorded observations and save as txt files.
# ===========================================================
#
# Joanneke E Jansen, April 2017

# Import the cbcbeat module
from cbcbeat import *
import numpy as np
import argparse

# Retrieve parameters
parser = argparse.ArgumentParser()
parser.add_argument("-gnaf", type=float, help="The factor of g_na")
parser.add_argument("-gk1f", type=float, help="The factor of g_k1")
parser.add_argument("-gkrf", type=float, help="The factor of g_kr")
parser.add_argument("-gcalf", type=float, help="The factor of g_cal")
parser.add_argument("-gna", type=int, help="The added percentage of noise to g_na")
parser.add_argument("-gk1", type=int, help="The added percentage of noise to g_k1")
parser.add_argument("-gkr", type=int, help="The added percentage of noise to g_kr")
parser.add_argument("-gcal", type=int, help="The added percentage of noise to g_cal")
parser.add_argument("-v", type=int, help="The added percentage of noise to the observed v")
parser.add_argument("-cai", type=int, help="The added percentage of noise to the observed cai")
args = parser.parse_args()

# Set parameter factors
global gnaf
global gk1f
global gkrf
global gcalf
gnaf = args.gnaf
gk1f = args.gk1f
gkrf = args.gkrf
gcalf = args.gcalf

# Set noise parameters
global noise_percentage_gna
global noise_percentage_gk1
global noise_percentage_gkr
global noise_percentage_gcal
global noise_percentage_v
global noise_percentage_cai

if args.gna:
  noise_percentage_gna = args.gna
else:
  noise_percentage_gna = 0
if args.gk1:
  noise_percentage_gk1 = args.gk1
else:
  noise_percentage_gk1 = 0
if args.gkr:
  noise_percentage_gkr = args.gkr
else:
  noise_percentage_gkr = 0
if args.gcal:
  noise_percentage_gcal = args.gcal
else:
  noise_percentage_gcal = 0
if args.v:
  noise_percentage_v = args.v
else:
  noise_percentage_v = 0
if args.cai:
  noise_percentage_cai = args.cai
else:
  noise_percentage_cai = 0
if args.gna or args.gk1 or args.gkr or args.gcal or args.v or args.cai:
   observed_times ='observed_times_{0}_{1}_{2}_{3}_12mm_strip_{4}_{5}_{6}_{7}_{8}_{9}.txt'.format(gnaf,gk1f,gkrf,gcalf, \
     noise_percentage_gna, noise_percentage_gk1, noise_percentage_gkr, noise_percentage_gcal, noise_percentage_v, noise_percentage_cai)
   observed_v ='observed_v_{0}_{1}_{2}_{3}_12mm_strip_{4}_{5}_{6}_{7}_{8}_{9}.h5'.format(gnaf,gk1f,gkrf,gcalf, \
     noise_percentage_gna, noise_percentage_gk1, noise_percentage_gkr, noise_percentage_gcal, noise_percentage_v, noise_percentage_cai)
   observed_cai ='observed_cai_{0}_{1}_{2}_{3}_12mm_strip_{4}_{5}_{6}_{7}_{8}_{9}.h5'.format(gnaf,gk1f,gkrf,gcalf, \
     noise_percentage_gna, noise_percentage_gk1, noise_percentage_gkr, noise_percentage_gcal, noise_percentage_v, noise_percentage_cai)
else:
  observed_times = "observed_times_{0}_{1}_{2}_{3}_12mm_strip.txt".format(gnaf,gk1f,gkrf,gcalf)
  observed_v = "observed_v_{0}_{1}_{2}_{3}_12mm_strip.h5".format(gnaf,gk1f,gkrf,gcalf)
  observed_cai = "observed_cai_{0}_{1}_{2}_{3}_12mm_strip.h5".format(gnaf,gk1f,gkrf,gcalf)

# For computing faster
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["cpp_optimize"] = True
flags = "-O3 -ffast-math -march=native"
parameters["form_compiler"]["cpp_optimize_flags"] = flags

# Define a [0,12]*[0,0.01] mm^2 rectangular domain
mesh = RectangleMesh(Point(0.0, 0.0), Point(12.0, 0.01), 200, 1)

# Load recorded data 
times = np.loadtxt(observed_times)
hdf_v = HDF5File(mesh.mpi_comm(), observed_v, "r")
hdf_cai = HDF5File(mesh.mpi_comm(), observed_cai, "r")
attr_v = hdf_v.attributes("v")
attr_cai = hdf_cai.attributes("cai")
N = attr_v['count']
Q = FunctionSpace(mesh, "CG", 1)
v_obs = Function(Q)
cai_obs = Function(Q)
gdim = mesh.geometry().dim()
dofmap = Q.dofmap()
dofs = dofmap.dofs()

# Get coordinates as len(dofs) x gdim array
dofs_x = Q.tabulate_dof_coordinates().reshape((-1, gdim))

# Create arrays for storage
values_v = np.zeros(2*N)
values_c = np.zeros(2*N)
local_values_v = np.zeros(N)
local_values_c = np.zeros(N)
number_of_added_values_v = np.zeros(2*N)
number_of_added_values_c = np.zeros(2*N)
summed_values=np.zeros(N)

# Set reference voltage peak
reference_peak=N

v_peak = 0
v_peak_time=0

# Extract waves from data
for dof, dof_x in zip(dofs, dofs_x):
  if dof_x[0]>5.0 and dof_x[0]<10.0:   
    for i in range(N):
        dataset_v = "v/vector_%d"%i
        dataset_cai = "cai/vector_%d"%i
        hdf_v.read(v_obs, dataset_v)
        hdf_cai.read(cai_obs, dataset_cai)
        v=v_obs.vector().array()
        cai=cai_obs.vector().array()
        local_values_v[i] = v[dof]
        local_values_c[i] = cai[dof]
        if dof_x[0] == 7.5 and dof_x[1] == 0 and v[dof]>v_peak:
            v_peak=v[dof]
            v_peak_time=times[i]

    local_shift = reference_peak - numpy.argmax(local_values_v)
    values_v[local_shift:local_shift+N] += local_values_v
    values_c[local_shift:local_shift+N] += local_values_c
    number_of_added_values_v[local_shift:local_shift+N] += 1
    number_of_added_values_c[local_shift:local_shift+N] += 1
del hdf_v, hdf_cai

# Find interval with most information
for i in range(N):
    summed_values[i]=np.sum(number_of_added_values_v[i:i+N]+number_of_added_values_c[i:i+N])
interval=numpy.argmax(summed_values)

# Truncate and average storage vectors
values_v=np.divide(values_v[interval:interval+N],number_of_added_values_v[interval:interval+N])
values_c=np.divide(values_c[interval:interval+N],number_of_added_values_c[interval:interval+N])

# Calculate characteristics
amplitude_v=max(values_v)-min(values_v)
amplitude_c=max(values_c)-min(values_c)
upstroke_velocity_c= max(values_c[5:(np.argmax(values_c))]-values_c[0:(np.argmax(values_c)-5)])*200 #mM/s   
upstroke_velocity_v= max(values_v[5:(np.argmax(values_v))]-values_v[0:(np.argmax(values_v)-5)])*200 #V/s  

temp=np.argmax(values_v>(min(values_v)+0.7*amplitude_v))
v_30_p=np.argmax(values_v[temp:-1]<(min(values_v)+0.7*amplitude_v))
temp=np.argmax(values_c>(min(values_c)+0.7*amplitude_c))
c_30_p=np.argmax(values_c[temp:-1]<(min(values_c)+0.7*amplitude_c))

temp=np.argmax(values_v>(min(values_v)+0.5*amplitude_v))
v_50_p=np.argmax(values_v[temp:-1]<(min(values_v)+0.5*amplitude_v))
temp=np.argmax(values_c>(min(values_c)+0.5*amplitude_c))
c_50_p=np.argmax(values_c[temp:-1]<(min(values_c)+0.5*amplitude_c))

temp=np.argmax(values_v>(min(values_v)+0.3*amplitude_v))
v_70_p=np.argmax(values_v[temp:-1]<(min(values_v)+0.3*amplitude_v))
temp=np.argmax(values_c>(min(values_c)+0.3*amplitude_c))
c_70_p=np.argmax(values_c[temp:-1]<(min(values_c)+0.3*amplitude_c))

temp=np.argmax(values_v>(min(values_v)+0.1*amplitude_v))
v_90_p=np.argmax(values_v[temp:-1]<(min(values_v)+0.1*amplitude_v))
temp=np.argmax(values_c>(min(values_c)+0.1*amplitude_c))
c_90_p=np.argmax(values_c[temp:-1]<(min(values_c)+0.1*amplitude_c))

# Save to file
np.savetxt('characteristics_observed_data_{0}_{1}_{2}_{3}_12mm_strip_{4}_{5}_{6}_{7}_{8}_{9}.txt'.format(gnaf,gk1f,gkrf,gcalf, \
     noise_percentage_gna, noise_percentage_gk1, noise_percentage_gkr, noise_percentage_gcal, noise_percentage_v, noise_percentage_cai), np.array([amplitude_v, \
      amplitude_c, upstroke_velocity_c, upstroke_velocity_v, v_30_p, c_30_p, v_50_p, c_50_p, v_70_p, c_70_p, v_90_p, c_90_p]))
np.savetxt('v_peak_time_{0}_{1}_{2}_{3}_12mm_strip_{4}_{5}_{6}_{7}_{8}_{9}.txt'.format(gnaf,gk1f,gkrf,gcalf, \
     noise_percentage_gna, noise_percentage_gk1, noise_percentage_gkr, noise_percentage_gcal, noise_percentage_v, noise_percentage_cai), np.array([v_peak_time]))
np.savetxt('averaged_observed_v_{0}_{1}_{2}_{3}_12mm_strip_{4}_{5}_{6}_{7}_{8}_{9}.txt'.format(gnaf,gk1f,gkrf,gcalf, \
     noise_percentage_gna, noise_percentage_gk1, noise_percentage_gkr, noise_percentage_gcal, noise_percentage_v, noise_percentage_cai), values_v)
np.savetxt('averaged_observed_cai_{0}_{1}_{2}_{3}_12mm_strip_{4}_{5}_{6}_{7}_{8}_{9}.txt'.format(gnaf,gk1f,gkrf,gcalf, \
     noise_percentage_gna, noise_percentage_gk1, noise_percentage_gkr, noise_percentage_gcal, noise_percentage_v, noise_percentage_cai), values_c)
np.savetxt('averaged_observed_times_{0}_{1}_{2}_{3}_12mm_strip_{4}_{5}_{6}_{7}_{8}_{9}.txt'.format(gnaf,gk1f,gkrf,gcalf, \
     noise_percentage_gna, noise_percentage_gk1, noise_percentage_gkr, noise_percentage_gcal, noise_percentage_v, noise_percentage_cai), times)

print "Success!"