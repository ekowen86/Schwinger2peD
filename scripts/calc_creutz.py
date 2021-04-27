import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.optimize as opt
import cmath
import sys
import warnings
from statistics import jackknifeMean
from statistics import jackknifeCreutz
from statistics import autocorrTime

warnings.simplefilter(action='ignore', category=FutureWarning)

L = int(sys.argv[1]) # lattice size
print("L: %d" % (L))

Lz = int(sys.argv[2])
print("Lz: %d" % (Lz))

beta = float(sys.argv[3])
print("beta: %f" % (beta))

eps3 = float(sys.argv[4])
print("eps3: %f" % (eps3))

m_fermion = float(sys.argv[5])
print("m_fermion: %f" % (m_fermion))

m_sign = "p"
if m_fermion < 0.0:
	m_sign = "m"

if (Lz == 1):
    id = "%d_%d_%s%d" % (L, round(beta * 1000), m_sign, round(abs(m_fermion) * 1000))
    path = "../jobs/2D/%s" % id
else:
    id = "%d_%d_%d_%d_%s%d" % (L, Lz, round(beta * 1000), round(eps3 * 1000), \
		m_sign, round(abs(m_fermion) * 1000))
    path = "../jobs/3D/%s" % id
print("id: %s" % (id))

first_id = int(sys.argv[6])
print("first_id: %d" % (first_id))

last_id = int(sys.argv[7])
print("last_id: %d" % (last_id))

id_inc = int(sys.argv[8])
print("id_inc: %d" % (id_inc))

t_max = float(sys.argv[9])
print("t_max: %f" % (t_max))

dt = float(sys.argv[10])
print("dt: %f" % (dt))

# r_min = int(sys.argv[11]) - 2
# print("r_min: %d" % (r_min + 2))
#
# r_max = int(sys.argv[12]) - 2
# print("r_max: %d" % (r_max + 2))

# maximum wilson loop size
loopMax = L / 2

# number of wilson loop values
nLoops = loopMax * 3 - 2

# number of trajectories
n_traj = (last_id - first_id) / id_inc + 1
print("n_traj = %d" % (n_traj))

# number of wilson flow values
n_wf = int(t_max / dt + 1)

# first index is wilson flow time
# second index is wilson loop size
# third index is trajectory id
wilsonLoops = np.zeros((n_wf, nLoops, n_traj))


def readWilsonLoops(file, n):
	lines = file.readlines()
	for (l, line) in enumerate(lines):
		tokens = line.split()
		for (t, token) in enumerate(tokens):
			if (t == 0):
				continue
			wilsonLoops[l, t - 1, n] = float(token)
	return


print("\nReading data files...")
nWilsonLoops = 0; # number of field strength values
for i in range(first_id, last_id + 1, id_inc):
	file = open("%s/wf/wilson_loops.%d" % (path, i), "r")
	readWilsonLoops(file, nWilsonLoops)
	nWilsonLoops += 1
	file.close()

# write wilson loop statistics to data file

file = open("%s/wilson_loops.dat" % (path), "w")
for i in range(0, n_wf):
	file.write("%.3f" % (i * dt))
	for w in range(0, nLoops):
		value = jackknifeMean(wilsonLoops[i,w,:])
		tau = autocorrTime(wilsonLoops[i,w,:])
		file.write(" %.12e %.12e %.12e" % (value[0], value[1], tau))
	file.write("\n")

file.close();

print("\nCalculating Creutz Ratios (smeared/raw)...")

# wilson loops are in this order
# (1,1), (1,2)
# (2,1), (2,2), (2,3)
#        (3,2), (3,3), (3,4) ...

# write creutz ratios to file

file = open("%s/creutz.dat" % (path), "w")
for i in range(0, n_wf):
	file.write("%.3f" % (i * dt))
	for r in range(0, loopMax - 1):

		# measure smeared (with wilson flow)
		w00 = wilsonLoops[i, r * 3]
		w01 = wilsonLoops[i, r * 3 + 1]
		w10 = wilsonLoops[i, r * 3 + 2]
		w11 = wilsonLoops[i, r * 3 + 3]
		value = jackknifeCreutz(w00, w01, w10, w11)

		file.write(" %.12e %.12e" % (value[0], value[1]))
	file.write("\n")

file.close();
