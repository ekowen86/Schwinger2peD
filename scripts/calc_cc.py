import numpy as np
import cmath
import sys
import warnings
from statistics import jackknifeMean

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

first_traj = int(sys.argv[6])
print("first_traj: %d" % (first_traj))

m_sign = "p"
if m_fermion < 0.0:
	m_sign = "m"

if (Lz == 1):
    id = "2D/%d_%d_%s%d" % (L, round(beta * 1000), m_sign, round(abs(m_fermion) * 1000))
else:
    id = "3D/%d_%d_%d_%d_%s%d" % (L, Lz, round(beta * 1000), round(eps3 * 1000), m_sign, round(abs(m_fermion) * 1000))
print("id: %s" % (id))


def parse_data_file(file):
	lines = file.readlines()

	first_l = 0
	for l, line in enumerate(lines):
		tokens = line.split()
		traj = int(tokens[0])
		if traj < first_traj:
			# skip lines before first id
			continue
		elif traj == first_traj:
			# create the array
			first_l = l;
			a = np.empty(len(lines) - l)

		# populate the array
		for t, token in enumerate(tokens):
			if (t == 0):
				continue # skip trajectory id
			a[l - first_l] = float(token)
	return a


cc_file = open("../jobs/%s/cc.dat" % (id), "r")
cc = jackknifeMean(parse_data_file(cc_file))
cc_file.close()

print("cc = %.12f (%.12f)" % (cc[0], cc[1]))

result_file = open("../jobs/cc.dat", "a")
result_file.write("%d %d %.3f %.3f %.3f %.12e %.12e\n" % (L, Lz, beta, eps3, m_fermion, cc[0], cc[1]))
result_file.close()
