import numpy as np
import sys
import warnings
from statistics import autocorrTime

warnings.simplefilter(action='ignore', category=FutureWarning)


L = int(sys.argv[1]) # lattice size
print("L: %d" % (L))

Lz = int(sys.argv[2])
print("Lz: %d" % (Lz))
zCenter = int((Lz - 1) / 2)
print("zCenter: %d" % (zCenter))

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
    id = "3D/%d_%d_%d_%d_%s%d" % (L, Lz, round(beta * 1000), round(abs(eps3) * 1000), m_sign, round(abs(m_fermion) * 1000))
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
		a[l - first_l] = int(tokens[zCenter + 1])
	return a


topo_charge_file = open("../jobs/%s/topo_charge.dat" % (id), "r")
topo_charge = parse_data_file(topo_charge_file)
topo_charge_file.close()

n_frozen = 0
q_prev = topo_charge[0]
for q in topo_charge[1:]:
	if (q == q_prev):
		n_frozen += 1
	q_prev = q
p_frozen = float(n_frozen) / len(topo_charge)

tau_autocorr = autocorrTime(topo_charge)

print("n_frozen = %d" % (n_frozen))
print("p_frozen = %.6f" % (p_frozen))
print("tau = %.12f" % (tau_autocorr))

result_file = open("../jobs/topo_freeze.dat", "a")
result_file.write("%d %d %.3f %.3f %.3f %.6f %.12e\n" % (L, Lz, beta, eps3, m_fermion, p_frozen, tau_autocorr))
result_file.close()
