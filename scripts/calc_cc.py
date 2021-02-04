import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import cmath
import sys
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)


L = int(sys.argv[1]) # lattice size
print("L: %d" % (L))
beta = float(sys.argv[2])
print("beta: %f" % (beta))
m_fermion = float(sys.argv[3])
print("m_fermion: %f" % (m_fermion))
first_id = 200 # first configuration id

m_sign = "p"
if m_fermion < 0.0:
    m_sign = "m"

id = "%d_%d_%s%d" % (L, beta * 1000, m_sign, abs(m_fermion) * 1000)
print("id: %s" % (id))

def jackknife_mean(a):
	n = len(a)
	f_n = float(n)
	a_bar = np.mean(a)
	d_a = 0.0

	for i in range(0, n):

		# copy the array, deleting the current value
		# and add to the error
		d_a += (np.mean(np.delete(a, i)) - a_bar)**2.0

	d_a = np.sqrt((f_n - 1) / f_n * d_a)
	return (a_bar, d_a)


def parse_data_file(file):
	lines = file.readlines()

	first_l = 0
	for l, line in enumerate(lines):
		tokens = line.split()
		id = int(tokens[0])
		if id < first_id:
			# skip lines before first id
			continue
		elif id == first_id:
			# create the array
			first_l = l;
			a = np.empty(len(lines) - l)

		# populate the array
		for t, token in enumerate(tokens):
			if (t == 0):
				continue # skip trajectory id
			a[l - first_l] = float(token)
	return a


cc_file = open("../jobs/2D/%s/cc.dat" % (id), "r")
cc = parse_data_file(cc_file)
cc_bar, d_cc = jackknife_mean(cc)
cc_file.close()

print("cc = %.12f (%.12f)" % (cc_bar, d_cc))

result_file = open("../jobs/2D/cc_%d_%d.dat" % (L, beta * 1000), "a")
result_file.write("%.3f %.12e %.12e\n" % (m_fermion, cc_bar, d_cc))
result_file.close()
