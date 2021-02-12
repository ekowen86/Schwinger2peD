import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.optimize as opt
import cmath
import scipy.special as special
import sys
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)


L = int(sys.argv[1]) # lattice size
print("L: %d" % (L))
Lz = int(sys.argv[2]) # lattice size
print("Lz: %d" % (Lz))
betaMin = float(sys.argv[3])
print("betaMin: %f" % (betaMin))
betaMax = float(sys.argv[4])
print("betaMax: %f" % (betaMax))
betaInc = float(sys.argv[5])
print("betaInc: %f" % (betaInc))
eps3 = float(sys.argv[6])
print("eps3: %f" % (eps3))

first_id = 200 # first configuration id

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
			a = np.empty((len(tokens), len(lines) - l))

		for (t, token) in enumerate(tokens):
			a[t, l - first_l] = float(tokens[t])
	return a

n = int(round((betaMax - betaMin) / betaInc)) + 1
beta = np.empty(n)
plaq = np.empty((n, 2))
w12 = np.empty((n, 2))
w21 = np.empty((n, 2))
w22 = np.empty((n, 2))

for i in range(0, n):

	beta[i] = betaMin + i * betaInc

	id = "%d_%d_%d_%d" % (L, Lz, round(beta[i] * 1000), round(eps3 * 1000))
# 	print("id: %s" % (id))

	plaq_file = open("plaq/plaq_%s.dat" % (id), "r")
	plaq_data = parse_data_file(plaq_file)
	plaq_file.close()

	plaq[i] = jackknife_mean(plaq_data[1])
	print("plaq(%.2f) = %.12f (%.12f)" % (beta[i], plaq[i,0], plaq[i,1]))

	w12[i] = jackknife_mean(plaq_data[3])
	print("w12(%.2f)  = %.12f (%.12f)" % (beta[i], w12[i,0], w12[i,1]))

	w21[i] = jackknife_mean(plaq_data[4])
	print("w21(%.2f)  = %.12f (%.12f)" % (beta[i], w21[i,0], w21[i,1]))

	w22[i] = jackknife_mean(plaq_data[5])
	print("w22(%.2f)  = %.12f (%.12f)" % (beta[i], w22[i,0], w22[i,1]))
	print("")

# calculate best fit curves
beta_A = np.linspace(0, betaMax + betaInc)
plaq_A = np.zeros(len(beta_A))
for i in range(0, len(beta_A)):
	plaq_A[i] = special.iv(1,beta_A[i]) / special.iv(0,beta_A[i])

# matplotlib.rc('text', usetex = False)
# matplotlib.rc('font', **{'family' : "sans-serif"})
# params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
# plt.rcParams.update(params)

# plot plaquette vs beta
# plt.figure()
# plt.xlim(2.0, betaMax + betaInc)
# plt.ylim(0.7, 0.9)
# plt.errorbar(beta, plaq[:,0], yerr=plaq[:,1], color="black", marker='o', ms=5, mew=0.5, mfc='none', linestyle='none', linewidth=0.5, capsize=2.5, capthick=0.5)
# plt.plot(beta_A, plaq_A, color="black", linewidth=0.5, linestyle='dashed')
# plt.xlabel(r"$\beta$")
# plt.ylabel(r"Plaq$(\beta)$")
# plt.savefig("plaq_%d_%d_%d.pdf" % (L, Lz, eps3 * 1000))
# plt.close()
