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
g = 1.0 / np.sqrt(beta)

mc = float(sys.argv[3])
print("mc: %f" % (mc))

id = "%d_%d" % (L, beta * 1000)
print("id: %s" % (id))


def not_comment(line):
	if line[0] == "#":
		return False
	else:
		return True


def parse_data_file(file):
	lines = file.readlines()

	# remove comments
	lines = filter(not_comment, lines)

	a = np.empty((3,len(lines)))
	for l, line in enumerate(lines):
		tokens = line.split()
		for t, token in enumerate(tokens):
			a[t,l] = float(token)
	return a


def cc_fit(m, p0):
	return p0 - 0.388 * np.sign(m + mc) * np.abs(m + mc)**(1.0/3.0) * g**(2.0/3.0)


cc_file = open("../jobs/2D/cc_%s.dat" % (id), "r")
cc = parse_data_file(cc_file)

# calculate best fit for cc vs m_fermion
popt, pcov = opt.curve_fit(cc_fit, cc[0], cc[1], [1.0], sigma=cc[2])
p0 = popt[0]
d_p0 = np.sqrt(pcov[0][0])

print("p0    = %.12f (%.12f)" % (p0, d_p0))

m_A = np.linspace(-mc, 0.5, 1000)
cc_A = np.empty(len(m_A))
for m in range(0, len(m_A)):
	cc_A[m] = cc_fit(m_A[m], p0)

plt.rcParams.update({
	"text.usetex": False,
	"font.family": "sans-serif"})

# plot cc vs m
plt.figure()
plt.xlim(-0.4, 0.2)
plt.errorbar(cc[0], cc[1], yerr=cc[2], color="blue", marker='o', ms=5, mew=0.5, mfc='none', linestyle='none', linewidth=0.5, capsize=2.5, capthick=0.5)
plt.plot(m_A, cc_A, color="blue", linewidth=0.5)
plt.xlabel(r"$m_0$")
plt.ylabel(r"$\langle \bar{\psi} \psi \rangle$")
plt.savefig("../jobs/2D/cc_%s.pdf" % (id))
plt.close()
