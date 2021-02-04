import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import cmath
import sys
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)


L = int(sys.argv[1])
print("L: %d" % (L))
beta = float(sys.argv[2])
print("beta: %f" % (beta))

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

	a = np.empty((5,len(lines)))
	for l, line in enumerate(lines):
		tokens = line.split()
		for t, token in enumerate(tokens):
			a[t,l] = float(token)
	return a


# fitting functions from Smilga
def m_fit(m, mc, g):
	return 2.008 * np.absolute(m + mc)**(2.0/3.0) * g**(1.0/3.0)


# read pion mass data
mass_file = open("../jobs/2D/m_pi_%s.dat" % (id), "r")
m_pi = parse_data_file(mass_file)

# fit pion mass to get critical mass and effective coupling
popt, pcov = opt.curve_fit(m_fit, m_pi[0], m_pi[1], [0.35, 1.0 / np.sqrt(beta)], sigma=m_pi[2])
mc = (popt[0], np.sqrt(pcov[0][0]))
g = (popt[1], np.sqrt(pcov[1][1]))

print("mc = %.12f (%.12f)" % (mc[0], mc[1]))
print("g  = %.12f (%.12f)" % (g[0], g[1]))


# read chiral condensate data
cc_file = open("../jobs/2D/cc_%s.dat" % (id), "r")
cc = parse_data_file(cc_file)

def cc_fit(m, p0):
	return p0 + 0.388 * np.sign(m + mc[0]) * np.abs(m + mc[0])**(1.0/3.0) * g[0]**(2.0/3.0)

# fit pion mass to get critical mass and effective coupling
popt, pcov = opt.curve_fit(cc_fit, cc[0], cc[1], [1.0], sigma=cc[2])
p0 = (popt[0], np.sqrt(pcov[0][0]))

print("p0 = %.12f (%.12f)" % (p0[0], p0[1]))

m_min = -0.4
m_max = 0.2

# create best fit curve data
m_A = np.linspace(m_min, m_max, 1000)
cc_A = np.empty(len(m_A))
M_A = np.empty(len(m_A))
for i in range(0, len(m_A)):
	M_A[i] = m_fit(m_A[i], mc[0], g[0])
	cc_A[i] = cc_fit(m_A[i], p0[0])


plt.rcParams.update({
	"text.usetex": False,
	"font.family": "sans-serif"})


# plot M vs m
plt.figure()
plt.xlim(m_min, m_max)
plt.errorbar(m_pi[0], m_pi[1], yerr=m_pi[2], color="black", marker='o', ms=5, mew=0.5, mfc='none', linestyle='none', linewidth=0.5, capsize=2.5, capthick=0.5)
plt.plot(m_A, M_A, color="black", linewidth=0.5, linestyle='dashed', label=(r"$M_\pi = 2.008 (m_0 + m_c)^{2/3} g^{1/3}$" + "\n $m_c = %.3f \pm %.3f$\n $g = %.3f \pm %.3f$" % (mc[0], mc[1], g[0], g[1])))
plt.legend(loc='upper left', fontsize=10, handlelength=2.5, frameon=False)
plt.xlabel("$m_0$")
plt.ylabel("$M_{\pi}$")
plt.savefig("../jobs/2D/m_pi_%s.pdf" % (id))
plt.close()

# plot cc vs m
plt.figure()
plt.xlim(m_min, m_max)
plt.errorbar(cc[0], cc[1], yerr=cc[2], color="black", marker='o', ms=5, mew=0.5, mfc='none', linestyle='none', linewidth=0.5, capsize=2.5, capthick=0.5)
plt.plot(m_A, cc_A, color="black", linewidth=0.5, linestyle='dashed', label=(r"$\langle \bar{\psi} \psi \rangle = p_c + 0.388 (m_0 + m_c)^{1/3} g^{2/3}$" + "\n $m_c = %.3f \pm %.3f$\n $g = %.3f \pm %.3f$\n $p_c = %.3f \pm %.3f$" % (mc[0], mc[1], g[0], g[1], p0[0], p0[1])))
plt.legend(loc='center right', fontsize=10, handlelength=2.5, frameon=False)
plt.xlabel(r"$m_0$")
plt.ylabel(r"$\langle \bar{\psi} \psi \rangle$")
plt.savefig("../jobs/2D/cc_%s.pdf" % (id))
plt.close()
