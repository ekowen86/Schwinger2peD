import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import cmath
import sys
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)


L = int(sys.argv[1])
print("L: %d" % (L))

Lz = int(sys.argv[2])
print("Lz: %d" % (Lz))

beta = float(sys.argv[3])
print("beta: %f" % (beta))

eps3 = float(sys.argv[4])
print("eps3: %f" % (eps3))

if (Lz == 1):
	id = "%d_%d" % (L, round(beta * 1000))
else:
	id = "%d_%d_%d_%d" % (L, Lz, round(beta * 1000), round(eps3 * 1000))
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
	begin = 0
	end = len(lines)
	for l, line in enumerate(lines):
		tokens = line.split()
		for t, token in enumerate(tokens):
			if token == '*b':
				begin = l
				continue
			elif token == '*e':
				end = l + 1
				continue
			a[t,l] = float(token)
	return (a, a[:,begin:end])


# read chiral condensate data
if (Lz == 1):
	cc_file = open("../jobs/2D/cc_%s.dat" % (id), "r")
else:
	cc_file = open("../jobs/3D/cc_%s.dat" % (id), "r")
(cc_all, cc) = parse_data_file(cc_file)

# read pion mass data
if (Lz == 1):
	mass_file = open("../jobs/2D/m_pi_%s.dat" % (id), "r")
else:
	mass_file = open("../jobs/3D/m_pi_%s.dat" % (id), "r")
(m_pi_all, m_pi) = parse_data_file(mass_file)


combined_m = np.append(cc[0,:], m_pi[0,:])
combined_data = np.append(cc[1,:], m_pi[1,:])
combined_err = np.append(cc[2,:], m_pi[2,:])

# fitting functions from Smilga
def cc_fit(m, mc, g, p0):
	return p0 + 0.388 * np.sign(m + mc) * np.absolute(m + mc)**(1.0/3.0) * g**(2.0/3.0)

def m_fit(m, mc, g):
	return 2.008 * np.absolute(m + mc)**(2.0/3.0) * g**(1.0/3.0)

def combined_fit(m, mc, g, p0):
	m1 = m[:np.shape(cc)[1]]
	m2 = m[np.shape(cc)[1]:]
	return np.append(cc_fit(m1, mc, g, p0), m_fit(m2, mc, g))


# fit pion mass to get critical mass and effective coupling
# popt, pcov = opt.curve_fit(m_fit, m_pi[0], m_pi[1], [0.13, 0.5], m_pi[2])
# popt, pcov = opt.curve_fit(cc_fit, cc[0], cc[1], [0.13, 0.5, 0.8], cc[2])
popt, pcov = opt.curve_fit(combined_fit, combined_m, combined_data, [0.32, 1.0, 0.8], combined_err)
mc = (popt[0], np.sqrt(pcov[0][0]))
g = (popt[1], np.sqrt(pcov[1][1]))
p0 = (popt[2], np.sqrt(pcov[2][2]))
# p0 = (0.795, 0.0)

print("mc = %.12f (%.12f)" % (mc[0], mc[1]))
print("g  = %.12f (%.12f)" % (g[0], g[1]))
print("p0 = %.12f (%.12f)" % (p0[0], p0[1]))


m_min = -0.6
m_max = 0.2

# create best fit curve data
m_A = np.linspace(m_min, m_max, 1000)
cc_A = np.empty(len(m_A))
M_A = np.empty(len(m_A))
for i in range(0, len(m_A)):
	M_A[i] = m_fit(m_A[i], mc[0], g[0])
	cc_A[i] = cc_fit(m_A[i], mc[0], g[0], p0[0])


plt.rcParams.update({
	"text.usetex": False,
	"font.family": "sans-serif"})


# plot M vs m
plt.figure()
plt.xlim(m_min, m_max)
plt.errorbar(m_pi_all[0], m_pi_all[1], yerr=m_pi_all[2], color="black", marker='o', ms=5, mew=0.5, mfc='none', linestyle='none', linewidth=0.5, capsize=2.5, capthick=0.5)
plt.plot(m_A, M_A, color="black", linewidth=0.5, linestyle='dashed', label=(r"$M_\pi = 2.008 (m_0 + m_c)^{2/3} g^{1/3}$" + "\n $m_c = %.3f \pm %.3f$\n $g = %.3f \pm %.3f$" % (mc[0], mc[1], g[0], g[1])))
plt.legend(loc='upper left', fontsize=10, handlelength=2.5, frameon=False)
plt.xlabel("$m_0$")
plt.ylabel("$M_{\pi}$")
if (Lz == 1):
	plt.savefig("../jobs/2D/m_pi_%s.pdf" % (id))
else:
	plt.savefig("../jobs/3D/m_pi_%s.pdf" % (id))
plt.close()

# plot cc vs m
plt.figure()
plt.xlim(m_min, m_max)
plt.errorbar(cc_all[0], cc_all[1], yerr=cc_all[2], color="black", marker='o', ms=5, mew=0.5, mfc='none', linestyle='none', linewidth=0.5, capsize=2.5, capthick=0.5)
plt.plot(m_A, cc_A, color="black", linewidth=0.5, linestyle='dashed', label=(r"$\langle \bar{\psi} \psi \rangle = p_c + 0.388 (m_0 + m_c)^{1/3} g^{2/3}$" + "\n $m_c = %.3f \pm %.3f$\n $g = %.3f \pm %.3f$\n $p_c = %.3f \pm %.3f$" % (mc[0], mc[1], g[0], g[1], p0[0], p0[1])))
plt.legend(loc='center right', fontsize=10, handlelength=2.5, frameon=False)
plt.xlabel(r"$m_0$")
plt.ylabel(r"$\langle \bar{\psi} \psi \rangle$")
if (Lz == 1):
	plt.savefig("../jobs/2D/cc_%s.pdf" % (id))
else:
	plt.savefig("../jobs/3D/cc_%s.pdf" % (id))
plt.close()
