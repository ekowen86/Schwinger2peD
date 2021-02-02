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
R_half = int(L / 2) # half lattice size
R_min = int(sys.argv[4]) # min R value for fit
print("R_min: %d" % (R_min))
R_max = R_half # max R value for fit
first_id = 100 # first configuration id

m_sign = "p"
if m_fermion < 0.0:
    m_sign = "m"

id = "%d_%d_%s%d" % (L, beta * 1000, m_sign, abs(m_fermion) * 1000)
print("id: %s" % (id))

def corr_fit(r, m, z):
	# a = 0.0
    #
	# # this includes periodic finite size effects from 5 lattices away (only affects z, not m)
	# for n in range(-5, 6):
	# 	a += z * np.exp(-m * abs(r - L * n))
    #
	# return a

	a = z * np.exp(-m * r)
	b = z * np.exp(-m * (L - r))
	return a + b


def pion_mass(corr):
	n = float(corr.shape[1])
	r = np.empty(R_half + 1)
	corr_bar = np.empty(R_half + 1)
	d_corr = np.empty(R_half + 1)

	for i in range(0, R_half + 1):
		r[i] = i;
		corr_bar[i] = np.mean(corr[i]);
		d_corr[i] = np.std(corr[i]) / np.sqrt(n);

	popt, pcov = opt.curve_fit(corr_fit, r[R_min-1:R_max], corr_bar[R_min-1:R_max], [1.0, 1.0], sigma=d_corr[R_min-1:R_max])
	m = popt[0]
	d_m = np.sqrt(pcov[0][0])
	z = popt[1]
	d_z = np.sqrt(pcov[1][1])
	return m, d_m, z, d_z


def jackknife_pion_mass(corr):
	n = corr.shape[1]
	f_n = float(n) # number of data subsets

	m_bar, _, z_bar, _ = pion_mass(corr)
	d_m = 0.0
	d_z = 0.0

	for i in range(0, n):
		# copy the array and delete the current id
		# calculate the pion mass and add to the error
		m_d, _, z_d, _ = pion_mass(np.delete(corr, i, axis=1))
		d_m += (m_d - m_bar)**2.0
		d_z += (z_d - z_bar)**2.0

	d_m = np.sqrt((f_n - 1) / f_n * d_m)
	d_z = np.sqrt((f_n - 1) / f_n * d_z)
	return m_bar, d_m, z_bar, d_z


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
			# 0 axis is distance, 1 axis is trajectory
			a = np.empty((R_half + 1, len(lines) - l))

		# populate the array
		for t, token in enumerate(tokens):
			if (t == 0):
				continue # skip trajectory id
			a[t-1, l - first_l] = float(token)
	return a


pion_file = open("./jobs/2D/%s/pion_corr/pion_corr.dat" % (id), "r")
C_pi = parse_data_file(pion_file)
m_bar, d_m, z_bar, d_z = jackknife_pion_mass(C_pi)
pion_file.close()

print("m = %.12f (%.12f)" % (m_bar, d_m))
print("z = %.12f (%.12f)" % (z_bar, d_z))

mass_file = open("../jobs/2D/pi_mass_%d_%d.dat" % (L, beta * 1000), "a")
mass_file.write("%.3f %.12e %.12e\n" % (m_fermion, m_bar, d_m))
mass_file.close()

n = R_half + 1
R = np.empty(n)
corr_bar = np.empty(n)
d_corr = np.empty(n)

for i in range(0, n):
	R[i] = i;
	corr_bar[i] = np.mean(C_pi[i]);
	d_corr[i] = np.std(C_pi[i]) / np.sqrt(n);

# calculate best fit curves
R_A = np.linspace(0, R_half, 1000)
corr_A = np.zeros(len(R_A))
for r in range(0, len(R_A)):
	corr_A[r] = corr_fit(R_A[r], m_bar, z_bar)

plt.rcParams.update({
	"text.usetex": False,
	"font.family": "sans-serif"})

# plot C vs r
plt.figure()
plt.xlim(0, R_half)
plt.yscale("log")
# plt.ylim(0.0, chi[1] + 0.1)
plt.errorbar(R, corr_bar, yerr=d_corr, color="blue", marker='o', ms=5, mew=0.5, mfc='none', linestyle='none', linewidth=0.5, capsize=2.5, capthick=0.5)
plt.plot(R_A, corr_A, color="blue", linewidth=0.5)
plt.xlabel(r"$t$")
plt.ylabel(r"$\langle C_{\pi}(0) C_{\pi}(t) \rangle $")
plt.savefig("../jobs/2D/%s/plots/pi_corr.pdf" % (id))
plt.close()
