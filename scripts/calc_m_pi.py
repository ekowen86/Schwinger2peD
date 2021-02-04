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
tmin = int(sys.argv[4]) # min t value for fit
print("tmin: %d" % (tmin))
tmax = R_half # max t value for fit
first_id = 200 # first configuration id

m_sign = "p"
if m_fermion < 0.0:
	m_sign = "m"

id = "%d_%d_%s%d" % (L, beta * 1000, m_sign, abs(m_fermion) * 1000)
print("id: %s" % (id))

def corr_fit(t, m, z):
	# a = 0.0
	#
	# # this includes periodic finite size effects from 5 lattices away (only affects z, not m)
	# for n in range(-5, 6):
	#	a += z * np.exp(-m * abs(t - L * n))
	#
	# return a

	a = z * np.exp(-m * t)
	b = z * np.exp(-m * (L - t))
	return a + b


def pion_mass(corr, tmin, tmax):
	n = float(corr.shape[1])
	r = np.empty(R_half + 1)
	corr_bar = np.empty(R_half + 1)
	d_corr = np.empty(R_half + 1)

	for i in range(0, R_half + 1):
		r[i] = i;
		corr_bar[i] = np.mean(corr[i]);
		d_corr[i] = np.std(corr[i]) / np.sqrt(n);

	popt, pcov = opt.curve_fit(corr_fit, r[tmin:tmax], corr_bar[tmin:tmax], [1.0, 1.0], sigma=d_corr[tmin:tmax])
	m = popt[0]
	# d_m = np.sqrt(pcov[0][0])
	z = popt[1]
	# d_z = np.sqrt(pcov[1][1])

	err = 0.0
	for t in range(tmin, tmax):
		err += (corr_fit(t, m, z) - corr_bar[t])**2.0

	return m, z, np.sqrt(err) / float(tmax - tmin)


def jackknife_pion_mass(corr, tmin, tmax):
	n = corr.shape[1]
	f_n = float(n) # number of data subsets

	m_bar, z_bar, _ = pion_mass(corr, tmin, tmax)
	d_m = 0.0
	d_z = 0.0

	for i in range(0, n):
		# copy the array and delete the current id
		# calculate the pion mass and add to the error
		m_del, z_del, _ = pion_mass(np.delete(corr, i, axis=1), tmin, tmax)
		d_m += (m_del - m_bar)**2.0
		d_z += (z_del - z_bar)**2.0

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


# parse the data file
pion_file = open("../jobs/2D/%s/pion_corr.dat" % (id), "r")
C_pi = parse_data_file(pion_file)

# find tmin and tmax which minimize the determinant of the covariance matrix
# tmin_best = 0
# tmax_best = R_half
# err_best = np.inf
# for tmin in range(0, R_half / 2):
# 	for tmax in range(tmin + R_half / 2, R_half):
# 		_, _, err = pion_mass(C_pi, tmin, tmax)
# 		print("%d, %d, %.12e" % (tmin, tmax - 1, err))
# 		if err < err_best:
# 			tmin_best = tmin
# 			tmax_best = tmax
# 			err_best = err


m_bar, d_m, z_bar, d_z = jackknife_pion_mass(C_pi, tmin, tmax)
pion_file.close()

# print("tmin = %d" % (tmin_best))
# print("tmax = %d" % (tmax_best - 1))
print("m = %.12f (%.12f)" % (m_bar, d_m))
print("z = %.12f (%.12f)" % (z_bar, d_z))

mass_file = open("../jobs/2D/m_pi_%d_%d.dat" % (L, beta * 1000), "a")
mass_file.write("%.3f %.12e %.12e %.12e %.12e\n" % (m_fermion, m_bar, d_m, z_bar, d_z))
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
R_A = np.linspace(0, R_half + 1, 1000)
corr_A = np.zeros(len(R_A))
for r in range(0, len(R_A)):
	corr_A[r] = corr_fit(R_A[r], m_bar, z_bar)

plt.rcParams.update({
	"text.usetex": False,
	"font.family": "sans-serif"})

# plot C vs r
plt.figure()
plt.xlim(0, R_half + 1)
plt.yscale("log")
# plt.ylim(0.0, chi[1] + 0.1)
plt.errorbar(R, corr_bar, yerr=d_corr, color="blue", marker='o', ms=5, mew=0.5, mfc='none', linestyle='none', linewidth=0.5, capsize=2.5, capthick=0.5)
plt.plot(R_A, corr_A, color="blue", linewidth=0.5)
plt.xlabel(r"$t$")
plt.ylabel(r"$\langle C_{\pi}(0) C_{\pi}(t) \rangle $")
plt.savefig("../jobs/2D/%s/plots/pi_corr.pdf" % (id))
plt.close()
