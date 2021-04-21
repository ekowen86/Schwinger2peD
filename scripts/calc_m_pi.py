import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.optimize as opt
import cmath
import sys
import warnings

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

R_half = int(L / 2) # half lattice size
tmin = int(sys.argv[6]) # min t value for fit
print("tmin: %d" % (tmin))
tmax = R_half # max t value for fit

first_traj = int(sys.argv[7])
print("first_traj: %d" % (first_traj))

m_sign = "p"
if m_fermion < 0.0:
	m_sign = "m"

if (Lz == 1):
    id = "2D/%d_%d_%s%d" % (L, round(beta * 1000), m_sign, round(abs(m_fermion) * 1000))
else:
    id = "3D/%d_%d_%d_%d_%s%d" % (L, Lz, round(beta * 1000), round(eps3 * 1000), m_sign, round(abs(m_fermion) * 1000))
print("id: %s" % (id))

def corr_fit(t, m, z):
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
	z = popt[1]

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
		# copy the array and delete the current trajectory
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
		traj = int(tokens[0])
		if traj < first_traj:
			# skip lines before first trajectory
			continue
		elif traj == first_traj:
			# create the array
			first_l = l;
			# 0 axis is distance, 1 axis is trajectory
			a = np.empty((R_half + 1, len(lines) - l))

		# populate the array
		for t, token in enumerate(tokens):
			if (t == 0):
				continue # skip trajectory index
			a[t-1, l - first_l] = float(token)
	return a


# parse the data file
pion_file = open("../jobs/%s/pion_corr.dat" % (id), "r")
C_pi = parse_data_file(pion_file)

m_bar, d_m, z_bar, d_z = jackknife_pion_mass(C_pi, tmin, tmax)
m = (m_bar, d_m)
z = (z_bar, d_z)
pion_file.close()

# print("tmin = %d" % (tmin_best))
# print("tmax = %d" % (tmax_best - 1))
print("m = %.12f (%.12f)" % (m[0], m[1]))
print("z = %.12f (%.12f)" % (z[0], z[1]))

# if (Lz == 1):
# 	mass_file = open("../jobs/2D/m_pi_%d_%d.dat" % (L, beta * 1000), "a")
# else:
# 	mass_file = open("../jobs/3D/m_pi_%d_%d_%d_%d.dat" % (L, Lz, round(beta * 1000), round(eps3 * 1000)), "a")
# mass_file.write("%.3f %.12e %.12e %.12e %.12e\n" % (m_fermion, m_bar, d_m, z_bar, d_z))
# mass_file.close()

result_file = open("m_pi.dat", "a")
result_file.write("%d %d %.3f %.3f %.3f %.12e %.12e %.12e %.12e\n" % (L, Lz, beta, eps3, m_fermion, m[0], m[1], z[0], z[1]))
result_file.close()

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
	corr_A[r] = corr_fit(R_A[r], m[0], z[0])

plt.rcParams.update({
	"text.usetex": False,
	"font.family": "sans-serif"})

# plot C vs r
plt.figure()
plt.xlim(0, R_half + 1)
plt.yscale("log")
plt.errorbar(R, corr_bar, yerr=d_corr, color="blue", marker='o', ms=5, mew=0.5, mfc='none', linestyle='none', linewidth=0.5, capsize=2.5, capthick=0.5)
plt.plot(R_A, corr_A, color="blue", linewidth=0.5)
plt.xlabel(r"$t$")
plt.ylabel(r"$\langle C_{\pi}(0) C_{\pi}(t) \rangle $")
plt.savefig("../jobs/%s/pi_corr.pdf" % (id))
plt.close()
