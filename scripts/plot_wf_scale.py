import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.optimize as opt
import cmath
import sys
import warnings
from statistics import jackknifeMean
from statistics import autocorrTime

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

m_sign = "p"
if m_fermion < 0.0:
	m_sign = "m"

if (Lz == 1):
    id = "%d_%d_%s%d" % (L, round(beta * 1000), m_sign, round(abs(m_fermion) * 1000))
    path = "../jobs/2D/%s" % id
else:
    id = "%d_%d_%d_%d_%s%d" % (L, Lz, round(beta * 1000), round(eps3 * 1000), m_sign, round(abs(m_fermion) * 1000))
    path = "../jobs/3D/%s" % id
print("id: %s" % (id))

first_id = int(sys.argv[6])
print("first_id: %d" % (first_id))

last_id = int(sys.argv[7])
print("last_id: %d" % (last_id))

id_inc = int(sys.argv[8])
print("id_inc: %d" % (id_inc))

t_max = float(sys.argv[9])
print("t_max: %f" % (t_max))

dt = float(sys.argv[10])
print("dt: %f" % (dt))

E0 = float(sys.argv[11])
print("E0: %f" % (E0))


def calc_t0(t, E):
	n = len(t)
	t2E = np.empty(n)
	for i in range(0, n):
		t2E[i] = t[i] * t[i] * np.mean(E[i,:])

	# find the index of the first time that t^2 E goes above E0
	findE0 = np.where(t2E > E0)
	if (len(findE0) == 0):
		return 0.0

	i = findE0[0][0]
	t2E_1 = (t2E[i] - t2E[i-1]) / dt # first derivative
	t0 = (E0 - t2E[i-1]) / t2E_1 + (t[i] - dt) # interpolate
	return t0


def jackknife_t0(t, E):
	n = E.shape[1]
	f_n = float(n)
	t0_bar = calc_t0(t, E)
	d_t0 = 0.0

	for i in range(0, n):

		# copy the array, deleting the current value
		# and add to the error
		E_del = np.delete(E, i, axis=1)
		d_t0 += (calc_t0(t, E_del) - t0_bar)**2.0

	d_t0 = np.sqrt((f_n - 1) / f_n * d_t0)
	return (t0_bar, d_t0)


def calc_w0(t, E):
	n = len(t)
	t2E = np.empty(n)
	for i in range(0, n):
		t2E[i] = t[i] * t[i] * np.mean(E[i,:])

	W = np.zeros(n)
	for i in range(1, n-1):
		W[i] = (t2E[i + 1] - t2E[i - 1]) / (2.0 * dt) * t[i]

	# find the index of the first time that W goes above E0
	findE0 = np.where(W > E0)
	if (len(findE0) == 0):
		return 0.0

	i = findE0[0][0]
	W_1 = (W[i] - W[i-1]) / dt # first derivative
	return (E0 - W[i-1]) / W_1 + (t[i] - dt) # interpolate


def jackknife_w0(t, E):
	n = E.shape[1]
	f_n = float(n)
	w0_bar = calc_w0(t, E)
	d_w0 = 0.0

	for i in range(0, n):

		# copy the array, deleting the current value
		# and add to the error
		E_del = np.delete(E, i, axis=1)
		d_w0 += (calc_w0(t, E_del) - w0_bar)**2.0

	d_w0 = np.sqrt((f_n - 1) / f_n * d_w0)
	return (w0_bar, d_w0)


def calc_W(Em, Ep, t):
	t2Ep = (t + dt)**2.0 * np.mean(Ep)
	t2Em = (t - dt)**2.0 * np.mean(Em)
	return (t2Ep - t2Em) / (2.0 * dt) * t


def jackknife_W(Em, Ep, t):
	n = len(Em)
	f_n = float(n)
	W_bar = calc_W(Em, Ep, t)
	d_W = 0.0

	for i in range(0, n):

		# copy the array, deleting the current value
		# and add to the error
		Em_del = np.delete(Em, i)
		Ep_del = np.delete(Ep, i)
		d_W += (calc_W(Em_del, Ep_del, t) - W_bar)**2.0

	d_W = np.sqrt((f_n - 1) / f_n * d_W)
	return (W_bar, d_W)


# number of trajectories
n_traj = (last_id - first_id) / id_inc + 1
print("n_traj = %d" % (n_traj))

# number of wilson flow values
n_wf = int(t_max / dt + 1)

# first index is wilson flow time
# second index is trajectory id
field_strength = np.zeros((n_wf, n_traj))

def read_field_strength(file, n):
	lines = file.readlines()
	for (l, line) in enumerate(lines):
		values = line.split()
		t = float(values[0])
		field_strength[l, n] = float(values[1])
	return


print("\nReading data files...")
n_field_strength = 0; # number of field strength values
for i in range(first_id, last_id + 1, id_inc):
	file = open("%s/wf/field_strength.%d" % (path, i), "r")
	read_field_strength(file, n_field_strength)
	n_field_strength += 1
	file.close()

T = np.linspace(0, t_max, n_wf)

print("\nCalculating t0 and w0...")

t0 = jackknife_t0(T, field_strength)
w0 = jackknife_w0(T, field_strength)

print("t0 = %.12f (%.12f)" % (t0[0], t0[1]))
print("w0 = %.12f (%.12f)" % (w0[0], w0[1]))

result_file = open("../jobs/wf_scale.dat", "a")
result_file.write("%d %d %.3f %.3f %.3f %.12e %.12e %.12e %.12e\n" % \
	(L, Lz, beta, eps3, m_fermion, t0[0], t0[1], w0[0], w0[1]))
result_file.close()

E = np.empty((n_wf, 2))
t2E = np.empty((n_wf, 2)) # t^2 * E(t)
W = np.zeros((n_wf, 2)) # t * d/dt(t^2 * E(t))

# write field strength statistics to data file

file = open("%s/field_strength.dat" % (path), "w")
for i in range(0, n_wf):
	file.write("%.3f" % (i * dt))
	value = jackknifeMean(field_strength[i,:])
	tau = autocorrTime(field_strength[i,:])

	t = T[i]
	t2E[i,0] = t * t * value[0]
	t2E[i,1] = t * t * value[1]
	if (i != 0) and (i != n_wf - 1):
		W[i] = jackknife_W(field_strength[i-1,:], field_strength[i+1,:], t)

	file.write(" %.12f %.12f %.12f %.12f %.12f %.12f %.12f\n" % \
		(value[0], value[1], t2E[i,0], t2E[i,1], W[i,0], W[i,1], tau))

file.close();


for i in range(0, n_wf):
	t = T[i]
	E[i] = jackknifeMean(field_strength[i,:])
	t2E[i] = t * t * E[i]

for i in range(1, n_wf - 1):
	t = T[i]
	W[i] = jackknife_W(field_strength[i-1,:], field_strength[i+1,:], t)


print("\nPlotting...")

plt.rcParams.update({
    "text.usetex": False,
    "font.family": "sans-serif"})

# plot t^2 * <E>
plt.figure()
if (t0[0] != 0.0):
	plt.xlim(0, t0[0] * 1.2)
else:
	plt.xlim(0, n_wf * 0.02)
plt.ylim(0.0, E0 * 1.5)
plt.errorbar(T, t2E[:,0], yerr=t2E[:,1], color="blue", marker='o', ms=5, mew=0.5, mfc='none', linestyle='none', linewidth=0.5, capsize=2.5, capthick=0.5)
if (t0[0] != 0.0):
	plt.plot([0, t0[0]], [E0, E0], color='k', linestyle='--', linewidth=0.5)
	plt.plot([t0[0], t0[0]], [0.0, E0], color='k', linestyle='--', linewidth=0.5)
	plt.annotate("$t_0 = %.3f \pm %.3f$" % (t0[0],t0[1]), (t0[0],E0), xytext=(t0[0],E0 * 1.05), ha='right')
plt.xlabel("$t$")
plt.ylabel("$t^2 \\langle E(t) \\rangle$")
plt.savefig("%s/t0.pdf" % (path))
plt.close()

# plot t * d/dt(t^2 * <E>) vs t
plt.figure()
if (w0[0] != 0.0):
	plt.xlim(0, w0[0] * 1.2)
else:
	plt.xlim(0, n_wf * 0.02)
plt.ylim(0.0, E0 * 1.5)
plt.errorbar(T, W[:,0], yerr=W[:,1], color="blue", marker='o', ms=5, mew=0.5, mfc='none', linestyle='none', linewidth=0.5, capsize=2.5, capthick=0.5)
if (w0[0] != 0.0):
	plt.plot([0, w0[0]], [E0, E0], color='k', linestyle='--', linewidth=0.5)
	plt.plot([w0[0], w0[0]], [0.0, E0], color='k', linestyle='--', linewidth=0.5)
	plt.annotate("$w_0^2 = %.3f \pm %.3f$" % (w0[0],w0[1]), (w0[0],E0), xytext=(w0[0],E0 * 1.05), ha='right')
plt.xlabel("$t$")
plt.ylabel("$d/dt \\left( t^2 \\langle E(t) \\rangle \\right) t$")
plt.savefig("%s/w0.pdf" % (path))
plt.close()
