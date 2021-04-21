import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.optimize as opt
import cmath
import sys
import warnings
from statistics import jackknifeMean
from statistics import jackknifeCreutz

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

# maximum wilson loop size
loopMax = L / 2

# number of wilson loop values
nLoops = loopMax * 3 - 2

# number of trajectories
n_traj = (last_id - first_id) / id_inc + 1
print("n_traj = %d" % (n_traj))

# number of wilson flow values
n_wf = int(t_max / dt + 1)

# first index is wilson flow time
# second index is wilson loop size
# third index is trajectory id
wilsonLoops = np.zeros((n_wf, nLoops, n_traj))


def readWilsonLoops(file, n):
	lines = file.readlines()
	for (l, line) in enumerate(lines):
		tokens = line.split()
		for (t, token) in enumerate(tokens):
			if (t == 0):
				continue
			wilsonLoops[l, t - 1, n] = float(token)
	return


print("\nReading data files...")
nWilsonLoops = 0; # number of field strength values
for i in range(first_id, last_id + 1, id_inc):
	file = open("%s/wf/wilson_loops.%d" % (path, i), "r")
	readWilsonLoops(file, nWilsonLoops)
	nWilsonLoops += 1
	file.close()


print("\nCalculating Creutz Ratios...")
R2inv = np.empty(loopMax - 2)
chiRaw = np.empty((len(R2inv), 2))
chiSmeared = np.empty((len(R2inv), 2))

# wilson loops are in this order
# (1,1), (1,2)
# (2,1), (2,2), (2,3)
#        (3,2), (3,3), (3,4) ...

for i in range(0, len(R2inv)):
	R = float(i + 2) # minimum R is 2
	R2inv[i] = 1 / R**2.0;

	# measure raw (without wilson flow)
	w00 = wilsonLoops[0, i * 3]
	w01 = wilsonLoops[0, i * 3 + 1]
	w10 = wilsonLoops[0, i * 3 + 2]
	w11 = wilsonLoops[0, i * 3 + 3]
	chiRaw[i][0], chiRaw[i][1] = jackknifeCreutz(w00, w01, w10, w11)
	print("chiRaw(%d) = %.12f (%.12f)" % (R, chiRaw[i][0], chiRaw[i][1]))

	# measure smeared (with wilson flow)
	T = (R - 2.0)**2.0 / 4; # smearing radius is sqrt(2dt) in d dimensions
	t = min(int(T * 50), n_wf - 1)
	w00 = wilsonLoops[t, i * 3]
	w01 = wilsonLoops[t, i * 3 + 1]
	w10 = wilsonLoops[t, i * 3 + 2]
	w11 = wilsonLoops[t, i * 3 + 3]
	chiSmeared[i][0], chiSmeared[i][1] = jackknifeCreutz(w00, w01, w10, w11)
	print("chiSmeared(%d) = %.12f (%.12f)" % (R, chiSmeared[i][0], chiSmeared[i][1]))

print("\nFitting string tension...")

# fit Creutz ratios to chi(R) = sigma + C / R^2
def sigmaFit(R2inv, C, sigma):
	return sigma + C * R2inv

popt, pcov = opt.curve_fit(sigmaFit, R2inv, chiRaw[:,0], [0.01, 0.1], sigma=chiRaw[:,1])
CFitRaw = (popt[0], np.sqrt(pcov[0][0]))
sigmaFitRaw = (popt[1], np.sqrt(pcov[1][1]))

popt, pcov = opt.curve_fit(sigmaFit, R2inv, chiSmeared[:,0], [0.01, 0.1], sigma=chiSmeared[:,1])
CFitSmeared = (popt[0], np.sqrt(pcov[0][0]))
sigmaFitSmeared = (popt[1], np.sqrt(pcov[1][1]))

R2invFit = np.linspace(0.0, R2inv[0] * 1.2, 1000)
chiFitRaw = sigmaFit(R2invFit, CFitRaw[0], sigmaFitRaw[0])
chiFitSmeared = sigmaFit(R2invFit, CFitSmeared[0], sigmaFitSmeared[0])

print("sigmaFitRaw = %.12f (%.12f)" % (sigmaFitRaw[0], sigmaFitRaw[1]))
print("sigmaFitSmeared = %.12f (%.12f)" % (sigmaFitSmeared[0], sigmaFitSmeared[1]))

print("\nPlotting...")

plt.rcParams.update({
    "text.usetex": False,
    "font.family": "sans-serif"})

chiMax = max(np.amax(chiRaw[:,0]), np.amax(chiSmeared[:,0]))

# plot raw Creutz ratio vs 1/R^2
plt.figure()
plt.xlim(0, R2inv[0] * 1.2)
plt.ylim(0.0, chiRaw[0,0] * 1.5)
# plt.ylim(0.0, chiMax * 1.5)
plt.errorbar(R2inv, chiRaw[:,0], yerr=chiRaw[:,1], color='b', mec='b', marker='o', ms=5, mew=0.5, mfc='none', linestyle='none', linewidth=0.5, capsize=2.5, capthick=0.5)
plt.plot(R2invFit, chiFitRaw, color='b', linewidth=0.5, linestyle='dashed', label='Raw')
plt.errorbar(R2inv, chiSmeared[:,0], yerr=chiSmeared[:,1], color='g', mec='g', marker='o', ms=5, mew=0.5, mfc='none', linestyle='none', linewidth=0.5, capsize=2.5, capthick=0.5)
plt.plot(R2invFit, chiFitSmeared, color='g', linewidth=0.5, linestyle='dashed', label='Smeared')
plt.legend(loc='lower right', fontsize=10, handlelength=2.5, frameon=False)
plt.xlabel("$1/R^2$")
plt.ylabel("$\chi(R)$")
plt.savefig("%s/creutz.pdf" % (path))
plt.close()
