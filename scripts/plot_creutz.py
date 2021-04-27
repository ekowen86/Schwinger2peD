import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.optimize as opt
import pandas as pd
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
    id = "%d_%d_%d_%d_%s%d" % (L, Lz, round(beta * 1000), round(eps3 * 1000), \
		m_sign, round(abs(m_fermion) * 1000))
    path = "../jobs/3D/%s" % id
print("id: %s" % (id))

t_max = float(sys.argv[6])
print("t_max: %f" % (t_max))

dt = float(sys.argv[7])
print("dt: %f" % (dt))

r_min = int(sys.argv[8]) - 2
print("r_min: %d" % (r_min + 2))

r_max = int(sys.argv[9]) - 2
print("r_max: %d" % (r_max + 2))

# maximum wilson loop size
loopMax = L / 2

# number of wilson flow values
n_wf = int(t_max / dt + 1)


# read creutz ratios from file
creutz = pd.read_csv("%s/creutz.dat" % path, header=None, delim_whitespace=True)

print("\nCalculating Creutz Ratios (smeared/raw)...")
R2inv = np.empty(loopMax - 1)
chiRaw = np.empty((len(R2inv), 2))
chiSmeared = np.empty((len(R2inv), 2))

for i in range(0, len(R2inv)):
	R = float(i + 2) # minimum R is 2
	R2inv[i] = 1 / R**2.0;

	# measure raw (without wilson flow)
	chiRaw[i][0] = creutz.iloc[0, i * 2 + 1]
	chiRaw[i][1] = creutz.iloc[0, i * 2 + 2]

	# measure smeared (with wilson flow)
    # wilson flow smearing radius is sqrt(2dt) in d dimensions
	# T = T_smear[i + 2]
	T = max(0.0, (R - 1.0) / 2.0)**2.0 / 4.0;
	t = min(int(T / dt), n_wf - 1)
	chiSmeared[i][0] = creutz.iloc[t, i * 2 + 1]
	chiSmeared[i][1] = creutz.iloc[t, i * 2 + 2]

print("\nFitting string tension...")

# fit Creutz ratios to chi(R) = sigma + C / R^2
def sigmaFit(R2inv, C, sigma):
	return sigma + C * R2inv

# popt, pcov = opt.curve_fit(sigmaFit, R2inv, chiRaw[:,0], [0.01, 0.1], sigma=chiRaw[:,1])
# CFitRaw = (popt[0], np.sqrt(pcov[0][0]))
# sigmaFitRaw = (popt[1], np.sqrt(pcov[1][1]))

popt, pcov = opt.curve_fit(sigmaFit, R2inv[r_min:r_max], chiSmeared[r_min:r_max,0], \
	[0.01, 0.1], sigma=chiSmeared[r_min:r_max,1])
CFitSmeared = (popt[0], np.sqrt(pcov[0][0]))
sigmaFitSmeared = (popt[1], np.sqrt(pcov[1][1]))

R2invFit = np.linspace(0.0, R2inv[0] * 1.2, 1000)
# chiFitRaw = sigmaFit(R2invFit, CFitRaw[0], sigmaFitRaw[0])
chiFitSmeared = sigmaFit(R2invFit, CFitSmeared[0], sigmaFitSmeared[0])

# print("sigmaFitRaw = %.12f (%.12f)" % (sigmaFitRaw[0], sigmaFitRaw[1]))
print("sigmaFitSmeared = %.12f (%.12f)" % (sigmaFitSmeared[0], sigmaFitSmeared[1]))

result_file = open("../jobs/string_tension.dat", "a")
result_file.write("%d %d %.3f %.3f %.3f %.12e %.12e\n" % \
	(L, Lz, beta, eps3, m_fermion, sigmaFitSmeared[0], sigmaFitSmeared[1]))
result_file.close()

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
plt.errorbar(R2inv, chiRaw[:,0], yerr=chiRaw[:,1], color='b', mec='b', \
	marker='o', ms=5, mew=0.5, mfc='none', linestyle='none', linewidth=0.5, \
	capsize=2.5, capthick=0.5, label='Raw')
# plt.plot(R2invFit, chiFitRaw, color='b', linewidth=0.5, linestyle='dashed')
plt.errorbar(R2inv, chiSmeared[:,0], yerr=chiSmeared[:,1], color='g', mec='g', \
	marker='o', ms=5, mew=0.5, mfc='none', linestyle='none', linewidth=0.5, \
	capsize=2.5, capthick=0.5, label='Smeared')
plt.plot(R2invFit, chiFitSmeared, color='g', linewidth=0.5, linestyle='dashed')
plt.legend(loc='lower right', fontsize=10, handlelength=2.5, frameon=False)
plt.xlabel("$1/R^2$")
plt.ylabel("$\chi(R,R)$")
plt.savefig("%s/creutz.pdf" % (path))
plt.close()
