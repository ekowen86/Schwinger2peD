import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import cmath
import sys
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)


L = int(16) # lattice size
# NZ = int(1) # number of slices in extra dimension
# beta = 5.00 # coupling
# epsZ = 1.00 # coupling ratio in extra dimension
R_half = int(L / 2) # half lattice size
R_min = 8 # min R value for fit
R_max = R_half # max R value for fit


def parse_data_file(file):
    lines = file.readlines()
    # a = np.empty((R_half+1,len(lines)))
    a = np.empty((R_half+1,500))
    i = 0
    for l, line in enumerate(lines):
        if (l >= len(a[0])):
            break;
        tokens = line.split()
        for t, token in enumerate(tokens):
            if (t == 0):
                continue # skip trajectory id
            a[t-1,l] = float(token)
    return a


def M_fit(m, mc, beta):
    return 2.008 * (m + mc)**(2./3.) * beta**(-1./6.)


# beta = 2.0
# m_array = np.array([0.0, 0.1, 0.2, 0.3])
# M = [0.377967182,0.589456907,0.766733317,0.924566522,1.048216867]
# d_M = [0.017962195,0.010350555,0.015477475,0.011885409,0.009947168]

# beta = 5.0
m_array = np.array([0.0, 0.05, 0.1, 0.2, 0.3, 0.5])
M = [0.334683013516,0.401868769909,0.492526228113,0.670103004281,0.845283075913,1.080287672987]
d_M = [0.012212786704,0.012625537171,0.008760219893,0.007238093454,0.008984428864,0.006494529616]


# calculate best fit for M vs m
popt, pcov = opt.curve_fit(M_fit, m_array, M, [0.2, 2.0], sigma=d_M)
mc = popt[0]
d_mc = np.sqrt(pcov[0][0])
beta_eff = popt[1]
d_beta_eff = np.sqrt(pcov[1][1])

print("mc       = %.12f (%.12f)" % (mc, d_mc))
print("beta_eff = %.12f (%.12f)" % (beta_eff, d_beta_eff))

m_A = np.linspace(-mc, 1.2, 1000)
M_A = np.zeros(len(m_A))
for m in range(1, len(m_A)):
	M_A[m] = M_fit(m_A[m], mc, beta_eff)

plt.rcParams.update({
    "text.usetex": False,
    "font.family": "sans-serif"})

# plot M vs m
plt.figure()
plt.xlim(-0.3, 0.4)
plt.errorbar(m_array, M, yerr=d_M, color="blue", marker='o', ms=5, mew=0.5, mfc='none', linestyle='none', linewidth=0.5, capsize=2.5, capthick=0.5)
plt.plot(m_A, M_A, color="blue", linewidth=0.5)
plt.xlabel("$m_0$")
plt.ylabel("$M_{\pi}$")
plt.savefig("beta_eff.pdf")
plt.close()
