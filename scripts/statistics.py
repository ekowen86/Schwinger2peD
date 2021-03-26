import numpy as np


def jackknifeMean(a):
	n = len(a) # number of data values (int)
	f_n = float(n) # number of data values (float)
	mean = np.mean(a) # average
	err = 0.0 # error

	for i in range(0, n):

		# copy the array, deleting the current value
		# and add to the error
		err += (np.mean(np.delete(a, i)) - mean)**2.0

	err = np.sqrt((f_n - 1) / f_n * err)
	return (mean, err)


def autocorrGamma(a, n):

	result = 0.0
	mean = np.mean(a)
	start = 0
	end = N - n
	
	# if n is negative, we have to start the loop at a later point
	if (n < 0):
		start = -n
		end = N

	for i in range(start, end):
		result += (a[i] - mean) * (a[i + n] - mean);

	return result / (end - start);


def autocorrTime(a):

	Gamma0 = autocorrGamma(a, 0);
	result = 0.5 * Gamma0

	for n in range(1, N):
		curGamma = autocorrGamma(a, n)
		if (curGamma < 0):
			break
		result += curGamma

	return result / Gamma0
