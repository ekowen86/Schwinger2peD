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


def creutzRatio(w00, w01, w10, w11):
	return -np.log(np.mean(w00) * np.mean(w11) / np.mean(w01) / np.mean(w10))


def jackknifeCreutz(w00, w01, w10, w11):
	n = len(w00) # number of data values (int)
	f_n = float(n) # number of data values (float)
	mean = creutzRatio(w00, w01, w10, w11)
	if mean != mean: # check for nan
		return (0.0, float('inf'))
	err = 0.0

	for i in range(0, n):

		# copy each array and delete the current value
		del00 = np.delete(w00, i)
		del01 = np.delete(w01, i)
		del10 = np.delete(w10, i)
		del11 = np.delete(w11, i)

		# calculate the creutz ratio and add to the error
		err += (creutzRatio(del00, del01, del10, del11) - mean)**2.0

	err = np.sqrt((f_n - 1) / f_n * err)
	return (mean, err)


def autocorrGamma(a, n):

	N = len(a)
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

	for n in range(1, len(a)):
		curGamma = autocorrGamma(a, n)
		if (curGamma < 0):
			break
		result += curGamma

	return result / Gamma0
