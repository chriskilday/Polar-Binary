""" plot.py """

import sys
import math
from array import *
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import bisect
from input import *

particle_file = "particles.txt"
treefile = "tree"

ecc_break = 0.1
PLOT_N = 1e2

colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']

def main():
	if (len(sys.argv) == 1): return

	global n_ecc
	commands = sys.argv[1:]
	for i in range(len(commands)):
		if (commands[i] == '-positions'):
			plot_positions(not (len(commands)-1 > i and commands[i+1] == '-b'))
		elif (commands[i] == '-orbits'):
			plot_orbits()
		elif (commands[i] == '-eccentricity'):
			plot_eccentricity()
		elif (commands[i] == '-density'):
			plot_density()
		elif (commands[i] == '-n'):
			plot_noft()

def plot_noft():
	(noft, t) = read_noft(particle_file)
	fig, nax = plt.subplots()

	nax.plot(t[:len(noft)], noft)

	nax.set(title="N of t")
	plt.show()

def plot_positions(particlesp):
	(n, x, y, z, bx, by, bz, t) = read_positions(particle_file)
	fig, (xaxs, yaxs, zaxs) = plt.subplots(3, 1, sharex=True)

	if (particlesp):
		for i in range(int(min(PLOT_N, n))):
			t_temp = t[:len(x[i])]
			xaxs.plot(t_temp, x[i]); yaxs.plot(t_temp, y[i]); zaxs.plot(t_temp, z[i])
	else:
		t = t[:len(bx[0])]
		xaxs.plot(t, bx[0])
		yaxs.plot(t, by[0])
		zaxs.plot(t, bz[0])

		if (len(bx[1]) != 0):
			xaxs.plot(t, bx[1])
			yaxs.plot(t, by[1])
			zaxs.plot(t, bz[1])

	xaxs.set_ylabel("X"); yaxs.set_ylabel("Y"); zaxs.set_ylabel("Z");
	fig.suptitle("Positions")
	plt.show()

def plot_orbits():
	(n, a, e, inc, Omega, t) = read_orbits(particle_file)
	fig, (aaxs, eaxs, iaxs, Oaxs) = plt.subplots(4, 1, sharex=True)

	for i in range(int(min(PLOT_N, n))):
		t_temp = t[:len(a[i])]
		aaxs.plot(t_temp, a[i]);
		eaxs.plot(t_temp, e[i]);
		iaxs.plot(t_temp, [inci / math.pi for inci in inc[i]]);
		Oaxs.plot(t_temp, [Omegai / math.pi for Omegai in Omega[i]])

	aaxs.set_ylabel("A")
	eaxs.set_ylabel("E")
	iaxs.set_ylabel("I/pi")
	Oaxs.set_ylabel("Omega/pi")
	fig.suptitle("Orbits")
	plt.show()

def plot_eccentricity():
	mr_last = 0
	i = -1
	for file in os.listdir():
		if (not (file.startswith("particles") and file.endswith(".txt"))): continue
		i += 1

		(n, a, e, inc, Omega, t) = read_orbits(file)
		a_avg = average_all(a, ecc_break)
		e_avg = average_all(e, ecc_break)

		(mr,ab,eb) = read_binary(file); e_last = eb
		mu = min(1,mr) / (1+mr)
		legend = "m1/m2 = " + str(1/mr)
		plt.plot(a_avg, e_avg, colors[i%len(colors)]+".", label=legend)

		a_sorted = list(filter(lambda item: item is not None, a_avg))
		a_sorted.sort()
		plt.plot(a_sorted, predicted_fit(a_sorted, ab, eb, mu), \
			 colors[i%len(colors)]+'--');

	plt.title("Binary Eccentricity: " + str(e_last))
	plt.xlabel("Semi_Major Axis")
	plt.ylabel("Eccentricity")
	plt.legend()
	plt.show()

"""
	for i in range(n_ecc):
		tag = "_ecc" + str(i+1) + ".txt"
		ascii_file = "ascii" + tag; orbit_file = "orbits" + tag; energy_file = "energy" + tag
		(n,t,x,v,a,e,I,Omega,omega) = read_files(ascii_file, orbit_file, energy_file)

		(m1, m2, ab, eb) = binary_initials(ascii_file)
		mu = min(m1,m2) / (m1 + m2)

		a_avg = average_all(a[1:], ecc_break); e_avg = average_all(e[1:], ecc_break)
		plt.plot(a_avg, e_avg, colors[i%len(colors)]+".")
		plt.plot(a_avg, best_fit(a_avg, e_avg, ab, eb, mu), colors[i%len(colors)]+'-')
		plt.plot(a_avg, predicted_fit(a_avg, ab, eb, mu), colors[i%len(colors)]+'--');

	plt.title("Polar: e = 0.5")
	plt.xlabel("Semi-Major Axis")
	plt.ylabel("Eccentricity")
	plt.show()
"""

tree_roots = []
depth = [8,8,8,8,8,7,7,7,7,6]
def plot_density():
	num_tree = sum(1 for file in os.listdir() \
			if (file.startswith(treefile+'0_') and file.endswith(".txt")))
	num_out = sum(1 for file in os.listdir() \
			if (file.startswith(treefile) and file.endswith('_0'+".txt")))

	fig, ax_vec = plt.subplots(math.ceil(num_out/5), min(5,num_out))

	for i in range(num_out):
		tree_roots = read_tree(treefile+str(i), num_tree)

		mass = []
		for root in tree_roots:
			add_masses(root, mass, depth[i], 1)

		densities = {}
		x, y = [], []
		max_density = 0
		for m in mass:
			if ((m[0], m[1]) in densities):
				densities[(m[0], m[1])] += m[2] % 1
			else:
				densities[(m[0], m[1])] = m[2] % 1
			if (densities[(m[0], m[1])] > max_density):
				max_density = densities[(m[0], m[1])]

			if (m[0] not in x): bisect.insort(x, m[0])
			if (m[1] not in y): bisect.insort(y, m[1])

		size = max(len(x), len(y))
		if (len(x) < size): x = y
		if (len(y) < size): y = x
		plot_arr = [[0 for i in range(size)] for j in range(size)]
		for row in range(size):
			for col in range(size):
				if (row < len(x) and col < len(y) and \
				    (x[row], y[col]) in densities):
					plot_arr[row][col] = densities[(x[row],y[col])] / max_density

		if (num_out == 1):
			ax_vec.pcolor(x, y, plot_arr, cmap="coolwarm", shading="auto")
		elif (num_out <= 5):
			ax_vec[i].pcolor(x, y, plot_arr, cmap="coolwarm", shading="auto")
		else:
			ax_vec[i//5][i%5].pcolor(x,y, plot_arr, cmap="coolwarm", shading="auto")

	plt.show()

def add_masses(node, mass, max_depth, cur_depth):
	if (cur_depth == max_depth):
		mass.append((node.x,node.y,node.mass))
	else:
		for i in range(8):
			if (node.oct[i] != None):
				add_masses(node.oct[i], mass, max_depth, cur_depth+1)

def fill_area(x, y, density, size):
	xs = [x - size/2.0, x - size/2.0, x + size/2.0, x + size/2.0]
	ys = [y - size/2.0, y + size/2.0, y + size/2.0, y - size/2.0]

	plt.fill(xs, ys, 'r', alpha=density)

def best_fit(x, y, ab, eb, mu):
	coefs, covar = curve_fit(fitting_function, x, y)

	print(average(coefs) / (ab * eb * (1 - 2*mu)))

	x.sort()
	return [coefs[0] / x[i] for i in range(len(x)) if x[i] != None]

def fitting_function(x, c):
	return c / x

def predicted_fit(x, ab, eb, mu):
	const = ab*eb * 5 * (1 - 2*mu) / 4;
	return [ const / xi for xi in x ]

def average_all(list, vbreak):
	return [average(i) for i in list[math.floor(vbreak*len(list)):] if (i != None)]

def average(list):
	sum = 0; len = 0
	for i in list:
		sum += i; len += 1
	return sum/len

main()
