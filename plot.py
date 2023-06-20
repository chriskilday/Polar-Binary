""" plot.py """

import sys
import math
#from array import *
import os
import numpy as np
import matplotlib.pyplot as plt
#from random import randint
#from random import random
from input import *

particle_file = "out_particles/particles.txt"		# Path to input file
break_point = 0.9	# Used for a_avg and e_avg. Averaging will ignore first (break_point*100)% of time integrated
nbuckets = 100
tree_plotn = 100	# Number of frames in tree.mp4
depth = 8		# Depth of gravity tree searched when plotting density
PLOT_N = 1e1		# Number of particles shown when plotting positions

ffmpeg_path = "ffmpeg"	# Path to ffmpeg executable

def main():
	global particle_file
	if (len(sys.argv) == 1): return

	commands = sys.argv[1:]
	for i in range(len(commands)):
		if (commands[i] == '-positions'):
			plot_positions(not (len(commands)-1 > i and commands[i+1] == '-b'))
		elif (commands[i] == '-n'):
			plot_noft()
		elif (commands[i] == '-e_avg'):
			plot_eavg()
		elif (commands[i] == '-a_avg'):
			plot_aavg()
		elif (commands[i] == '-tree'):
			plot_tree()
		elif (commands[i] == '-i'):
			particle_file = commands[i+1]

def plot_noft():
	(noft, t) = read_noft(particle_file)
	fig, nax = plt.subplots()

	nax.plot(t[:len(noft)], noft)

	nax.set(title="N of t")
	plt.show()

def plot_positions(particlesp):
	(n_initial, t, noft, particles, stars) = read_particles(particle_file)
	fig, (xaxs, yaxs, zaxs, raxs) = plt.subplots(4, 1, sharex=True)

	if (particlesp):
		count = 0
		for particle in particles:
			t_particle = t[:particle.max_t()]
			xaxs.plot(t_particle, particle.get_x())
			yaxs.plot(t_particle, particle.get_y())
			zaxs.plot(t_particle, particle.get_z())
			raxs.plot(t_particle, particle.get_a())

			count += 1
			if (count >= PLOT_N): break
	else:
		for star in stars:
			t_star = t[:star.max_t()]
			xaxs.plot(t_star, star.get_x())
			yaxs.plot(t_star, star.get_y())
			zaxs.plot(t_star, star.get_z())

		if (len(stars) == 2):
			bx = [stars[0].get_x(), stars[1].get_x()]
			by = [stars[0].get_y(), stars[1].get_y()]
			bz = [stars[0].get_z(), stars[1].get_z()]
			raxs.plot(t[:min(stars[0].max_t(), stars[1].max_t())], \
				[math.sqrt((bx[0][i]-bx[1][i])**2 + (by[0][i]-by[1][i])**2 + (bz[0][i]-bz[1][i])**2) \
				for i in range(min(stars[0].max_t(), stars[1].max_t()))])

	xaxs.set_ylabel("X"); yaxs.set_ylabel("Y"); zaxs.set_ylabel("Z"); raxs.set_ylabel("A")
	fig.suptitle("Positions")
	plt.show()

def plot_aavg():
	(n_initial, t, noft, particles, stars) = read_particles(particle_file)
	break_index = math.floor(break_point * len(t))

	mina,maxa = float('inf'),0
	maxn = 0
	buckets = [0]*nbuckets
	for particle in particles:
		for i in range(break_index, particle.max_t()):
			if (particle.get_e()[i] >= 1 or particle.get_a()[i] > 30): continue
			if (particle.get_a()[i] > maxa): maxa = particle.get_a()[i]
			if (particle.get_a()[i] < mina): mina = particle.get_a()[i]

	for particle in particles:
		for i in range(break_index, particle.max_t()):
			if (particle.get_e()[i] >= 1 or particle.get_a()[i] > 30): continue
			buckets[math.floor((particle.get_a()[i]-mina) / (1.01*(maxa-mina)) * nbuckets)] += 1

	for i in range(nbuckets):
		buckets[i] /= len(t) - break_index
		if (buckets[i] > maxn): maxn = buckets[i]

	plt.scatter(np.linspace(mina, maxa, nbuckets), buckets, c= "red", s=10)
	plt.xscale('log'); plt.yscale('log')
	plt.gca().set(xlabel='Semi-major Axis', ylabel='N Particles', xlim=(0.9*mina,1.1*maxa), ylim=(0.9,1.1*maxn))
	plt.show()

def plot_eavg():
	(n_initial, t, noft, particles, stars) = read_particles(particle_file)
	break_index = math.floor(break_point * len(t))

	mina,maxa = 0,0
	mine,maxe = 0,-float('inf')
	for particle in particles:
		if (particle.max_t() == 0): continue
		for i in range(break_index, particle.max_t()):
			if (particle.get_e()[i] > 0.9): continue
			if (particle.get_a()[i] > maxa): maxa = particle.get_a()[i]
			if (particle.get_e()[i] > maxe): maxe = particle.get_e()[i]

	buckets = [[] for i in range(nbuckets)]
	for particle in particles:
		for i in range(break_index, particle.max_t()):
			if (particle.get_e()[i] > 0.9): continue
			buckets[math.floor((particle.get_a()[i]-mina) / (1.01*(maxa-mina)) * nbuckets)].append(particle.get_e()[i])

	e_out = [0]*nbuckets
	percentile_low, percentile_high = [0]*nbuckets, [0]*nbuckets
	for i in range(nbuckets):
		if (len(buckets[i]) == 0): continue

		e_out[i] = sum(buckets[i]) / len(buckets[i])
		percentile_low[i] = np.percentile(buckets[i], 10)
		percentile_high[i] = np.percentile(buckets[i], 90)

	plt.scatter(np.linspace(mina,maxa,nbuckets), e_out, c="red", s=10)
	plt.errorbar(np.linspace(mina,maxa,nbuckets), e_out, yerr=[percentile_low,percentile_high], fmt='none', ecolor='red')
	plt.gca().set(xlabel="Semi-major Axis", ylabel="Eccentricity", xlim=(mina,maxa), ylim=(mine,maxe))
	plt.show()

def plot_tree(plot_disc=True, plot_lost=False, plot_binary=False):
	(n_initial, t, noft, particles, stars) = read_particles(particle_file)
	lost_particles = read_lost_particles(particles, len(t))

	for i in range(tree_plotn):
		curt = math.floor(i * len(t) / (1.01*(tree_plotn - 1))) if (tree_plotn > 1) else 0

		if (plot_disc):
			root = read_tree(particles, curt)
			xyplane = [[0 for d1 in range(2**depth)] for d2 in range(2**depth)]
			zyplane = [[0 for d1 in range(2**depth)] for d2 in range(2**depth)]
			nodes = [root]
			next = []
			for d in range(depth):
				for node in nodes:
					for child in node.children:
						if (isinstance(child, Tree_Node)): next.append(child)
						elif (isinstance(child, Particle)):
							x = round((child.get_x()[curt] - tree_width/2**(depth+1)) / (tree_width / 2**depth) + 2**(depth-1))
							y = round((child.get_y()[curt] - tree_width/2**(depth+1)) / (tree_width / 2**depth) + 2**(depth-1))
							z = round((child.get_z()[curt] - tree_width/2**(depth+1)) / (tree_width / 2**depth) + 2**(depth-1))

							if (x < 0 or y < 0 or x >= 2**depth or y >= 2**depth):
								continue

							plt.figure(1); xyplane[x][y] = 1
							plt.figure(2); zyplane[z][y] = 1

				nodes = [node for node in next]
				next = []

			for node in nodes:
				x = int((node.x - tree_width/2**(depth+1)) / (tree_width / 2**depth) + 2**(depth-1))
				y = int((node.y - tree_width/2**(depth+1)) / (tree_width / 2**depth) + 2**(depth-1))
				z = int((node.z - tree_width/2**(depth+1)) / (tree_width / 2**depth) + 2**(depth-1))

				xyplane[x][y] += node.mass
				zyplane[z][y] += node.mass

			xymaxm = np.amax(xyplane)
			zymaxm = np.amax(zyplane)
			for d1 in range(2**depth):
				for d2 in range(2**depth):
					xyplane[d1][d2] /= xymaxm
					zyplane[d1][d2] /= zymaxm

		if (plot_lost):
			for particle in lost_particles:
				if (particle.max_t() > curt):
					x = (particle.get_x()[curt] - tree_width/2**(depth+1)) / (tree_width / 2**depth) + 2**(depth-1)
					y = (particle.get_y()[curt] - tree_width/2**(depth+1)) / (tree_width / 2**depth) + 2**(depth-1)
					z = (particle.get_z()[curt] - tree_width/2**(depth+1)) / (tree_width / 2**depth) + 2**(depth-1)

					plt.figure(1); plt.plot(x, y, 'go', zorder=1, markersize=1)
					plt.figure(2); plt.plot(y, z, 'go', zorder=1, markersize=1)

		if (plot_binary):
			for star in stars:
					x = (star.get_x()[curt] - tree_width/2**(depth+1)) / (tree_width / 2**depth) + 2**(depth-1)
					y = (star.get_y()[curt] - tree_width/2**(depth+1)) / (tree_width / 2**depth) + 2**(depth-1)
					z = (star.get_z()[curt] - tree_width/2**(depth+1)) / (tree_width / 2**depth) + 2**(depth-1)

					plt.figure(1); plt.plot(x, y, 'co', zorder=1, markersize=1)
					plt.figure(2); plt.plot(y, z, 'co', zorder=1, markersize=1)
				

		if (plot_disc):
			plt.figure(1); plt.imshow(xyplane, cmap = "coolwarm"); plt.xlim(0,2**depth); plt.ylim(0,2**depth)
			plt.figure(2); plt.imshow(zyplane, cmap = "coolwarm"); plt.ylim(2**(depth-1)-2**depth/8,2**(depth-1)+2**depth/8)
		else:
			plt.figure(1); plt.xlim(2**(depth-1)-2**depth/8,2**(depth-1)+2**depth/8); plt.ylim(2**(depth-1)-2**depth/8,2**(depth-1)+2**depth/8)
			plt.figure(2); plt.xlim(2**(depth-1)-2**depth/8,2**(depth-1)+2**depth/8); plt.ylim(2**(depth-1)-2**depth/8,2**(depth-1)+2**depth/8)
		plt.figure(1); plt.tight_layout(); plt.savefig("treexy" + (str(i) if (i > 99) else '0' + str(i) if (i > 9) else '00' + str(i)) + ".png", dpi=300)
		plt.figure(2); plt.tight_layout(); plt.savefig("treeyz" + (str(i) if (i > 99) else '0' + str(i) if (i > 9) else '00' + str(i)) + ".png", dpi=300)
		plt.figure(1); plt.cla()
		plt.figure(2); plt.cla()

	os.system("rm -f treexy.mp4 && " + ffmpeg_path + " -hide_banner -loglevel error -pattern_type glob -r 5 -i 'treexy*.png' treexy.mp4 && rm treexy*.png")
	os.system("rm -f treeyz.mp4 && " + ffmpeg_path + " -hide_banner -loglevel error -pattern_type glob -r 5 -i 'treeyz*.png' treeyz.mp4 && rm treeyz*.png")

main()
