""" input.py """

import os

class Particle:
	def __init__(self):
		self.position_vec = []
		self.velocity_vec = []
		self.vtotal_vec = []
		self.a_vec = []	
		self.e_vec = []

	def add_position(self, x, y, z):
		self.position_vec.append((x,y,z))

	def add_velocity(self, x, y, z, v):
		self.velocity_vec.append((x,y,z))
		self.vtotal_vec.append(v)

	def add_a(self, a):
		self.a_vec.append(a)

	def add_e(self, e):
		self.e_vec.append(e)

	def max_t(self):
		return len(self.position_vec)

	def get_x(self):
		return [pos[0] for pos in self.position_vec]

	def get_y(self):
		return [pos[1] for pos in self.position_vec]

	def get_z(self):
		return [pos[2] for pos in self.position_vec]

	def get_a(self):
		return self.a_vec

	def get_e(self):
		return self.e_vec

def read_particles(particle_filename):
	particle_file = open(particle_filename, 'r')

	n_initial = get_n_initial(particle_file)
	t = []; noft = []; particles = []; stars = []
	for i in range(n_initial):
		particles.append(Particle())
	stars.append(Particle()); stars.append(Particle())

	for line in particle_file.readlines():
		if (line.startswith("***") or line.startswith("@@@") or line.startswith("###")):
			if (line.startswith("@@@")):
				t.append(t_val(line))
			if (line.startswith("###")):
				noft.append(n_val(line))
			continue

		splitline = line.split()
		if (len(splitline) != 10): splitline.append('0')
		hash = splitline[0]
		x, y, z = float(splitline[1]), float(splitline[2]), float(splitline[3])
		vx, vy, vz = float(splitline[4]), float(splitline[5]), float(splitline[6])
		vtotal = float(splitline[7])
		a, e = float(splitline[8]), float(splitline[9])

		if (hash.startswith("Star")):
			hash = int(hash.strip("Star")) - 1

			stars[hash].add_position(x,y,z)
			stars[hash].add_velocity(vx,vy,vz,vtotal)
			stars[hash].add_a(a)
			stars[hash].add_e(e)
		else:
			hash = int(hash)

			particles[hash].add_position(x,y,z)
			particles[hash].add_velocity(vx,vy,vz,vtotal)
			particles[hash].add_a(a)
			particles[hash].add_e(e)

	particle_file.close()
	return (n_initial, t, noft, particles, stars)

def read_lost_particles(particles, maxt):
	out_particles = []
	for particle in particles:
		if (particle.max_t() < maxt):
			out_particles.append(particle)

	return out_particles

tree_width = 60
def read_tree(particles, time):
	root = Tree_Node(0,0,0,tree_width)

	for particle in particles:
		if (particle.max_t() <= time): continue
		root.add_particle(particle, time)

	return root

def read_noft(particle_filename):
	(n_initial, t, noft, particles, stars) = read_particles(particle_filename)
	return (noft, t)

def get_n_initial(file):
	file.readline()
	return int(file.readline().strip(" *N="))

def t_val(line):
	return float(line.strip(" @\n"))

def n_val(line):
	return int(line.strip(" #\n"))

class Tree_Node:
	def __init__(self, x, y, z, w):
		self.x = x
		self.y = y
		self.z = z
		self.width = w

		self.mass = 0
		self.children = [None] * 8

	def add_particle(self, particle, t):
		self.mass += 1
		index = 2**2 * (1 if (particle.get_x()[t] > self.x) else 0) \
			+ 2**1 * (1 if (particle.get_y()[t] > self.y) else 0) \
			+ 2**0 * (1 if (particle.get_z()[t] > self.z) else 0)

		if (self.children[index] == None): self.children[index] = particle
		elif (isinstance(self.children[index], Particle)):
			node = Tree_Node(self.x + self.width / 4 * (1 if (particle.get_x()[t] > self.x) else -1), \
					 self.y + self.width / 4 * (1 if (particle.get_y()[t] > self.y) else -1), \
					 self.z + self.width / 4 * (1 if (particle.get_z()[t] > self.z) else -1), \
					 self.width / 2)

			node.add_particle(self.children[index], t)
			node.add_particle(particle, t)

			self.children[index] = node
		elif (isinstance(self.children[index], Tree_Node)):
			self.children[index].add_particle(particle, t)

def contains(lst, elm):
	for i in lst:
		if (i == elm): return True
	return False
