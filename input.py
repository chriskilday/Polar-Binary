""" input.py """

def read_files(ascii_file, orbit_file, energy_file):
	fascii = open(ascii_file, 'r')
	forbit = open(orbit_file, 'r')
	fenergy = open(energy_file, 'r')

	n = get_n(fascii); t = []
	x = [None]*n; v = [None]*n
	a = [None]*(n-1); e = [None]*(n-1); I = [None]*(n-1)
	Omega = [None]*(n-1); omega = [None]*(n-1)

	for i in range(n):
		x[i] = []; v[i] = []
		if (i<n-1):
			a[i]=[]; e[i]=[]; I[i]=[]; Omega[i]=[]; omega[i]=[]

	line_num = 0
	order = []
	for line in forbit.readlines():
		if (line.startswith("***")): continue
		if (line.startswith("^^^")): order = parse_order(line.strip("^ "), n); continue
		if (line.startswith("---")):
			a[get_removed(line)+1] = None
			e[get_removed(line)+1] = None
			continue

		while (line_num % (n-1) >= 1 and order[(line_num % (n-1))-1] == None):
			line_num += 1
		if (len(t) == 0 or get_t(line) != t[len(t)-1]):
			t.append(get_t(line))

		if (line_num % (n-1) == 0):
			a[0].append(get_a(line))
			e[0].append(get_e(line))
			I[0].append(get_i(line))
			Omega[0].append(get_Omega(line))
			omega[0].append(get_omega(line))
		else:
			a[order[(line_num % (n-1))-1]+1].append(get_a(line))
			e[order[(line_num % (n-1))-1]+1].append(get_e(line))
			I[order[(line_num % (n-1))-1]+1].append(get_i(line))
			Omega[order[(line_num % (n-1))-1]+1].append(get_Omega(line))
			omega[order[(line_num % (n-1))-1]+1].append(get_omega(line))
		line_num += 1

	fascii.close(); forbit.close(); fenergy.close()
	return (n,t,x,v,a,e,I,Omega,omega)

def read_positions(particle_file):
	pt_file = open(particle_file, 'r')

	n = get_n(pt_file)
	x = []; y = []; z = []; t = []
	bx = [[],[]]; by = [[],[]]; bz = [[],[]]
	for i in range(n):
		x.append([]); y.append([]); z.append([])

	for line in pt_file.readlines():
		if (line.startswith("***")): continue
		if (line.startswith("@@@")):
			t.append(t_val(line))
			continue

		hash = hash_val(line)
		if (hash.startswith("Star")):
			index = int(hash[-1]) - 1
			bx[index].append(x_val(line))
			by[index].append(y_val(line))
			bz[index].append(z_val(line))
		else:
			hash = int(hash)

			x[hash].append(x_val(line))
			y[hash].append(y_val(line))
			z[hash].append(z_val(line))

	pt_file.close()

	for i in range(n):
		if (len(x[i]) < len(t)):
			x[i] = None
			y[i] = None
			z[i] = None

	return (n, x, y, z, bx, by, bz, t)


def read_orbits(particle_file):
	pt_file = open(particle_file, 'r')

	n = get_n(pt_file)
	a, e, i, Omega, t = [], [], [], [], []
	for j in range(n):
		a.append([]); e.append([]); i.append([]); Omega.append([])

	for line in pt_file.readlines():
		if (line.startswith("***")): continue
		if (line.startswith("@@@")):
			t.append(t_val(line))
			continue

		hash = hash_val(line)
		if (hash.startswith("Star")): continue
		hash = int(hash)

		a[hash].append(a_val(line))
		e[hash].append(e_val(line))
		i[hash].append(i_val(line))
		Omega[hash].append(Omega_val(line))

	pt_file.close()

	for j in range(n):
		if (len(a[j]) < len(t)):
			a[j] = None
			e[j] = None
			i[j] = None
			Omega[j] = None

	return (n, a, e, i, Omega, t)

def read_binary(particle_file):
	pt_file = open(particle_file, 'r')
	for line in pt_file.readlines():
		if (line.startswith("*** MR")):
			split = line.split()
			return (float(split[3]), float(split[6]), float(split[9]))

	return 0

def read_tree(treefile, n_nodes):
	tree_roots = [None] * n_nodes
	for i in range(n_nodes):
		filename = treefile + '_' + str(i)+ ".txt" 
		infile = open(filename, 'r')

		tree_roots[i] = read_node(infile.readline())
		for line in infile.readlines():
			tree_roots[i].add_child(read_node(line))

		infile.close()

	return tree_roots

def read_node(line):
	vals = line.split()[1::2]
	return Tree_Node(float(vals[0]), float(vals[1]), float(vals[2]), \
			 float(vals[3]), float(vals[4]))

class Tree_Node:
	x,y,z = 0.,0.,0.
	width = 0.
	mass = 0.

	oct = []

	def __init__(self, x, y, z, width, mass):
		self.x = x; self.y = y; self.z = z
		self.width = width
		self.mass = mass
		self.oct = [None] * 8

	def __repr__(self):
		return f"x: {self.x} y: {self.y} z: {self.z} " + \
			f"w: {self.width} m: {self.mass}"

	def add_child(self, node):
		if (node.width < 0.51 * self.width and \
		    node.width > 0.49 * self.width):
			self.oct[self.get_oct(node)] = node
		else:
			self.oct[self.get_oct(node)].add_child(node)

	def get_oct(self, child):
		oct = int("".join([ "1" if child.x > self.x else "0", \
		       		     "1" if child.y > self.y else "0", \
				     "1" if child.z > self.z else "0"]), \
			   2)

		return oct

	def get_mass(self, x, y, width):
		if (disjoint(x,y,width, self.x,self.y,self.width)):
			return 0
		if (subregion(self.x,self.y,self.width, x,y,width)):
			return self.mass

		out = 0
		for child in self.oct:
			if (child != None): out += child.get_mass(x,y,width)

#		if (out == 0): return self.mass * (width**2 / self.width**2)
		return out

def get_n(file):
	file.readline()
	return int(file.readline().strip(" *N="))

def get_t(line):
	return float(line.split()[0])

def get_a(line):
	return float(line.split()[1])

def get_e(line):
	return float(line.split()[2])

def get_i(line):
	return float(line.split()[3])

def get_Omega(line):
	return float(line.split()[4])

def get_omega(line):
	return float(line.split()[5])

def get_removed(line):
	return int(line.strip("- "))

def parse_order(line, n):
	out = [None] * (n-2)
	split = line.split()
	for i in range(len(split)):
		out[i] = int(split[i])
	return out

def binary_initials(filename):
	f = open(filename, 'r')
	line = f.readline()
	while(not line.startswith('#')): line = f.readline()

	vals = line.strip('# ').split()
	return (float(vals[0]), float(vals[1]), float(vals[2]), float(vals[3]))

def t_val(line):
	return float(line.strip(" @\n"))

def hash_val(line):
	return line.split()[0]

def x_val(line):
	return float(line.split()[1])

def y_val(line):
	return float(line.split()[2])

def z_val(line):
	return float(line.split()[3])

def a_val(line):
	return float(line.split()[4])

def e_val(line):
	return float(line.split()[5])

def i_val(line):
	return float(line.split()[6])

def Omega_val(line):
	return float(line.split()[7])
