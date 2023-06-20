C:
The makefile creates an executable called "pb.exe" which can be run with any combination of the following flags.
	Ex. mpirun -np 9 pb.exe -pol -e 0.5 -discm 0.001 -o out_particles
Will output to a file named "particles.txt" to the specified directory. Any existing directory of the same name will be deleted.

Flags:
	-t <double>: set the max time integrated to
	-mr <double>: set the mass ratio of the binary
	-a <double>: set the sem-major axis of the binary
	-e <double>: set the eccentricity of the binary
	-n <integer>: set initial number of particles in disc
	-discm <double>: set the total mass of the disc
	-inner <double>: set inner edge of disc
	-outer <double>: set outer edge of disc
	-dt <double>: set integration time step size
	-o <string>: set the name of the output directory
	-oi <integer>: set the output interval
	-copl: make the system coplanar
	-pol: make the system polar
	-!: shorthand for "-discm 0.0"

The default values at the head of polar_binary.c will be used for any parameters not set by the flags.
The total mass of the binary is always 1.0. m1 = mr. m2 = 1-mr
Setting a mass ratio of 0 or 1 will replace the binary with a single star of mass 1.0.
"outer" and "inner" are set in units of the binary aphelion.

BOXN controlls the size of the gravity tree. I have been using a 9x9 tree which allows me to run with 9 processors for testing on elgato and 27 when running on ocelote. If you want to use a different number of processors, you may have to change this value. The binary orbit must remain on a single processor, so only odd numbers work.


Python:
Use plot.py with the following flags to plot output. If multiple flags are used, they will be plot sequentially, not at the same time.
	Ex. python3 plot.py -i out_particles/particles.txt -positions

Flags:
	-n: plots number of particles in simulations as a function of time
	-positions: plots positions of particles in disc as a fuction of time. (x, y, z, and semi-major axis)
	-positions -b: plots positions of stars in binary as a fuction of time. (x, y, z, and distance apart)
	-tree: creates 2 mp4 files: treexy.mp4, treeyz.mp4. Shows 2D projection of simlation.
	-a_avg: plots number of particles vs. semi-major axis
	-e_avg: plots eccentricity vs. semi-major axis
	-i: path to input file. must be used before any other flags, or default specified in plot.py will be used

plot_tree() has 3 parameters. plot_disc will show the density of the disc. plot_lost will show the particles lost by the end of the simulation. plot_binary will show the stars.
tree_plotn controlls the number of snapshots plotted for the mp4. These snapshots will be evenly spaced from beginning to end of the simulation. i.e. if there are 100 total time units output by the integration, tree_plotn=100 will show every time unit, tree_plotn=25 will show every 4th time unit...
