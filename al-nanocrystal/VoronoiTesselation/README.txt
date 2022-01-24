Program for constructing nanocrystalline samples
Author: Pablo Piaggi
Should you have any questions or suggestions, please send an 
email to: ppiaggi@hotmail.com

Instructions:
1) Download source code with the instruction: git clone \
https://ppiaggi@bitbucket.org/ppiaggi/voronoitesselation.git \
VoronoiTesselation
2) Compile source: Type make
3) Edit input file input.in
4) Run main.x input.in
5) Use the created data file (md.data) and the grain size histogram
file (graindistrib.dat)

Parameters in input file:
1) avgGrainSize	= Desired average grain size.
2) latticeConst	= Lattice constant, a for hcp.
3) caRatio = c/a ratio for the hexagonal lattice.
4) bravaisLat	= Bravais lattice (b c c, h c p, or s c - mind the
spaces! Sorry for that.)
5) regionCells	A	B	C  = Box size in lattice constants
6) cells	A'	B'	C' = Number of cells for cell lists. The
recommended values are half the unprimed values (i.e. A'=A/2).
7) rCut		= Minimum distance between two particles. This parameter
controles the porosity.
8) rCutGrain	= Minimum grain size. In case you want to truncate the 
grain size distribution.
9) grainCellsFactor	= This factor determines the size of the crystal
constructed arround each seed. If too small, there will be porosity or 
empty space between grains. If to big, the program will take longer.
10) stepGD		= Step between bins in grain size histogram.
11) rangeGD		= Range around the average grain size where
the histogram is made. min/max = avgGrainSize -/+ rangeGD/2
12) randSeedP		= Random seed for positions and angles 
