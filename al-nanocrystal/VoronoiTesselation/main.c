// Program by Pablo Piaggi
// January-July 2014

#include "mddefs.h"

#define NDIM 3

int nMol, nMolMax, numAtomType, grainNumber, nBasis, *cellList, *atomsPerGrain, randSeedP;
real density, temperature, velMag, *mass, avgGrainSize, latticeConst, rot[3][3], rCut, rCutGrain, grainCellsFactor, stepGD, rangeGD, massIn, caRatio;
Mol *mol;
VecR region, primVect[3], *basisVect, grainRegion;
Nuclei *nuclei;
VecI regionCells, cells, grainCells;
Prop kinEnergy, totEnergy, potEnergy;
char bravaisLat[3];
FILE *pos, *data, *grain;

NameList nameList[] = {
	NameI (regionCells),
	NameR (avgGrainSize),
	NameR (latticeConst),
	NameR (caRatio),
	NameC (bravaisLat),
	NameR (rCutGrain),
	NameI (cells),
	NameR (rCut),
	NameR (grainCellsFactor),
	NameR (stepGD),
	NameR (rangeGD),
	NameR (massIn),
	NameI (randSeedP),
};

int main (int argc, char **argv)
{
	GetNameList (argc, argv);
	PrintNameList (stdout);
	SetParams ();
	SetupJob ();
	CreateNuclei ();
	CreateCrystal ();
	EraseOverlap ();
	data = fopen ( "md.data", "w");
	WriteData (data);
	fclose (data);
	grain = fopen ( "graindistrib.dat", "w");
	WriteGrainDistrib (grain);
	fclose (grain);
}

void SetParams ()
{
	int n;
	real c, newGrainCellsFactor;

	// New variable to circumvent problems with HCP
	newGrainCellsFactor = grainCellsFactor;
	if ( strcmp (bravaisLat, "bcc") == 0 ) {
			VSet (primVect[0], latticeConst , 0., 0.);
			VSet (primVect[1], 0., latticeConst, 0.);
			VSet (primVect[2], 0., 0., latticeConst);
			nBasis = 2;
			AllocMem (basisVect, nBasis, VecR);
		       	VSetAll (basisVect[0], 0.); 
		       	VSetAll (basisVect[1], latticeConst / 2.); 
			density = 2. / Cube (latticeConst);
	} else if ( strcmp (bravaisLat, "fcc") == 0 ) {
			VSet (primVect[0], latticeConst , 0., 0.);
			VSet (primVect[1], 0., latticeConst, 0.);
			VSet (primVect[2], 0., 0., latticeConst);
			nBasis = 4;
			AllocMem (basisVect, nBasis, VecR);
		       	VSetAll (basisVect[0], 0.); 
		       	VSet (basisVect[1], latticeConst / 2., latticeConst / 2., 0.); 
		       	VSet (basisVect[2], 0., latticeConst / 2., latticeConst / 2.); 
		       	VSet (basisVect[3], latticeConst / 2., 0., latticeConst / 2.); 
			density = 4. / Cube (latticeConst);
	} else if ( strcmp (bravaisLat, "sc") == 0 ) {
			VSet (primVect[0], latticeConst , 0., 0.);
			VSet (primVect[1], 0., latticeConst, 0.);
			VSet (primVect[2], 0., 0., latticeConst);
			nBasis = 1;
			AllocMem (basisVect, nBasis, VecR);
		       	VSetAll (basisVect[0], 0.); 
			density = 1. / Cube (latticeConst);
	} else if ( strcmp (bravaisLat, "hcp") == 0 ) {
			c=caRatio*latticeConst;
			VSet (primVect[0], latticeConst , 0., 0.);
			VSet (primVect[1], latticeConst/2. , pow (3.,1./3.)*
				       latticeConst/2.	, 0.);
			VSet (primVect[2], 0., 0., c);
			nBasis = 2;
			AllocMem (basisVect, nBasis, VecR);
		       	VSetAll (basisVect[0], 0.); 
		       	VSet (basisVect[1], latticeConst/2. , pow (3.,1./3.)*
				       latticeConst/6.	, c/2.);
			density = 2. / ( pow (3.,1./3.)* Sqr (latticeConst) * c/2.);
			// New variable to circumvent problems with HCP
			newGrainCellsFactor = grainCellsFactor * 1.3;
	}
	grainNumber=6. * Cube (latticeConst) * VProd (regionCells) / (Pi * 
			Cube (avgGrainSize));
	printf("Grain number %d \n", grainNumber);
	VSCopy (region, latticeConst, regionCells);
	printf("Region size %f %f %f \n", WriteVec (region));
	//velMag = sqrt (NDIM * (1. - 1. /nMol) * temperature);
	VSet (grainCells, avgGrainSize * newGrainCellsFactor  / VLen (primVect[0]),
			 avgGrainSize * newGrainCellsFactor  / VLen (primVect[1]),
			 avgGrainSize * newGrainCellsFactor  / VLen (primVect[2])
			);
	printf("grain Cells %d %d %d \n", WriteVec (grainCells));
	// VSCopy (grainRegion, latticeConst, grainCells);
	// printf("Grain region size %f %f %f \n", WriteVec (grainRegion));
	numAtomType = 1;
	nMolMax = VProd (region) * density * 1.1;
	printf("Max number of atoms is %d \n", nMolMax);
}

void SetupJob ()
{
	AllocArrays ();
	SetMass ();
	/*
	InitVels ();
	InitAccels ();
	AccumProps (0);
	*/
}

void SetMass ()
{
	int n;

	mass[0] = massIn;
	for (n = 0; n < nMolMax; n++) mol[n].atomType = 1;
}

#define IADD 453806245
#define IMUL 314159269
#define MASK 2147483647
#define SCALE 0.4656612873e-9

//int randSeedP = 3571;
//int randSeedP = 17;

real RandR ()
{
	randSeedP = (randSeedP * IMUL + IADD) & MASK;
	return (randSeedP * SCALE);
}

void CreateNuclei ()
{
	int n, error, k;
	real rrCutGrain, rr;
	VecR dr;

	rrCutGrain = Sqr (rCutGrain);
	n = 0;
	while ( n < grainNumber) {
		VSCopy (nuclei[n].r, RandR () - 0.5, region);
		//VSCopy (nuclei[n].r, 0., region);
		error = 0; // No error
		k = 0;
		while ( (k < n) && ( error == 0)) {
			VSub (dr, nuclei[n].r, nuclei[k].r);
			VWrapAll (dr);
			rr = VLenSq (dr);
			if ( rr < rrCutGrain) error = 1; //Error found
			++k;
		}
		if (error == 0) ++n; 
	}
	for (n = 0; n < grainNumber; n++) {
		// http://www.cognitive-antics.net/mw/index.php?title=Uniform_random_orientation
		VSet (nuclei[n].angles, 2. * Pi * RandR () , asin (RandR ()), Pi * (2. * RandR () - 1. ));
		//VSet (nuclei[n].angles, 0.0,  0.0, 0.0 );
		//printf("nuclei %d \t pos %f %f %f \n", n, WriteVec (nuclei[n].r));
		//printf("nuclei %d \t angle %f %f %f \n", n, WriteVec (nuclei[n].angles));
	}
}

void CreateCrystal ()
{
	int n, nx, ny, nz, i, j, grainMin, k, flag;
	VecR rlat, rrot, dr, rotPrimVect[3], rotBasisVect[nBasis], rotRegion, rRegion;
	VecI halfGrainCells;
	real dist2, min, rotX[3][3], rotY[3][3], rotZ[3][3], rlat2, rlatCutOff, rlatCutOff2;

	n = 0;
	for (i = 0; i < grainNumber; i++) {
		// Rotate primitive and basis vectors
		RotMatX (rotX, nuclei[i].angles.x);
		RotMatY (rotY, nuclei[i].angles.y);
		RotMatZ (rotZ, nuclei[i].angles.z);
		for (k = 0; k < 3; k++) {
				// Matriz de rotación
				MatMul (rotPrimVect[k], primVect[k], rotX);
				rrot = rotPrimVect[k];
				MatMul (rotPrimVect[k], rrot, rotY);
				rrot = rotPrimVect[k];
				MatMul (rotPrimVect[k], rrot, rotZ);		
		}
		for (k = 0; k < nBasis; k++) {
				// Matriz de rotación
				MatMul (rotBasisVect[k], basisVect[k], rotX);
				rrot = rotBasisVect[k];
				MatMul (rotBasisVect[k], rrot, rotY);
				rrot = rotBasisVect[k];
				MatMul (rotBasisVect[k], rrot, rotZ);		
		}
		// Rotate region vectors
		// MatMul (rotRegion, grainRegion, rotX);
		// rrot = rotRegion;
		// MatMul (rotRegion, rrot, rotY);
		// rrot = rotRegion;
		// MatMul (rotRegion, rrot, rotZ);	
		// Initialize atomsPerGrain
		atomsPerGrain[i] = 0;
		// Start lattice loop
		DO_LATTICE (grainCells) {
			for (k = 0; k < nBasis; k++) {
				// Lattice vector
				VLat (rlat, nx, ny, nz, rotPrimVect);
				// Add basis vector
				VVAdd (rlat, rotBasisVect[k]);
				// Substract half region
				VSAdd2 (halfGrainCells, -0.5, grainCells);
				VLat (rRegion, halfGrainCells.x,
					halfGrainCells.y, halfGrainCells.z , rotPrimVect);
				VVAdd (rlat,rRegion);
				/* Porción vieja
				RotationMatrix (rot, nuclei[i].angles.x, nuclei[i].angles.y, nuclei[i].angles.z);
				// Rotate
				MatMul (rrot, rlat, rot);
				*/
				// Check if atom is within a cutoff radius from the grain
				// This will generate spherical crystals :)
				rlat2 = VLenSq (rlat);
				rlatCutOff = (avgGrainSize * grainCellsFactor / 2.);
				rlatCutOff2 = Sqr (rlatCutOff);
				if ( rlat2 < rlatCutOff2) {

				// Translate to grain center
				VVAdd (rlat, nuclei[i].r);
				// Set minimum distance in grain 0
				VSub (dr, rlat, nuclei[0].r);
				VWrapAll (dr);
				min = VLenSq (dr);
				grainMin = 0;
				for (j = 1; j < grainNumber; j++) {
					VSub (dr, rlat, nuclei[j].r);
					VWrapAll (dr);
					dist2 = VLenSq (dr);
					if ( min >= dist2 ) {
						min = dist2;
						grainMin = j;
					}
				}
				if ( grainMin == i) {
					// Apply Periodic Boundary Conditions
					VWrapAll (rlat);
					mol[n].r = rlat;
					mol[n].grain = i;
					++atomsPerGrain[i];
					if ( n >= nMolMax) {
						ErrExit (ERR_TOO_MANY_MOL);
					} else {
					++ n;
					}
				}
				}
			}
		END_LATTICE 
		}	
	}
	nMol = n;
	printf("Total atoms prior to overlap check is %d \n", nMol);
}


void EraseOverlap ()
{
	VecR dr, invWid, rs, shift;
	VecI cc, m1v, m2v, vOff[] = OFFSET_VALS;
	real rr, rrCut;
	int c, j1, j2, m1, m1x, m1y, m1z, m2, n, offset, numErase, *erase, counter;

	AllocMem (erase, nMol, int);
	AllocMem (cellList, VProd (cells) + nMol, int);
	rrCut = Sqr (rCut);
	VDiv (invWid, cells, region);
	for (n = nMol; n < nMol + VProd (cells); n++) cellList[n]=-1;
	DO_MOL {
		VSAdd (rs, mol[n].r, 0.5, region);
		VMul (cc, rs, invWid);
		c = VLinear (cc, cells) + nMol;
		cellList[n] = cellList[c];
		cellList[c] = n;
	}
	DO_MOL erase[n] = 0;
	for (m1z = 0; m1z < cells.z; m1z ++) {
		for (m1y = 0; m1y < cells.y; m1y ++) {
			for (m1x = 0; m1x < cells.x; m1x ++) {
				VSet (m1v, m1x, m1y, m1z);
				m1 = VLinear (m1v, cells) + nMol;
				for (offset = 0; offset < N_OFFSET; offset ++) {
					VAdd (m2v, m1v, vOff[offset]);
					VZero (shift);
					VCellWrapAll ();
					m2 = VLinear (m2v, cells) + nMol;
					DO_CELL (j1, m1) {
						DO_CELL (j2, m2) {
							if (m1 != m2 || j2 < j1) {
								if ( erase[j1] == 0 && erase[j2] == 0) {
									VSub (dr, mol[j1].r, mol[j2].r);
									VVSub (dr, shift);
									rr = VLenSq (dr);
									if (rr < rrCut) {
										if (RandR () < 0.5) erase[j1] = 1;
										else erase[j2] = 1;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	counter = 0;
	for (n = 0; n < nMol; n++) {
		if (erase[n] == 0) {
			mol[counter].r = mol[n].r;
			mol[counter].atomType = mol[n].atomType;
			mol[counter].grain = mol[n].grain;
			++ counter;
		} else if (erase[n] == 1) {
			--atomsPerGrain[mol[n].grain];
		}
	}
	printf("Porosity fraction is %f \n", 1.-(real) counter / ( VProd (region) * density) );
	nMol = counter;
	printf("Total atoms is %d \n", nMol);
}


/*
void SingleStep () 
{
	++ stepCount;
	timeNow = stepCount * deltaT;
	LeapfrogStep (1);
	ApplyBoundaryCond ();
	if (nebrNow) {
		nebrNow = 0;
		dispHi = 0.;
		BuildNebrList ();
	}
	ComputeForces ();
	LeapfrogStep (2);
	EvalProps ();
	AccumProps (1);
	if (stepCount % stepAvg == 0) {
		TakeSnapshot (pos);
		AccumProps (2);
		PrintSummary (stdout);
		AccumProps (0);
	}
}
*/

int GetNameList (int argc, char **argv)
{
	int id, j, k, match, ok;
	char buff[80], *token, str[20];
	float val;
	FILE *fp;

	strcpy (buff, argv[1]);
	//strcat (buff,".in"); // Erased. It was a bane.
	printf ("%s \n",buff);
	if ((fp = fopen (buff, "r")) == 0) return (0);
	// Status 0 : Initialized to no data
	for (k = 0; k < sizeof (nameList) / sizeof (NameList); k++) 
		nameList[k].vStatus = 0;
	ok = 1;
	while (1) {
		fgets (buff, 80, fp);
	       	// If EOF exit while loop
		if ( feof (fp) ) break;
		token = strtok (buff, " \t\n");
	       	// If token is void exit while loop
		if (! token) break;
	       	// Initializes match
		match = 0;
		for (k = 0; k < sizeof (nameList) / sizeof (NameList); k++) {
			if ( strcmp (token,nameList[k].vName) == 0 ) {
			       	// A coincidence was found between nameList and the file
			       	match = 1;
				// If the coincidence has not yet been found proceed
				if (nameList[k].vStatus == 0 ) {
				       	// 1 if everything is OK
					nameList[k].vStatus = 1;
					// Reads every element of array
					for (j = 0; j < nameList[k].vLen; j++) {
						token = strtok (NULL, " \t\n");
						if (token) {
							switch (nameList[k].vType) {
							case N_I:
								*NP_I = atol (token);
								break;
							case N_R:
								*NP_R = atof (token);
								break;
							case N_C:
								*NP_C = *token;
								//printf("%s \n", token);
								//printf("%s",  NP_C);
								break;
							}
						} else {
				       			// 2 if missing data
							nameList[k].vStatus = 2;
							ok=0;
						}
					}
				}
			}
		}
		if (! match) ok = 0;
	}
	fclose (fp);
	for (k = 0; k < sizeof (nameList) / sizeof (NameList); k++) {
		if (nameList[k].vStatus != 1) ok = 0;
	}
	return (ok);
}
	

#define Bar_Length 70

void PrintNameList (FILE *fp)
{
	int j, k;

	for (k = 0; k < Bar_Length; k++) {
	       	fprintf (fp, "-");
       	}
	fprintf (fp, "\n");
	fprintf (fp, "NameList - data \n");	
	for (k = 0; k < sizeof (nameList) / sizeof (NameList) ; k++) {
		fprintf (fp, "%s \t", nameList[k].vName);
		for (j = 0; j < nameList[k].vLen; j++) {
			switch (nameList[k].vType) {
				case N_I:
					fprintf (fp, "%i ", *NP_I);
					break;
				case N_R:
					fprintf (fp, "%lf ", *NP_R);
					break;
				case N_C:
					fprintf (fp, "%c", *NP_C);
					break;
		}
		}
		switch (nameList[k].vStatus) {
			case 0:
				fprintf (fp, "** no data");
				break;
			case 1:
				break;
			case 2:
				fprintf (fp, "** missing data");
				break;
		}
		fprintf (fp, "\n");
       	}
	for (k = 0; k < Bar_Length; k++) {
		fprintf (fp, "-");
       	}	
	fprintf (fp, "\n");
}


void AllocArrays ()
{
	AllocMem (mol, nMolMax, Mol);
	AllocMem (mass, numAtomType, real);
	AllocMem (nuclei, grainNumber, Nuclei);
	AllocMem (atomsPerGrain, grainNumber, int);
}

void VRand (VecR *p)
{
	real s, x, y, z;

	s=2.;
	while (s > 1.) {
		x = 2. * RandR () - 1.;
		y = 2. * RandR () - 1.;
		s = Sqr (x) + Sqr (y);
	}
	p->z = 1. - 2. * s;
	s = 2. * sqrt (1. - s);
	p->x = s * x;
	p->y = s * y;
}

/*
void PrintSummary (FILE *fp)
{
	fprintf (fp, 
	"%5d %8.4f %9.6f %9.6f %9.6f %9.7f %9.6f\n",
	stepCount, timeNow, VCSum (vSum) / nMol, totEnergy.sum,
	kinEnergy.sum, potEnergy.sum, pressure.sum);
}
*/

void WriteData (FILE *fp)
{
	int n;

	fprintf ( fp, "LAMMPS data file by Pablo Piaggi, timestep = 0 \n\n");
	fprintf ( fp, "%d atoms\n", nMol);
	fprintf ( fp, "1 atom types\n\n");
	fprintf ( fp, "%17.16e %17.16e xlo xhi\n", -region.x / 2., region.x / 2.);
	fprintf ( fp, "%17.16e %17.16e ylo yhi\n", -region.y / 2., region.y / 2.);
	fprintf ( fp, "%17.16e %17.16e zlo zhi\n\n", -region.z / 2., region.z / 2.);
	fprintf ( fp, "Masses\n\n");
	for (n = 0; n < numAtomType; n++) fprintf ( fp,"%d %f\n\n", n + 1, mass[n]);
	fprintf ( fp, "Atoms\n\n");
	DO_MOL fprintf (fp, "%d %d %17.16e %17.16e %17.16e %d %d %d\n", n+1, mol[n].atomType,
	WriteVec (mol[n].r), 0, 0, 0);	
	fprintf ( fp, "\n");
	fprintf ( fp, "Velocities\n\n");
	DO_MOL fprintf (fp, "%d %17.16e %17.16e %17.16e\n", n+1, WriteVec (mol[n].rv));	
	fprintf ( fp, "\n");
}

void ErrExit (int code)
{
	printf ("Error: %s\n", errorMsg[code]);
	exit (0);
}

void TakeSnapshot (FILE *fp)
{
	int n;

	fprintf (fp," %d \n \n", nMol);
	DO_MOL fprintf (fp," %17.16e %17.16e %17.16e \n", mol[n].r.x, 
	mol[n].r.y, mol[n].r.z);
}

void WriteGrainDistrib (FILE *fp)
{
	real grainDiameter, dmin, dmax, d1, d2;
	int n, i, steps = rangeGD / stepGD, grains[steps];
	
	dmin = avgGrainSize-rangeGD/2.;	
	dmax = avgGrainSize+rangeGD/2.;	
	for (i = 0; i < steps; i++) {
		d1 = dmin + i * (dmax-dmin) / (real) steps;
		d2 = dmin + (i+1) * (dmax-dmin) / (real) steps;
		grains[i] = 0;
		for (n = 0; n < grainNumber; n++) {	
			grainDiameter = pow (6. * atomsPerGrain[n] / ( Pi * density ), 1. / 3.);
			if (grainDiameter > d1 && grainDiameter < d2) ++grains[i];
		}
		fprintf (fp,"%f %f \n", (d1 + d2) / 2., (real) grains[i] / ( (real) grainNumber * stepGD));
	}
}       	


