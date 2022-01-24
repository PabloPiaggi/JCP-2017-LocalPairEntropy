#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef double real;

typedef enum {N_I, N_R, N_C} VType;

typedef struct {
	char *vName;
	void *vPtr;
	VType vType;
	int vLen, vStatus;
} NameList;

typedef struct {
	real x, y, z;
} VecR;

typedef struct {
	int x, y, z;
} VecI;

typedef struct {
	VecR r, rv, ra;
	int atomType;
	int grain;
} Mol;

typedef struct {
	VecR r;
	VecR angles;
} Nuclei;

typedef struct {
	real val, sum, sum2;
} Prop;

#define NameI(x) {#x, &x, N_I, sizeof (x) / sizeof (int)}
#define NameR(x) {#x, &x, N_R, sizeof (x) / sizeof (real)}
#define NameC(x) {#x, &x, N_C, sizeof (x) / sizeof (char)}

#define DO_MOL	for (n = 0; n < nMol; n++)

#define DO_LATTICE(v1)	for (nx = 0; nx < (v1).x; nx++) { \
				for (ny = 0; ny < (v1).y; ny++) { \
					for (nz = 0; nz < (v1).z; nz++)

#define END_LATTICE	} }

#define VProd(v)	((v).x * (v).y * (v).z)

#define VAdd(v1, v2, v3)  \
	(v1).x = (v2).x + (v3).x, \
	(v1).y = (v2).y + (v3).y, \
	(v1).z = (v2).z + (v3).z

#define VVAdd(v1, v2)	VAdd(v1, v1, v2)

#define VScale(v, s) \
	(v).x *= s, \
	(v).y *= s, \
	(v).z *= s

#define VSub(v1, v2, v3)  \
	(v1).x = (v2).x - (v3).x, \
	(v1).y = (v2).y - (v3).y, \
	(v1).z = (v2).z - (v3).z

#define VDot(v1, v2) \
	((v1).x * (v2).x + (v1).y * (v2).y + (v1).z * (v2).z )

#define VSAdd(v1, v2, s3, v3)  \
	(v1).x = (v2).x + (s3) * (v3).x, \
	(v1).y = (v2).y + (s3) * (v3).y, \
	(v1).z = (v2).z + (s3) * (v3).z

#define VSAdd2(v1, s3, v3)  \
	(v1).x = (s3) * (v3).x, \
	(v1).y = (s3) * (v3).y, \
	(v1).z = (s3) * (v3).z

#define VSet(v, sx, sy, sz) \
	(v).x = sx, \
	(v).y = sy, \
	(v).z = sz

#define VSetAll(v, s)	VSet (v, s, s, s)

#define VZero(v)	VSetAll (v, 0)

#define VVSAdd(v1, s2, v2) VSAdd (v1, v1, s2, v2)

#define VLenSq(v)	VDot (v, v)

#define VLen(v)		pow ( VLenSq (v), 1./2.)

#define NP_I ((int *) (nameList[k].vPtr) + j )
#define NP_R ((real *) (nameList[k].vPtr) + j )
#define NP_C ((char *) (nameList[k].vPtr) + j )

#define AllocMem(a, n, t)  a = (t *) malloc ((n) * sizeof (t) )

#define VDivInt(v1, s2) \
	(v1).x = ( (v1).x / (real)s2 ), \
	(v1).y = ( (v1).y / (real)s2 ), \
	(v1).z = ( (v1).z / (real)s2 )

#define PropZero(v) \
	v.sum = 0., \
	v.sum2 = 0.

#define PropAccum(v) \
	v.sum += v.val, \
	v.sum2 += Sqr (v.val)

#define PropAvg(v, n) \
	v.sum /= n, \
	v.sum2 = sqrt (Max (v.sum2 / n - Sqr (v.sum), 0.))

#define PropEst(v) \
	v.sum, v.sum2 

#define Sqr(x) ((x) * (x))

#define Cube(x) ((x) * (x) * (x))

#define Sgn(x, y) (((y) >= 0) ? (x) : (- (x)))

#define IsEven(x) ((x) & ~1)

#define IsOdd(x) ((x) & 1)

#define Nint(x) \
(((x) < 0.) ? (- (int) (0.5 - (x))): ((int) (0.5 + (x))))

#define Min(x1, x2) \
(((x1) < (x2)) ? (x1) : (x2))

#define Max(x1, x2) \
(((x1) > (x2)) ? (x1) : (x2))

#define Min3(x1, x2, x3) \
(((x1) < (x2)) ? (((x1) < (x3)) ? (x1) : (x3)) : \
(((x2) < (x3)) ? (x2) : (x3)))

#define Max3(x1, x2, x3) \
(((x1) > (x2)) ? (((x1) > (x3)) ? (x1) : (x3)) : \
(((x2) > (x3)) ? (x2) : (x3)))

#define Clamp(x, lo, hi) \
(((x) >= (lo) && (x) <= (hi)) ? (x) : (((x) < (lo)) ? \
(lo) : (hi)))

#define VCSum(v1) \
	((v1).x + (v1).y + (v1).z)

#define VWrap(v, t) \
	if (v.t >= 0.5 * region.t) v.t -= region.t; \
	else if (v.t < -0.5 * region.t) v.t += region.t

#define VWrapAll(v) \
	{VWrap (v, x); \
	VWrap (v, y); VWrap (v, z);}

#define VMul(v1, v2, v3) \
	(v1).x = (v2).x * (v3).x, \
	(v1).y = (v2).y * (v3).y, \
	(v1).z = (v2).z * (v3).z

#define VDiv(v1, v2, v3) \
	(v1).x = (v2).x / (v3).x, \
	(v1).y = (v2).y / (v3).y, \
	(v1).z = (v2).z / (v3).z

#define VSCopy(v2, s1, v1) \
	 (v2).x = (s1) * (v1).x, \
	 (v2).y = (s1) * (v1).y, \
	 (v2).z = (s1) * (v1).z

#define VVSub(v1, v2)	VSub (v1, v1, v2)

#define VLinear(p, s) \
	(((p).z * (s).y + (p).y) * (s).x + (p).x)

#define VCellWrap(t) \
	if (m2v.t >= cells.t) { \
		m2v.t = 0; \
		shift.t = region.t; \
	} else if (m2v.t < 0) { \
		m2v.t = cells.t - 1; \
		shift.t = - region.t; \
	}

#define VCellWrapAll() \
	{VCellWrap (x); \
	VCellWrap (y); \
	VCellWrap (z);}

#define DO_CELL(j, m) \
	for (j = cellList[m]; j >= 0; j = cellList[j])

#define N_OFFSET 14

#define OFFSET_VALS \
	{{0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {-1,1,0}, {0,0,1}, \
	{1,0,1}, {1,1,1}, {0,1,1}, {-1,1,1}, {-1,0,1}, \
	{-1,-1,1}, {0,-1,1}, {1,-1,1}}

enum {ERR_NONE, ERR_TOO_MANY_NEBRS, ERR_TOO_MANY_MOL};

char *errorMsg[] = {"", "too many neighbors", "too many atoms"};

#define WriteVec(v1)	(v1).x, (v1).y, (v1).z

#define RotationMatrix(A, alfa, beta, gamma) \
(A)[0][0] = cos (gamma) * cos (alfa) - cos (beta) * sin (alfa) * sin (gamma); \
(A)[0][1] = cos (gamma) * sin (alfa) + cos (beta) * cos (alfa) * sin (gamma); \
(A)[0][2] = sin (beta) * sin (gamma); \
(A)[1][0] = - sin (gamma) * cos (alfa) - cos (beta) * sin (alfa) * cos (gamma); \
(A)[1][1] = - sin (gamma) * sin (alfa) + cos (beta) * cos (alfa) * cos (gamma); \
(A)[1][2] = sin (beta) * cos (gamma); \
(A)[2][0] = sin (beta) * sin (alfa); \
(A)[2][1] = - sin (beta) * cos (alfa); \
(A)[2][2] = cos (beta)

#define RotMatX(A, alpha) \
(A)[0][0] = 1. ; \
(A)[0][1] = 0. ; \
(A)[0][2] = 0. ; \
(A)[1][0] = 0. ; \ 
(A)[1][1] = cos (alpha) ; \
(A)[1][2] = sin (alpha) ; \ 
(A)[2][0] = 0. ; \
(A)[2][1] = - sin (alpha) ; \
(A)[2][2] = cos (alpha)  

#define RotMatY(A, alpha) \
(A)[0][0] = cos (alpha) ; \
(A)[0][1] = 0. ; \
(A)[0][2] = -sin (alpha); \
(A)[1][0] = 0. ; \ 
(A)[1][1] = 1. ; \
(A)[1][2] = 0. ; \ 
(A)[2][0] = sin (alpha) ; \
(A)[2][1] = 0. ; \
(A)[2][2] = cos (alpha) 

#define RotMatZ(A, alpha) \
(A)[0][0] = cos (alpha) ; \
(A)[0][1] = sin (alpha) ; \
(A)[0][2] = 0. ; \
(A)[1][0] = - sin (alpha) ; \ 
(A)[1][1] = cos (alpha) ; \
(A)[1][2] = 0. ; \ 
(A)[2][0] = 0. ; \
(A)[2][1] = 0. ; \
(A)[2][2] = 1. 

#define MatMul(v1,v2,A) \
	(v1).x = A[0][0] * (v2).x + A[0][1] * (v2).y + A[0][2] * (v2).z , \
	(v1).y = A[1][0] * (v2).x + A[1][1] * (v2).y + A[1][2] * (v2).z , \
	(v1).z = A[2][0] * (v2).x + A[2][1] * (v2).y + A[2][2] * (v2).z 

#define Pi 3.14159265359

#define VLat(v1, s1, s2, s3, pv) \
	(v1).x = s1 * (pv[0]).x + s2 * (pv[1]).x + s3 * (pv[2]).x; \
	(v1).y = s1 * (pv[0]).y + s2 * (pv[1]).y + s3 * (pv[2]).y; \
	(v1).z = s1 * (pv[0]).z + s2 * (pv[1]).z + s3 * (pv[2]).z
