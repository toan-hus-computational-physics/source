
// This is a header file for common usefull functions, macros used in 
// various simulation

#ifndef NTTCOMMON_H
#define NTTCOMMON_H

#include <ctype.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "amdlibm.h"

#include <x86intrin.h>
#include <sys/time.h>

// for cl_double3 vector types
#include "CL/cl.h"

/* use GSL for random number generator */
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

extern int randomSeed;
extern gsl_rng * GCMCr;
void initRand(void);
/**
 * A random number between [0,1) 
 */
#define arand()  gsl_rng_uniform (GCMCr)
/** 
 * a random integer between [0, N-1] 
 */
#define arandN(N) ( (int) gsl_rng_uniform_int (GCMCr, N) )
/** 
 * Creat a random vector on the surface of an unit sphere. Using
 * Marsaglia algorithm, Allen & Tildesley, p. 349, G.4 
 */
void rand_unitvector(double *x, double *y, double *z);

/** 
 * Search the string STR for ORIG and replace it by REP 
 * DEST holds the new str.
 * return DEST if success, otherwise (no ORIG in STR) returns NULL
 */
char *replace_str(char *dest, const char *str, const char *orig, const char *rep);


typedef unsigned long long timestamp_t;

/**
 * Get current timestamp as a unsigned long long
 */
static inline timestamp_t get_timestamp () {
	struct timeval now;
	gettimeofday (&now, NULL);
	return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}

//------------------------------ common abbreviations ----------------------------

/* 2^(-1/6) */
#define TWO_1_6	0.89089871814
#define TWO_1_3 0.79370052598

#define	DOUBLE_PRECISION

typedef double DOUBLE ;
#define zero	0.0
#define half	0.5
#define one	1.0
#define two	2.0
//#define Infinity	1e300

#define TRUE	1
#define FALSE	0

void vectormathinfo();

//----------------------------- short inline functions --------------------------

static inline void triangleStore ( const double value, double *E, const int N, 
    const int i, const int j) {
    /* store the value into the matrix element E[i,j] using 
     the upper triangle storage */
    if ( i<j ) 	*(E + N*i + j) = value;
    else	*(E + N*j + i) = value;
}

static inline void triangleAccess ( double *value, const double *E, const int N, 
    const int i, const int j) {
    /* store the value into the matrix element E[i,j] using 
     the upper triangle storage */
    if ( i<j ) 	*value = E[N*i+j];
    else	*value = E[N*j+i];
}

//static inline __attribute__ ((__gnu_inline__, __always_inline__,const))
static inline double sqr(const double x) {
    return x*x;
}

//static inline __attribute__ ((__gnu_inline__, __always_inline__,const))
static inline double cubic(const double x) {
    return x*x*x;
}

/**
 * PBC stands for Periodic Boundary Condition
 */
static inline double PBC_position(const double x, const double L) {
	/* put x in the periodic box from 0 to L */
	if ( x < 0 ) return x+L;
	if ( x >= L ) return x-L;
	return x;
}

static inline double PBC_distance(const double x, const double L_2) {
	/* put x in the periodic box from -L_2 to L_2 */
	if ( x < -L_2 ) return x+2*L_2;
	if ( x >= L_2 ) return x-2*L_2;
	return x;
}

static inline void PBC3_position(cl_double3 *x, const cl_double3 L) {
	x->x = PBC_position(x->x, L.x);
	x->y = PBC_position(x->y, L.y);
	x->z = PBC_position(x->z, L.z);
}

static inline void PBC3_distance(cl_double3 *x, const cl_double3 L_2) {
	x->x = PBC_distance(x->x, L_2.x);
	x->y = PBC_distance(x->y, L_2.y);
	x->z = PBC_distance(x->z, L_2.z);
}

/**
 * 3D Vector operations
 */
static inline void vecScale(const double alpha, cl_double3 *a) {
    a->x = alpha * a->x;
    a->y = alpha * a->y;
    a->z = alpha * a->z;
}

static inline void vecScaleTo(const double alpha, const cl_double3 a, cl_double3 *b) {
    b->x = alpha * a.x;
    b->y = alpha * a.y;
    b->z = alpha * a.z;
}

static inline void vecScaleAdd(const double alpha, const cl_double3 a, cl_double3 *b) {
    b->x += alpha * a.x;
    b->y += alpha * a.y;
    b->z += alpha * a.z;
}

static inline void vecSub (const cl_double3 a, const cl_double3 b, cl_double3 *c) {
    c->x = a.x - b.x;
    c->y = a.y - b.y;
    c->z = a.z - b.z;
}

static inline void vecAdd (const cl_double3 a, const cl_double3 b, cl_double3 *c) {
    c->x = a.x + b.x;
    c->y = a.y + b.y;
    c->z = a.z + b.z;
}

static inline void vecAddTo (const cl_double3 a, cl_double3 *b) {
    b->x += a.x;
    b->y += a.y;
    b->z += a.z;
}

static inline void vecSubTo (const cl_double3 a, cl_double3 *b) {
    b->x -= a.x;
    b->y -= a.y;
    b->z -= a.z;
}

static inline void crossProduct (const cl_double3 a, const cl_double3 b, cl_double3 *c) {
    c->x = a.y*b.z - a.z*b.y;
    c->y = a.z*b.x - a.x*b.z;
    c->z = a.x*b.y - a.y*b.x;
}

static inline double dotProduct(const cl_double3 a, const cl_double3 b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

static inline double distance2(const cl_double3 r0, const cl_double3 r1) {
    return sqr(r1.x-r0.x)+sqr(r1.y-r0.y)+sqr(r1.z-r0.z);
}

static inline double signum (double x) {
    if (x > 0) return 1.0;
    if (x < 0) return -1.0;
    return 0.0;
}

/* torsional angle */
static inline double dihedralAngle(const cl_double3 r0, const cl_double3 r1, const cl_double3 r2, const cl_double3 r3) {
	/* calculate the dihedral angle between the plane
	 * formed by the points (r0,r1,r2) and (r1,r2,r3) */
	/* use the formular found in wikipedia.org */

	cl_double3 b1, b2, b3, b12, b23, b123;
	
	vecSub (r2, r1, &b2);
	vecSub (r1, r0, &b1);	
	vecSub (r3, r2, &b3);
	
	crossProduct(b1, b2, &b12);
	crossProduct(b2, b3, &b23);
	crossProduct(b12, b23, &b123);

	return  atan2( dotProduct(b123, b2)/sqrt(dotProduct(b2,b2)), dotProduct(b12, b23) );

}

// static inline double dihedralAngle(const cl_double3 r0, const cl_double3 r1, const cl_double3 r2, const cl_double3 r3) {
// 	/* calculate the dihedral angle between the plane
// 	 * formed by the points (r0,r1,r2) and (r1,r2,r3) */
// 	/* use the formular found in Bekker's thesis.
// 	 * Seems that this and the previous code gives the same result.
// 	 */
// 
// 	double r01[3], r21[3], r32[3], R[3], S[3],n[3];
// 	
// 	vecSub (r0, r1, r01);
// 	vecSub (r2, r1, r21);
// 	vecSub (r3, r2, r32);
// 	
// 	double r21_2 = dotProduct(r21,r21);
// 	R.x = r01.x - dotProduct(r01,r21)*r21.x/r21_2;
// 	R.y = r01.y - dotProduct(r01,r21)*r21.y/r21_2;
// 	R.z = r01.z - dotProduct(r01,r21)*r21.z/r21_2;
// 	S.x = r32.x - dotProduct(r32,r21)*r21.x/r21_2;
// 	S.y = r32.y - dotProduct(r32,r21)*r21.y/r21_2;
// 	S.z = r32.z - dotProduct(r32,r21)*r21.z/r21_2;
// 
// 	double theta = dotProduct(R,S)/sqrt( dotProduct(R,R)*dotProduct(S,S) );
// 	if (theta > 1.0) theta = 1.0;
// 	if (theta < -1.0) theta = -1.0;
// 	theta = acos(theta);
// 	crossProduct(r21, r32, n);
// 	r21_2 = -dotProduct(r01, n);
// 	if (r21_2 < 0) return -theta;
// 	if (r21_2 > 0) return theta;
// 	return 0;
// 
// }

static inline double cosThetaBend ( const cl_double3 r0, const cl_double3 r1, const cl_double3 r2 ) {
	/* calculate the cosine of the angle made of r1->r0 and r1-> r2 */
	cl_double3 a, b;
	vecSub (r0, r1, &a);
	vecSub (r2, r1, &b);

	double res = dotProduct(a, b) / sqrt( dotProduct(a,a)*dotProduct(b,b) );
	if ( fabs(res) <= 1 ) return res;
	if (res > 0) return 1;
	
	return -1;
}

static inline double cosTheta4 ( const cl_double3 r0, const cl_double3 r1, const cl_double3 r2, const cl_double3 r3 ) {
	/* calculate the cosine of the angle made of r10 and r23 */
	cl_double3 a, b;
	vecSub (r0, r1, &a);
	vecSub (r3, r2, &b);
	return dotProduct(a, b) / sqrt( dotProduct(a,a)*dotProduct(b,b) );
}

double erfc_positive_arg(double x) __attribute__((const));

#define xDOTy(x,y) ( (x).x*(y).x+(x).y*(y).y+(x).z*(y).z )
// extern inline __attribute__ ((__gnu_inline__, __always_inline__))
// double xDOTy13(const cl_double3 x, const double x2, const double y2, const double z2)
// {
//     return x.x*x2+x.y*y2+x.z*z2;
// }
#define xDOTy13(x,x2,y2,z2) ( (x).x*(x2)+(x).y*(y2)+(x).z*(z2) )


/* =================================
 * 	Other utility functions 
 * ================================= */

/* Read a line from fp to string, drop comment line, and blank lines */
int read_a_line(FILE *fp, char *s, int maxlength);

void error(const char *s);

/* Memcpy(dest, src, length) copy length items of 8byte each 
 * (double precision numbers) from src to dest */

#define Memcpy(dest,src,length)	memcpy((dest),(src),(length)*sizeof(double))

#define cossin(a,c,s)	sincos( a,&(s),&(c) )

//------------ OpenMP scheduler functions -------------------

void omp_manual_static_schedule_getlimits(int t, int tmax, \
	int N, int offset, int *startpos, int *endpos);
void omp_manual_static_schedule_getlimits_even(int threadnum, int threadtotal, \
    int N, int offset, int *startpos, int *endpos);
void omp_manual_folding_schedule_getlimits(int threadnum, int threadtotal, \
	int N, int offset, int limits[4]);

//--------------------------   common memory routines   ---------------------------

/* using the following directive to replace all compiler math functions
with AMD LibM functions */

// #ifdef USE_ACML
// 
// #include <acml.h>
//     #include <acml_mv.h>
//     #include <omp.h>
// 
//     #define mydcopy(n,x,dx,y,dy)	{ dcopy_( &(n), (x), &(dx), (y), &(dy) ); }
//     #define myzfft3d(mode, l,m,n, X,comm,info) { zfft3d_(&(mode),&(l),&(m),&(n),(X),(comm),&(info)); }
//     #define myddot(len,x,dx,y,dy)	ddot_(&(len),(x),&(dx),(y),&(dy))
// 
// #elif defined(USE_MKL)
//     #include <math.h>
//     #include <omp.h>
//     #include <mkl.h>
//     
// 
//     #define mydcopy(n,x,dx,y,dy)	{ dcopy( &(n), (x), &(dx), (y), &(dy)  ); }
// 
// #else
//     #include <math.h>
//     #include <string.h>
//     #define mydcopy(n,x,dx,y,dy)	{ memcpy( (y), (x), (n)*sizeof(double) ); }
// 
// #endif


//void Memcpy(void *dest,void *src, int length);

#ifdef NTT_OPENCL
size_t PrintDeviceInfo(cl_device_id device);
const char * get_CL_error_string(cl_int err);
#endif


#endif
