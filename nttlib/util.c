/* ------------------------- Some utilitities -----------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <time.h>

#include "./nttcommon.h"

gsl_rng * GCMCr;
int randomSeed;

void initRand(void) {

    // random number generator initialization
    if (randomSeed < 0)
	// use the time function to seed the random number generator
	randomSeed = time(NULL);

    printf ("The random seed is %d \n", randomSeed );

    GCMCr = gsl_rng_alloc (gsl_rng_mt19937);
    gsl_rng_set (GCMCr, (unsigned int)randomSeed);
    
    //printf ("Random integer range %lu %lu.\n", gsl_rng_min(GCMCr), gsl_rng_max(GCMCr) );

    /* we run about 10^6 random number to avoid initial non-randomness */
    int i;
    for (i=0; i<1000000; i++)  arand();
    
    printf ("The first random number is %.4f\n", arand() );

}

void rand_unitvector(DOUBLE *x, DOUBLE *y, DOUBLE *z)
/* creat a random vector on the surface of an unit sphere. Using
 * Marsaglia algorithm, Allen & Tildesley, p. 349, G.4		*/
{
	DOUBLE zeta1, zeta2, r2;
	
	while(1) {
		zeta1	= 1.0 - arand()*2;
		zeta2	= 1.0 - arand()*2;
		
		r2	= zeta1*zeta1+zeta2*zeta2;
		if (r2 < 1.0)	break;
	}
	
	*x	= zeta1*2*sqrt(1.0-r2);
	*y	= zeta2*2*sqrt(1.0-r2);
	*z	= 1.0-r2-r2;
}

/* search the string STR for ORIG and replace it by REP */
char *replace_str(char *dest, const char *str, const char *orig, const char *rep) {
	char *p;

	/* p holds the the start of orig in str */
	if(!(p = strstr(str, orig)))  // Is 'orig' even in 'str'? 
	{
		return NULL;
	}

// 	char *buffer;
// 	size_t length = strlen(str) + strlen(rep) - strlen(orig) + 16;
// 	buffer = (char *)malloc(length);

	strncpy(dest, str, p-str); // Copy characters from 'str' start to 'orig' start
	dest[p-str] = '\0';

	sprintf(dest+(p-str), "%s%s", rep, p+strlen(orig));

	return dest;
}


void omp_manual_static_schedule_getlimits(int threadnum, int threadtotal, \
    int N, int offset, int *startpos, int *endpos)
/* Used for manually setting up the static schedule for "for" loop in omp parallel.
 * give the total iterations, threadnum, returns the starting
 * and ending positions of the thread's iterations.
 * 
 * E.g., 3 threads, 7 iterations
 * 	thread 0 does 0,1,2
 * 	thread 1 does 3,4
 * 	thread 2 does 5,6
 */
{
    int min_iter = N / threadtotal;
    int less = N % threadtotal;
   
    if (threadnum < less) {
      *startpos = (min_iter+1)*threadnum + offset;
      *endpos = (min_iter+1)*(threadnum+1) + offset;
    } else {
      *startpos = less+min_iter*threadnum + offset;
      *endpos = less+min_iter*(threadnum+1) + offset;
    }
    return;
}

void omp_manual_static_schedule_getlimits_even(int threadnum, int threadtotal, \
    int N, int offset, int *startpos, int *endpos)
/* Used for manually setting up the static schedule for "for" loop in omp parallel.
 * give the total iterations, threadnum, returns the starting
 * and ending positions of the thread's iterations.
 * 
 * Same as the one above but we want startpos to be even (for 2xdata alignment)
 */
{
  
/* some threads get 2x+2 iterations, some get 2x.
      m*(2x+2) + n*2x = N ( if N is even, N-1 if odd )
      m + n = Nt (threadtotal)
      0 <= m < Nt
      
      2x * Nt + 2*m = N < (2x+2) * Nt 
      
 */
    int x = N / (2*threadtotal);
    int m = ( N % (2*threadtotal) ) / 2;

    if (threadnum < m) {
	*startpos = (2*x+2)*threadnum + offset;
	*endpos = *startpos +2*x+2;
	return;
    } 
    
    *startpos = offset + m*(2*x+2) + (threadnum-m)*2*x;
    *endpos = *startpos + 2*x;
    
    if (threadnum == threadtotal-1)
	*endpos = N + offset;

}


void omp_manual_folding_schedule_getlimits(int threadnum, int threadtotal, \
	int N, int offset, int limits[4])
/* Used for manually setting up the folding schedule for the "for" loop in omp parallel.
 * give the total iterations, threadnum, returns the starting
 * and ending positions of the thread's iterations.
 * 
 * E.g., 2 threads, 7 iterations
 * 	loop 0 does 0,1,6
 * 	loop 1 does 2,3,4,5
 */
{
    int doubleNthreads = 2 * threadtotal;

    /* we do this by divide N to 2*threadtotal bins.
     * threadnum will process bin [threadnum] and bin [2*threadtotal-1-threadnum]
     */
    omp_manual_static_schedule_getlimits(threadnum, doubleNthreads, \
	N, offset, limits, limits+1);
    omp_manual_static_schedule_getlimits(doubleNthreads-threadnum-1, doubleNthreads, \
	N, offset, limits+2, limits+3);
    
//     if (threadnum == threadtotal-1) {
// 	limits[1] = limits[3];
// 	limits[2] = limits[3] = 0;
//     }
    
    return;
}


// #define COMMENT_CHAR	'#'
// 
// int read_a_line(FILE *fp, char *s, int maxlength)
// /* Read a line from fp to string, drop comment line, and blank lines */
// {
// 	int i, finished;
// 	char *sp;
// 
// 	finished = FALSE;
// 	while(finished == FALSE) {
// 		sp = fgets(s, maxlength,fp);
// 		if(sp==NULL) return FALSE;
// 		if(s[0]!=COMMENT_CHAR) {
// 			for(i=0; i<strlen(s);i++)
// 				if(isalnum(s[i])) finished=TRUE;
// 		}
// 	}
// 	for(i=0; i<strlen(s);i++)
// 		if(s[i]==COMMENT_CHAR) s[i]='\0';
// 	
// 	return TRUE;
// }

void error(const char *s)
{
	fprintf(stderr,"%s\n",s);
}

#ifdef NTT_OPENCL

size_t PrintDeviceInfo(cl_device_id device) {

	char queryBuffer[1024];
	int queryInt;
	size_t queryDim[3] = {0,0,0};
	// cl_int clError;
	// clError = 
	clGetDeviceInfo(device, CL_DEVICE_NAME, sizeof(queryBuffer), &queryBuffer, NULL);
	printf("CL_DEVICE_NAME: %s\n", queryBuffer);
	queryBuffer[0] = '\0';
	//clError = 
	clGetDeviceInfo(device, CL_DEVICE_VENDOR, sizeof(queryBuffer), &queryBuffer, NULL);
	printf("CL_DEVICE_VENDOR: %s\n", queryBuffer);
	queryBuffer[0] = '\0';
	//clError = 
	clGetDeviceInfo(device, CL_DRIVER_VERSION, sizeof(queryBuffer), &queryBuffer,	NULL);
	printf("CL_DRIVER_VERSION: %s\n", queryBuffer);
	queryBuffer[0] = '\0';
	//clError = 
	clGetDeviceInfo(device, CL_DEVICE_VERSION, sizeof(queryBuffer), &queryBuffer,	NULL);
	printf("CL_DEVICE_VERSION: %s\n", queryBuffer);
	queryBuffer[0] = '\0';
	//clError = 
	clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(int), &queryInt, NULL);
	printf("CL_DEVICE_MAX_COMPUTE_UNITS: %d\n", queryInt);
	//clError = 
	clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), queryDim, NULL);
	printf("CL_DEVICE_MAX_WORK_GROUP_SIZE: %lu\n", queryDim[0]);
	//clError = 
	clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_ITEM_SIZES, 3*sizeof(size_t), queryDim, NULL);
	printf("CL_DEVICE_MAX_WORK_ITEM_SIZES: (%lu %lu %lu)\n", queryDim[0], queryDim[1], queryDim[2]);
	
	clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(size_t), queryDim+1, NULL);
	printf("CL_DEVICE_LOCAL_MEM_SIZE: %lu\n", queryDim[1]);

	return queryDim[0];

}

const char * get_CL_error_string(cl_int err) {

         switch(err){
             case 0: return "CL_SUCCESS";
             case -1: return "CL_DEVICE_NOT_FOUND";
             case -2: return "CL_DEVICE_NOT_AVAILABLE";
             case -3: return "CL_COMPILER_NOT_AVAILABLE";
             case -4: return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
             case -5: return "CL_OUT_OF_RESOURCES";
             case -6: return "CL_OUT_OF_HOST_MEMORY";
             case -7: return "CL_PROFILING_INFO_NOT_AVAILABLE";
             case -8: return "CL_MEM_COPY_OVERLAP";
             case -9: return "CL_IMAGE_FORMAT_MISMATCH";
             case -10: return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
             case -11: return "CL_BUILD_PROGRAM_FAILURE";
             case -12: return "CL_MAP_FAILURE";

             case -30: return "CL_INVALID_VALUE";
             case -31: return "CL_INVALID_DEVICE_TYPE";
             case -32: return "CL_INVALID_PLATFORM";
             case -33: return "CL_INVALID_DEVICE";
             case -34: return "CL_INVALID_CONTEXT";
             case -35: return "CL_INVALID_QUEUE_PROPERTIES";
             case -36: return "CL_INVALID_COMMAND_QUEUE";
             case -37: return "CL_INVALID_HOST_PTR";
             case -38: return "CL_INVALID_MEM_OBJECT";
             case -39: return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
             case -40: return "CL_INVALID_IMAGE_SIZE";
             case -41: return "CL_INVALID_SAMPLER";
             case -42: return "CL_INVALID_BINARY";
             case -43: return "CL_INVALID_BUILD_OPTIONS";
             case -44: return "CL_INVALID_PROGRAM";
             case -45: return "CL_INVALID_PROGRAM_EXECUTABLE";
             case -46: return "CL_INVALID_KERNEL_NAME";
             case -47: return "CL_INVALID_KERNEL_DEFINITION";
             case -48: return "CL_INVALID_KERNEL";
             case -49: return "CL_INVALID_ARG_INDEX";
             case -50: return "CL_INVALID_ARG_VALUE";
             case -51: return "CL_INVALID_ARG_SIZE";
             case -52: return "CL_INVALID_KERNEL_ARGS";
             case -53: return "CL_INVALID_WORK_DIMENSION";
             case -54: return "CL_INVALID_WORK_GROUP_SIZE";
             case -55: return "CL_INVALID_WORK_ITEM_SIZE";
             case -56: return "CL_INVALID_GLOBAL_OFFSET";
             case -57: return "CL_INVALID_EVENT_WAIT_LIST";
             case -58: return "CL_INVALID_EVENT";
             case -59: return "CL_INVALID_OPERATION";
             case -60: return "CL_INVALID_GL_OBJECT";
             case -61: return "CL_INVALID_BUFFER_SIZE";
             case -62: return "CL_INVALID_MIP_LEVEL";
             case -63: return "CL_INVALID_GLOBAL_WORK_SIZE";
             default: return "Unknown OpenCL error";
         }
}

#endif