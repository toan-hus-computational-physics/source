#include <time.h>
#define S 			1e9	//(nm) dien tich mo phong
double  L=sqrt(S);			//(nm) kich thuoc mang
#define nNuclSite 		100	//so site tren mang L*L

#define nAtomPerParticle 	247	//
#define nParticlePerStep 	100	//
#define n_steps   		2000	//so buoc thoi gian

#define t1 			5.33e5	//(ns.nm2) thoi gian de co 1 atom Ag tren dien tich 1 nm2

double  dt=(t1/S*nAtomPerParticle*nParticlePerStep);	//(ns) time step
#define D 100				//(nm2/ns) do linh dong cua Ag tren Pt
#define a 0.409				//(nm) hang so mang tinh the cua Ag

#define N_mean ( (int) ( 3.16e3 / sqrt(D*dt) * nParticlePerStep ) ) //Average value cua so atom tu do
//#define R_Ag 0.172	//Ban kinh Van der Waals cua Ag: 0.172 nm
//#define V_Ag (4.0/3.0*M_PI*R_Ag*R_Ag*R_Ag) //The tich Van der Waals cua Ag

// cac thong tin gan voi moi atom
typedef struct {
	double pos[2];		//2 toa do cua atom
	int n_catch;		//so luong cac site co the bat atom do
	double catch_idx[nNuclSite];	//index cua cac site co the bat atom
	double r_to_site[nNuclSite];	//khoang cach den cac site co the bat atom
} atom;

atom   free_atom[2*N_mean];	//mang luu tru cac atom tu do
double nucl_site[nNuclSite][4]; 	//luu tru 2 toa do, ban kinh khoi bac tai Nucleation site va so atom bi bat
double R=sqrt(D*dt);		//khoang cach bat cua tam
double sigma=R/sqrt(2);		//P(x)=1/(sigma*sqrt(2*pi)).exp((-x^2)/(D.dt)) => 2*sigma^2=D*dt => sigma=sqrt(D*dt/2)
int N=0;			//so hat tu do hien tai
gsl_rng *rng=gsl_rng_alloc(gsl_rng_mt19937);
unsigned long int Seed = time(NULL);
  
