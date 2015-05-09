#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "define.h"

double R=sqrt(D*dt);

void site_position_setup() {
	int i;
	for (i=0;i<nNuclSite;i++) {
		nucl_site[i][0]=L*gsl_rng_uniform(rng);
		nucl_site[i][1]=L*gsl_rng_uniform(rng);
		nucl_site[i][2]=0; //ban kinh hat Ag tai site ban dau =0
		nucl_site[i][3]=0; //so atom Ag bi bat tai site ban dau =0
	}
}	

void add_particle() {
   int i;
   for (i=0;i<nParticlePerStep;i++) {
	N++;
	free_atom[N].pos[0]=L*gsl_rng_uniform(rng);
	free_atom[N].pos[1]=L*gsl_rng_uniform(rng);
   }
}

void Brown() {
   int idx;
   for (idx=0;idx<N;idx++) {
	double x=free_atom[idx].pos[0];
	x+=gsl_ran_gaussian(rng,sigma);
	//dieu kien bien tuan hoan
	if (x>L) x-= ((int) (x/L)) * L;
	if (x<0) x-= ((int) (x/L) - 1) * L ;
	free_atom[idx].pos[0]=x;
	
	x=free_atom[idx].pos[1];
	x+=gsl_ran_gaussian(rng, sigma); 
	//dieu kien bien tuan hoan
	if (x>L) x-= ((int) (x/L)) * L;
	if (x<0) x-= ((int) (x/L) - 1) * L ;
	free_atom[idx].pos[1]=x;
   }
}

void atom_catching(int idx, int site_idx) {
	int i;
	//xoa atom idx khoi mang 
	for (i=idx;i<N-1;i++) free_atom[i]=free_atom[i+1];
	N--;
	//ban kinh cua site bat duoc atom tang len khi bat 1 atom Ag
	double n=nucl_site[site_idx][3]+nAtomPerParticle;	
	if (nucl_site[site_idx][3]==0) nucl_site[site_idx][2]=0.5*a*pow(3/M_PI*nAtomPerParticle,1.0/3);  
	 else nucl_site[site_idx][2]*=pow(n/(nucl_site[site_idx][3]),1.0/3);
	nucl_site[site_idx][3]=n;
}

void catching_list(int idx) {
	int i, n=0;
//	free_atom[idx].n_catch=0;
	for (i=0;i<nNuclSite; i++) {
		double dx=free_atom[idx].pos[0]-nucl_site[i][0];
		double dy=free_atom[idx].pos[1]-nucl_site[i][1];
		double dr=sqrt(dx*dx+dy*dy); //khoang cach tu atom idx den site i
		if (dr<=R+nucl_site[i][2]) {
			n++;
			free_atom[idx].catch_idx[n-1]=i;
			free_atom[idx].r_to_site[n-1]=dr;
		}
	}
	
	free_atom[idx].n_catch=n;	
}

void catching_process() {
	int i;
	for (i=0;i<N;i++) {
		int n=free_atom[i].n_catch;	
		catching_list(i);
		if (n>0) {
		 int catch_site_idx; // index cua site bat duoc atom i
		 if (n==1) catch_site_idx=free_atom[i].catch_idx[0]; 
		 else { 
			double *P, P_total; //xac suat bat atom
			P=(double*)calloc(n,sizeof(double));
			int j;
			for (j=0;j<n;j++) { 
				double r2=free_atom[i].r_to_site[j];
				r2=r2*r2;	
				P[j]=exp(-r2/(R*R));
				P_total+=P[j];
			}
			for (j=0;j<n;j++) P[j]/=P_total; //normalization
		 	double x=gsl_rng_uniform(rng);
			double p=0;
			for (j=0;j<n;j++) {
			 p+=P[j];
			 if (x<=p) {
				catch_site_idx=free_atom[i].catch_idx[j]; 
				break;				
			 }
			}
		 }
		 atom_catching(i,catch_site_idx);		
		}
	}		
}

int main() {
	site_position_setup();
	gsl_rng_set(rng,Seed);
	
	int i;
	for (i=0;i<n_steps;i++) {
		add_particle();
		catching_process();
		Brown();
	}	
	printf("Tong so hat Ag den tam Pt: %d\n", n_steps*nParticlePerStep );
	printf("So hat bi bat: %d\n", n_steps*nParticlePerStep-N);
	FILE *f;
	f=fopen("site_size.txt","w");
	for (i=0;i<nNuclSite;i++) fprintf(f, "%d\t%lf\t%d\n", i, nucl_site[i][2], (int)nucl_site[i][3]);
	fclose(f);
	return 0;
}


