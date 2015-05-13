/*
   ###################################################
   # Ag particles formation process on Pt electrode  #		
   # 	Phuc Doan - K57 Honor Program of Physics     #
   # 	Hanoi University of Science - VNU            #
   #    Hanoi May-13				     #
   ###################################################
*/

#ifndef __DEFINE_H
#define __DEFINE_H

#define S 			1e6	//(nm2) dien tich mo phong

#define nSiteMax 		1000	//so site Max tren mang L*L
#define A			150e9   //nSite = nSiteMax * ( 1- exp(-At) )

#define nAtomPerParticle 	247	//so atom Ag trong moi hat bac 2nm
#define nParticlePerStep 	100	//so hat bac 2nm sinh ra trong moi time step
#define n_steps   		1000	//so buoc thoi gian

#define t1 			5.33e5	//(ns.nm2) thoi gian de co 1 atom Ag tren dien tich 1 nm2

#define D 			10	//(nm2/ns) do linh dong cua Ag tren Pt
#define a 			0.409	//(nm) hang so mang tinh the cua Ag

#define dMax			20	//Maximum size of site

#define N_mean ( (int) ( 3.16e3 / sqrt(D*dt) * nParticlePerStep ) ) //Average value cua so atom tu do

// cac thong tin gan voi moi atom
typedef struct {
	double pos[2];		//2 toa do cua atom
	int n_catch;		//so luong cac site co the bat atom do
	double catch_idx[nSiteMax];	//index cua cac site co the bat atom
	double r_to_site[nSiteMax];	//khoang cach den cac site co the bat atom
} atom;

void site_position_setup();			//setup toa do cua cac nucleation site
void add_particle(); 				//them 1 Ag vao mang
void Brown();					//chuyen dong Brown cua atom tu do tren Pt
void atom_catching(int idx, int site_idx);	//atom idx bi bat boi tam site_idx
void catching_list(int idx);			//liet ke cac site co the bat atom idx
void catching_process();			//qua trinh bat atom tren ca tam Pt


#endif
