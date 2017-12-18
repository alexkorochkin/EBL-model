#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <omp.h>
#include "gsl/gsl_math.h"
#include "gsl/gsl_const_mksa.h"
#include "spline.h"
#include "titles.h"
using namespace std;
//pr 241 21 500 1000 121 1000 10
struct Supp{
	static const int nspec = 241;
	static const int ndust = 21;
	static const int nr = 500;
	static const int nm = 1000;
	static const int bcnr = 1000;
	static const int nt = 1000;
	static const int nt_dust = 50;
	static const int npah = 8;
	double 	z_sfr_begin, z_ini_point, z_end_point, star_m_min, star_m_max,
			imf_x0, imf_D, imf_a, imf_norm, bc_size, bc_lifetime, bc_dust_density,
			bc_tau, bc_tau_slope, bc_efficiency, dust_size_min,
			dust_size_max, dust_dist_slope, sfr_scaling, N_C, N_H, pah_density,
			calzetti_A_v, gra_part;
	Supp(const vector<double> &);					
	vector<double> l;
	tk::spline spline_sfr;
	double a[ndust], gra[ndust][nspec], sil[ndust][nspec];
	double pah_cw[npah], pah_dc[npah], pah_ics[npah];};
Supp::Supp(const vector<double> &x){
	z_sfr_begin = x[0];
	z_ini_point = x[1];
	z_end_point = x[2];
	star_m_min = x[3];
	star_m_max = x[4];
	imf_x0 = x[5];
	imf_D = x[6];
	imf_a = x[7];
	imf_norm = x[8];
	bc_size = x[9];
	bc_lifetime = x[10];
	bc_dust_density = x[11];
	bc_tau = x[12];
	bc_tau_slope = x[13];
	bc_efficiency = x[14];
	dust_size_min = x[15];
	dust_size_max = x[16];
	dust_dist_slope = x[17];
	sfr_scaling = x[18];
	calzetti_A_v = x[19];
	pah_density = x[20];
	gra_part = 0.6;

	ifstream fin1, fin2;
	double w;
	char ch;
	string str;
	fin1.open("in/Gra_21.txt");
	fin2.open("in/Sil_21.txt");
	for(int h = 0; h < 5; h++) getline(fin1, str);
	for(int h = 0; h < 5; h++) getline(fin2, str); 	 
	for(int i = 0; i < ndust; i++){
		l.clear();
		fin1 >> a[i];
		fin2 >> w;
		a[i] = a[i]*1e-6;  //meter	
		for(int h = 0; h < 76; h++) fin1 >> ch;
		for(int h = 0; h < 61; h++) fin2 >> ch;	
		for(int j = nspec - 1; j > -1; j--){						
 			fin1 >> w;
			l.insert(l.begin(), w*1e-6);
			fin1 >> gra[i][j];
			fin1 >> w;
			fin1 >> w;
			fin2 >> w;
			fin2 >> sil[i][j];
			fin2 >> w;
			fin2 >> w;
		}}
	fin1.close();
	fin2.close();

	N_C = 60;
	N_H = 20;
	pah_cw[0] = 3.3e-6;
	pah_cw[1] = 6.2e-6;
	pah_cw[2] = 7.7e-6;
	pah_cw[3] = 8.6e-6;
	pah_cw[4] = 11.3e-6;
	pah_cw[5] = 11.9e-6;
	pah_cw[6] = 12.8e-6;
	pah_cw[7] = 13.6e-6;
	pah_dc[0] = 20e12;
	pah_dc[1] = 14e12;
	pah_dc[2] = 22e12;
	pah_dc[3] = 6e12;
	pah_dc[4] = 4e12;
	pah_dc[5] = 7e12;
	pah_dc[6] = 3.5e12;
	pah_dc[7] = 4e12;
	pah_ics[0] = 12*N_H;
	pah_ics[1] = 14*N_C;
	pah_ics[2] = 51*N_C;
	pah_ics[3] = 27*N_H;
	pah_ics[4] = 41*N_H;
	pah_ics[5] = 27*N_H;
	pah_ics[6] = 47*N_H;
	pah_ics[7] = 12*N_H;

	//imf
	double m, nm, dm, y1, y2, k;
	nm = 10000;
	w = 0;
	k = exp(-pow((log10(imf_x0)), 2)/(2*pow(imf_D, 2)));
	dm = (star_m_max - star_m_min)/(nm - 1);
	for(int i = 0; i < nm - 1; i++){
		m = star_m_min + i*dm;
		if (m < 1) 	y1 = (1/(m*imf_norm))*exp(-pow((log10(m) - log10(imf_x0)), 2)/
					(2*pow(imf_D, 2)));
		else 		y1 = (k/imf_norm)*pow(m, -imf_a);
		if (m + dm < 1) y2 = (1/((m + dm)*imf_norm))*exp(-pow((log10(m + dm) - log10(imf_x0)), 2)/
					(2*pow(imf_D, 2)));
		else 		y2 = (k/imf_norm)*pow(m + dm, -imf_a);
		w += (m*y1 + (m + dm)*y2)/2*dm;}
	imf_norm = imf_norm*w;
	//sfr
	vector<double> sfr1, sfr2; 
	fin1.open("in/sfr.txt");
	fin1 >> w;
	while(!fin1.eof()){
		sfr1.insert(sfr1.begin(), tz(w));
		fin1 >> w;
		sfr2.insert(sfr2.begin(), sfr_scaling*pow(10., w));
		fin1 >> w;}
	fin1.close();
	spline_sfr.set_points(sfr1, sfr2);}

double zt(const double t){
	double temp;
	double omega_lambda = 0.7;
	double omega_matter = 0.3;
	double H_rev_y = 1.46e10;
	temp = 	pow(omega_lambda/omega_matter, 1./3.)*
		   	pow(sinh(1.5*pow(omega_lambda, 0.5)*t/H_rev_y), -2./3.) - 1;
	return temp;}
double tz(const double z){
	double temp;
	double omega_lambda = 0.7;
	double omega_matter = 0.3;
	double H_rev_y = 1.46e10;
	temp = 	2./(3.*pow(omega_lambda, 0.5))*H_rev_y*
		   	log(pow(omega_lambda/(omega_matter*pow(1 + z, 3)), 0.5) + 
		   	pow(omega_lambda/(omega_matter*pow(1 + z, 3)) + 1, 0.5));
	return temp;} 
double rz(const double z){
	double nz, dz, temp;
	double omega_lambda = 0.69;
	double omega_matter = 0.31;
	double H_rev_y = 1.46e10;
	double secinyear = 31536000;
	dz = 0.0001;
	nz = z/dz;
	temp = 0;
	for(int i = 0; i < nz; i++)	
		temp +=	GSL_CONST_MKSA_SPEED_OF_LIGHT*secinyear*H_rev_y/
				(GSL_CONST_MKSA_PARSEC*1e6) * dz/pow(omega_matter* 
				pow(1 + (i + 0.5)*dz, 3) + omega_lambda, 0.5);
	return temp;}
double zr(const double r){
	double dz, z, temp;
	double omega_lambda = 0.76;
	double omega_matter = 0.24;
	double H_rev_y = 1.46e10;
	double secinyear = 31536000;
	dz = 0.0001;
	z = dz/2;
	temp = 0;
	while (r > temp) {
		temp +=	GSL_CONST_MKSA_SPEED_OF_LIGHT*secinyear*H_rev_y/
				(GSL_CONST_MKSA_PARSEC*1e6)* 
				dz/pow(omega_matter*pow(1 + z, 3) + omega_lambda, 0.5);
		z += dz;} 
	return z;}
double imf(const double m, const Supp *s){
	double y, k;
	k = exp(-pow((log10(s->imf_x0)), 2)/(2*pow(s->imf_D, 2)));
	if (m < 1) 	y = (1/(s->imf_norm*m))*exp(-pow((log10(m) - log10(s->imf_x0)), 2)/
					(2*pow(s->imf_D, 2)));
	else 		y = (k/s->imf_norm)*pow(m, -s->imf_a);
	return y;}
double sfr(const double t, const Supp *s){return s->spline_sfr(t);}
double calzetti(double lambda, const Supp *s){
	double R_v, A_lambda, k, k1, k2;
	lambda = lambda*1e6;	
	R_v = 4.05;
	k = 0;
	if(lambda < 0.12){
		k1 = 2.659*(-2.156 + 1.509/0.12 - 0.198/pow(0.12, 2) + 0.011/pow(0.12, 3)) + R_v;
		k2 = 2.659*(-2.156 + 1.509/0.11 - 0.198/pow(0.11, 2) + 0.011/pow(0.11, 3)) + R_v;
		k = k1 + (0.12 - lambda)*(k2 - k1)/(0.12 - 0.11);
	}
	if((0.12 <= lambda) && (lambda <= 0.63))
		k = 2.659*(-2.156 + 1.509/lambda - 0.198/pow(lambda, 2) + 0.011/pow(lambda, 3)) + R_v;
	if((0.63 < lambda) && (lambda <= 2.20))
		k = 2.659*(-1.857 + 1.040/lambda) + R_v; 
	A_lambda = k*s->calzetti_A_v/R_v;
	return pow(10, -0.4*A_lambda);}
void pah(const double *ss, double *ss1, const Supp *s){
	double flux, lambda_cut;
	flux = 0;
	lambda_cut = 1630e-10 + 370e-10*sqrt(s->N_C);
	for(int i = 0; i < s->nspec - 1; i++)
		if(s->l[i] < lambda_cut)
			flux += ss[i]*(s->l[i + 1] - s->l[i]);
	for(int i = 0; i < s->nspec; i++) ss1[i] = 0;
	for(int k = 0; k < s->npah; k++)
		for(int i = 0; i < s->nspec; i++) 
			ss1[i] += flux*s->pah_ics[k]*1e-26*
					(GSL_CONST_MKSA_SPEED_OF_LIGHT/pow(s->l[i], 2))*(s->pah_dc[k]*pow(GSL_CONST_MKSA_SPEED_OF_LIGHT/s->l[i], 2))/
					(pow(M_PI*(pow(GSL_CONST_MKSA_SPEED_OF_LIGHT/s->l[i], 2) - pow(GSL_CONST_MKSA_SPEED_OF_LIGHT/s->pah_cw[k], 2)), 2) 
					+ pow(s->pah_dc[k]*GSL_CONST_MKSA_SPEED_OF_LIGHT/(s->l[i]*2), 2));
	}

void print(const double *fs, const Supp *s){
	ofstream fout("out/spec.txt");
	for(int i = 0; i < s->nspec; i++) 
		fout << s->l[i]*1e10 << "   " << fs[i]*1e9*s->l[i] << endl;
	fout.close();}
void starspec(const double r, const double T, const double rad, double *ss, const Supp *s){
	for(int i = 0; i < s->nspec; i++) ss[i] =
		pow(rad/r, 2)*2*M_PI*GSL_CONST_MKSA_PLANCKS_CONSTANT_H*
		pow(GSL_CONST_MKSA_SPEED_OF_LIGHT, 2)*pow(s->l[i], -5)/
		(exp(GSL_CONST_MKSA_PLANCKS_CONSTANT_H*GSL_CONST_MKSA_SPEED_OF_LIGHT/
			(GSL_CONST_MKSA_BOLTZMANN*T*s->l[i])) - 1);}
void basicspec(const double t, const double r, double *bs, const Supp *s){
	double m, dm, met, Temp, Rad, am, tau, dd, vol, gtd, m_H, sunrad, 
			graph_density, hydr_mass, dust_density, norm;
	double ss[s->nspec], ss1[s->nspec], ss2[s->nspec];
	int imin, imax;
	Coef track;
	sunrad = 2.26e-14;
	met = 0.02;
	dm = (s->star_m_max - s->star_m_min)/(s->nm - 1);
	imin = 0;
	imax = 0;
	while (s->a[imin] < s->dust_size_min) imin++;
	while (s->a[imax] < s->dust_size_max) imax++;

	for(int i = 0; i < s->nm - 1; i++){
		m = s->star_m_min + (i + 0.5)*dm;
		if(t < track.lt(m, met)){
			Temp = track.Tt(t/1e6, m, met);
			Rad = pow(10, track.Radt(t/1e6, m, met))*sunrad;
			starspec(r, Temp, Rad, ss1, s);
			for(int j = 0; j < s->nspec; j++)
				bs[j] += imf(m, s)*ss1[j]*dm; 
		}
	}

	// double m, dm, dt, met, Temp1, Temp2, Rad1, Rad2, sunrad, tmax, tmin; 
	// double time[350];
	// int i_MS, i_BGB, i_HeI, i_He, i_DU, i_AGB, imin, imax;
	// i_MS = 350;
	// i_BGB = 350;
	// i_HeI = 350;
	// i_He = 350;
	// i_DU = 350;
	// i_AGB = 300;
	// imin = 0;
	// imax = 0;
	// Coef track;
	// sunrad = 2.26e-14;
	// met = 0.02;
	// dt = (tz(0) - tz(s->z_sfr_begin) - s->bc_lifetime)/(s->nt - 1);
	// dm = (s->star_m_max - s->star_m_min)/(s->nm - 1);
	// tmin = t*1e-6;
	// dt = dt*1e-6;
	// tmax = tmin + dt;
	
	// double **tmp;
	// tmp = new double*[s->nm];
	// for(int i = 0; i < s->nm; i++) tmp[i] = new double[s->nspec];
	// for(int i = 0; i < s->nm; i++) 
	// 	for(int j = 0; j < s->nspec; j++) tmp[i][j] = 0;		
	// double *cf;
	// cf = new double[350];
	// for(int i = 0; i < 350; i++) cf[i] = 1;
	// double *ss1;
	// ss1 = new double[s->nspec];
	// for(int i = 0; i < s->nspec; i++) ss1[i] = 0;
	// double *ss2;
	// ss2 = new double[s->nspec];
	// for(int i = 0; i < s->nspec; i++) ss2[i] = 0;
	// double *ss3;

	// for(int j = 0; j < s->nm - 1; j++){	 	
	// 	m = s->star_m_min + j*dm;
	// 	if(t*1e-6 < track.t_MS(m, met)){
	// 		if(m >= track.M_FGB(met)) i_BGB = 350;
	// 		if(2.25 < track.M_cBAGB(m, met)) i_DU = 350;
	// 		for(int i = 0; i < i_MS; i++) 			
	// 			time[i] = track.t_MS(m, met)*i/(i_MS - 1.);
	// 		for(int i = i_MS; i < i_BGB; i++)
	// 			time[i] = track.t_MS(m, met) + (track.t_BGB(m, met) - track.t_MS(m, met))*(i + 1. - i_MS)/(i_BGB - i_MS);
	// 		for(int i = i_BGB; i < i_HeI; i++) 
	// 			time[i] = track.t_BGB(m, met) + (track.t_HeI(m, met) - track.t_BGB(m, met))*(i + 1. - i_BGB)/(i_HeI - i_BGB);			
	// 		for(int i = i_HeI; i < i_He; i++) 
	// 			time[i] = track.t_HeI(m, met) + track.t_He(m, met)*(i + 1. - i_HeI)/(i_He - i_HeI);	
	// 		for(int i = i_He; i < i_DU; i++)
	// 			time[i] = track.t_HeI(m, met) + track.t_He(m, met) + (track.t_DU(m, met) - track.t_HeI(m, met) - track.t_He(m, met))*
	// 					  (i + 1. - i_He)/(i_DU - i_He);				
	// 		for(int i = i_DU; i < i_AGB; i++) 
	// 			time[i] = track.t_DU(m, met) + (track.t_AGB(m, met) - track.t_DU(m, met))*(i + 1. - i_DU)/(i_AGB - i_DU);
		
	// 		imin = 0;
	// 		imax = 0;
	// 		while((time[imin] <= tmin) && (imin < 349)) imin++;
	// 		imin--;
	// 		while((time[imax] < tmax) && (imax < 349)) imax++;
	// 		for(int i = 0; i < 350; i++) cf[i] = 1;
	// 		if((imax - imin) == 1) cf[imin] = (tmax - tmin)/(time[imax] - time[imin]);
	// 		else {
	// 			cf[imin] = (time[imin + 1] - tmin)/(time[imin + 1] - time[imin]);
	// 			cf[imax - 1] = (tmax - time[imax - 1])/(time[imax] - time[imax - 1]);
	// 		}
	// 		for(int i = 0; i < 35; i++){
	// 			cout << i << "   " << cf[i] << endl;
	// 			cin.get();
	// 		}			

	// 		Temp1 = track.Tt(tmin, m, met);
	// 		Rad1 = pow(10, track.Radt(tmin, m, met))*sunrad;
	// 		starspec(r, Temp1, Rad1, ss1, s);
	// 		for(int i = imin; i < imax; i++){
	// 			Temp2 = track.Tt(time[i + 1], m, met);
	// 			Rad2 = pow(10, track.Radt(time[i + 1], m, met))*sunrad;
	// 			starspec(r, Temp2, Rad2, ss2, s);
	// 			for(int k = 0; k < s->nspec; k++)
	// 				tmp[j][k] += (ss1[k] + ss2[k])/2*cf[i]*(time[i + 1] - time[i])/dt;
	// 			// for(int k = 0; k < s->nspec; k++)
	// 			// 	cout << m << "    " << s->l[k] << "    " << ss1[k] << "   " << ss2[k] << "    " << i << "    " << cf[i]*(time[i + 1] - time[i])/dt << endl;
	// 			// cin.get();
	// 			ss3 = ss1;
	// 			ss1 = ss2;
	// 			ss2 = ss3;
	// 		}	
	// 		// for(int k = 0; k < s->nspec; k++)
	// 		// 	cout << m << "    " << s->l[k] << "    " << tmp[j][k] << endl;
	// 		// cin.get();
	// 	}
	// }
	// for(int i = 0; i < s->nm - 2; i++){	
	// 	m = s->star_m_min + i*dm;
	// 	for(int j = 0; j < s->nspec; j++){
	// 		bs[j] += (imf(m, s)*tmp[i][j] + imf(m + dm, s)*tmp[i + 1][j])/2*dm;
	// 	}
	// }
	// // for(int j = 0; j < s->nspec; j++) cout << 1e10*s->l[j] << "   " << bs[j]*1e9*s->l[j] << endl;

	// for(int i = 0; i < s->nm; i++) delete tmp[i];
	// delete tmp;
	// delete ss1;
	// delete ss2;
	// delete cf;
	}
void basicspecdust(const double t, const double r, double *bs, const Supp *s){
	double Temp1, Temp2, Rad1, Rad2, tau, dt, norm, part;
	double ss[s->nspec], ss1[s->nspec], ss2[s->nspec];
	int imin, imax;

	double **ss3;
	ss3 = new double*[s->bcnr];
	for(int i = 0; i < s->bcnr; i++) ss3[i] = new double[s->nspec];
	
	double **q_abs;
	q_abs = new double*[s->ndust];
	for(int i = 0; i < s->ndust; i++) q_abs[i] = new double[s->nspec];
	
	for(int i = 0; i < s->bcnr; i++)
		for(int j = 0; j < s->nspec; j++) ss3[i][j] = 0;
	
	int bcnr_min = 3;
	part = s->gra_part;
	dt = s->bc_lifetime/(s->nt_dust - 1);
	imin = 0;
	imax = 0;
	while (s->a[imin] < s->dust_size_min) imin++;
	while (s->a[imax] < s->dust_size_max) imax++;
	imax--;
	for(int j = 0; j < s->nspec; j++){	
		ss[j] = 0;
		ss1[j] = 0;
		ss2[j] = 0;
	}
	basicspec(0, r, ss1, s);		
	for(int i = 1; i < int(t/dt) + 1; i++){
		basicspec(i*dt, r, ss2, s);
		for(int j = 0; j < s->nspec; j++){ 
			bs[j] += s->bc_efficiency*(dt/s->bc_lifetime)*(ss1[j] + ss2[j])/2;
		 	ss1[j] = ss2[j];
		 	ss2[j] = 0;
		 }
	}
	if(int(t/dt) == 0) 	for(int j = 0; j < s->nspec; j++) 
			bs[j] += s->bc_efficiency*(t/s->bc_lifetime)*ss1[j];



	norm = s->bc_dust_density*(s->dust_dist_slope - 1)/
		   (pow(s->dust_size_min, -s->dust_dist_slope + 1) - pow(s->dust_size_max, -s->dust_dist_slope + 1));


	for(int type = 0; type < 2; type++){
		if(type == 0) for(int i = 0; i < s->ndust; i++)
						for(int j = 0; j < s->nspec; j++) 
							q_abs[i][j] = s->gra[i][j]; 
		else		  {for(int i = 0; i < s->ndust; i++)
						for(int j = 0; j < s->nspec; j++) 
							q_abs[i][j] = s->sil[i][j];
					   part = 1 - part;}
		for(int i = bcnr_min; i < s->bcnr; i++){
			tau = s->bc_tau*(i - bcnr_min)/(s->bcnr - 1 - bcnr_min);
			for(int j = 0; j < s->nspec; j++)
				ss[j] = pow((r*(s->bcnr - 1.))/(s->bc_size*i), 2)*bs[j]*exp(-tau*pow(s->l[j]*1e10/5500, -s->bc_tau_slope));
			Temp1 = dusttemp(type, imin - 1, ss, s);
			Temp2 = dusttemp(type, imin, ss, s);
			Rad1 = s->a[imin - 1];
			Rad2 = s->a[imin];
			starspec(Rad1, Temp1, Rad1, ss1, s);
			starspec(Rad2, Temp2, Rad2, ss2, s);
			for(int j = 0; j < s->nspec; j++)
				ss3[i][j] += part*norm*(q_abs[imin - 1][j]*pow(Rad1, -s->dust_dist_slope)*pow(Rad1/r, 2)*ss1[j] +
					         q_abs[imin][j]*pow(Rad2, -s->dust_dist_slope)*pow(Rad2/r, 2)*ss2[j])/2*(Rad2 - s->dust_size_min);
			for(int q = imin; q < imax; q++){
				Temp1 = Temp2;
				Temp2 = dusttemp(type, q + 1, ss, s);
				Rad1 = s->a[q];
				Rad2 = s->a[q + 1];
				starspec(Rad1, Temp1, Rad1, ss1, s);
				starspec(Rad2, Temp2, Rad2, ss2, s);
				for(int j = 0; j < s->nspec; j++)
					ss3[i][j] += part*norm*(q_abs[q][j]*pow(Rad1, -s->dust_dist_slope)*pow(Rad1/r, 2)*ss1[j] +
							   	 q_abs[q + 1][j]*pow(Rad2, -s->dust_dist_slope)*pow(Rad2/r, 2)*ss2[j])/2*(Rad2 - Rad1);
			}
			Temp1 = Temp2;
			Temp2 = dusttemp(type, imax + 1, ss, s);
			Rad1 = s->a[imax];
			Rad2 = s->a[imax + 1];
			starspec(Rad1, Temp1, Rad1, ss1, s);
			starspec(Rad2, Temp2, Rad2, ss2, s);
			for(int j = 0; j < s->nspec; j++)
				ss3[i][j] += part*norm*(q_abs[imax][j]*pow(Rad1, -s->dust_dist_slope)*pow(Rad1/r, 2)*ss1[j] +
						  	 q_abs[imax + 1][j]*pow(Rad2, -s->dust_dist_slope)*pow(Rad2/r, 2)*ss2[j])/2*(s->dust_size_max - Rad1);
		}
	}
	for(int i = bcnr_min; i < s->bcnr; i++){
		tau = s->bc_tau*(i - bcnr_min)/(s->bcnr - 1 - bcnr_min);
		for(int j = 0; j < s->nspec; j++)
			ss[j] = pow((r*(s->bcnr - 1))/(s->bc_size*i), 2)*bs[j]*exp(-tau*pow(s->l[j]*1e10/5500, -s->bc_tau_slope));
		pah(ss, ss1, s);
		for(int j = 0; j < s->nspec; j++)
			ss3[i][j] += s->pah_density*ss1[j];
	}

	for(int i = bcnr_min; i < s->bcnr - 1; i++)
		for(int j = 0; j < s->nspec; j++)
			bs[j] += (4./3.)*M_PI*GSL_CONST_MKSA_PARSEC*1e6*
					 (pow(s->bc_size*(i + 1)/(s->bcnr - 1), 3) - pow(s->bc_size*i/(s->bcnr - 1), 3))*(ss3[i + 1][j] + ss3[i][j])/2;
	for(int i = 0; i < s->nspec; i++)
		bs[i] = bs[i]*exp(-s->bc_tau*pow(s->l[i]*1e10/5500, -s->bc_tau_slope));
	for (int i = 0; i < s->bcnr; i++) delete ss3[i];
	delete ss3;
	for (int i = 0; i < s->ndust; i++) delete q_abs[i];
	delete q_abs;}
double dusttemp(const int type, const int q, const double *fs, const Supp *s){
	double Tmin, Tmax, T, delta, N1, N2, area;			 
	double q_abs[s->nspec];
	Tmin = 0;
	Tmax = 4000;
	T = (Tmin + Tmax)/2;
	delta = 0.05;
	N1 = 0;
	N2 = 0;
	area = 4;
	if(type == 0) 	for(int i = 0; i < s->nspec; i++) q_abs[i] = s->gra[q][i];
	else 			for(int i = 0; i < s->nspec; i++) q_abs[i] = s->sil[q][i];

	for(int i = 0; i < s->nspec - 1; i++)
		N1 += (fs[i]*q_abs[i] + fs[i + 1]*q_abs[i + 1])/2*(s->l[i + 1] - s->l[i]);
	while (Tmax - Tmin > delta){
		for(int i = 0; i < s->nspec - 1; i++)			
			N2 += area*2*M_PI*GSL_CONST_MKSA_PLANCKS_CONSTANT_H*pow(GSL_CONST_MKSA_SPEED_OF_LIGHT, 2)* 
				(pow(s->l[i], -5)/(exp(GSL_CONST_MKSA_PLANCKS_CONSTANT_H*GSL_CONST_MKSA_SPEED_OF_LIGHT/
									(GSL_CONST_MKSA_BOLTZMANN*T*s->l[i])) - 1)*q_abs[i] + 
				 pow(s->l[i + 1], -5)/(exp(GSL_CONST_MKSA_PLANCKS_CONSTANT_H*GSL_CONST_MKSA_SPEED_OF_LIGHT/
									(GSL_CONST_MKSA_BOLTZMANN*T*s->l[i + 1])) - 1)*q_abs[i + 1])/2*(s->l[i + 1] - s->l[i]);
		if(N2 > N1) Tmax = T;
		else 		Tmin = T;
		T = (Tmin + Tmax)/2;
		N2 = 0;
	}
	return T;}
void spec(const vector<double> &x){
	Supp* s = new Supp(x);
	double dr = rz(s->z_sfr_begin)/s->nr;
	double dt = (tz(0) - tz(s->z_sfr_begin) - s->bc_lifetime)/(s->nt - 1);
	double dt_dust = s->bc_lifetime/(s->nt_dust - 1);
	double fs[s->nspec];
	int imin = rz(s->z_end_point)/dr + 1;
	int imax = rz(s->z_sfr_begin)/dr;
	double **tmp;
	tmp = new double*[s->nt];
	for (int i = 0; i < s->nt; i++) tmp[i] = new double[s->nspec];
	double **tmp_dust;
	tmp_dust = new double*[s->nt_dust];
	for (int i = 0; i < s->nt_dust; i++) tmp_dust[i] = new double[s->nspec];
	for (int j = 0; j < s->nspec; j++) fs[j] = 0;


 //********* Cloud Test *********//

	double ro, m_H, mass, norm, full_mass, gas_density;
	ro = 2.23*1e6;   //   g/m^3
	m_H = 1.67*1e-24;
	norm = s->bc_dust_density*(s->dust_dist_slope - 1)/
		   (pow(s->dust_size_min, -s->dust_dist_slope + 1) - pow(s->dust_size_max, -s->dust_dist_slope + 1));
	mass = fabs(ro*(4./3.)*M_PI*norm*(pow(s->dust_size_max, -s->dust_dist_slope + 3 + 1) - pow(s->dust_size_min, -s->dust_dist_slope + 3 + 1))
			/(-s->dust_dist_slope + 3 + 1));
	gas_density = mass*100/(m_H*1e6);
	full_mass = 4/3*M_PI*pow(s->bc_size*GSL_CONST_MKSA_PARSEC*1e6, 3)*mass/(2*1e33)*100;
	// cout << "cloud mass = " << full_mass << endl;
	// cout << "star mass = " << s->bc_efficiency << endl;
	// cout << "sfe = " << s->bc_efficiency/full_mass << endl;
	// cout << "gas density = " << gas_density << endl;
	// cin.get();
 /////////////////////////////////

 //********* Calzetti extinction law ***********//
	// for(int i = 0; i < s->nspec; i++){
	// 	cout << s->l[i]*1e10 << "   " << log(calzetti(s->l[i], s)) << "   " 
	// 		 << log(exp(-s->bc_tau*pow(s->l[i]*1e10/5500, -s->bc_tau_slope))) << endl;
	// 	cin.get();
	// }
 /////////////////////////////////////////////////

	#pragma omp parallel
	{	
		#pragma omp for schedule(dynamic)
		for (int i = 0; i < s->nt; i++){
			for(int j = 0; j < s->nspec; j++) tmp[i][j] = 0.0;
			basicspec(s->bc_lifetime + i*dt, dr, tmp[i], s);
			for(int j = 0; j < s->nspec; j++)
				tmp[i][j] = calzetti(s->l[j], s)*tmp[i][j];}
			// cout << i << endl;
		#pragma omp for schedule(dynamic)
		for (int i = 0; i < s->nt_dust; i++){
			for(int j = 0; j < s->nspec; j++) tmp_dust[i][j] = 0.0;
			basicspecdust((i + 0.5)*dt_dust, dr, tmp_dust[i], s);
		// cout << i << endl;
		}
		#pragma omp for schedule(dynamic)
		for (int i = imin; i < imax; i++){
			vector<double> bs;
			bs.clear();
			tk::spline splspec;
			double z, localtime, starage;
			z = zr(i*dr);
			for (int j = 0; j < s->nspec; j++) bs.push_back(0);
			for (int j = 0; j < s->nt_dust - 1; j++){
				localtime = tz(z) - j*dt_dust;
				if(localtime > tz(s->z_sfr_begin)) for(int k = 0; k < s->nspec; k++){
					bs[k] += (sfr(localtime, s)*tmp_dust[j][k] + sfr(localtime - dt_dust, s)*tmp_dust[j + 1][k])/2*
							 (1./s->bc_efficiency)*dt_dust*pow(dr, 3);
				}
			}
			for (int j = 0; j < s->nt - 1; j++){
				localtime = tz(z) - j*dt - s->bc_lifetime;	
				if(localtime > tz(s->z_sfr_begin)) for(int k = 0; k < s->nspec; k++) 
					bs[k] += (sfr(localtime, s)*tmp[j][k] + sfr(localtime - dt, s)*tmp[j + 1][k])/2*
							  dt*pow(dr, 3);
			}
			splspec.set_points(s->l, bs);
			#pragma omp critical
			{ 
				for(int j = 0; j < s->nspec; j++) fs[j] +=
					splspec(s->l[j]/(1 + z))/pow(1 + z, 3);
			}
		}
	}

	print(fs, s);
	delete s;
	for (int i = 0; i < s->nt; i++) delete tmp[i];
	delete tmp;
	for (int i = 0; i < s->nt_dust; i++) delete tmp_dust[i];
	delete tmp_dust;}
  
// int main(){
// 	vector<double> x, dx, xmin, xmax, propx;
// 	ifstream fin("in/param.txt");
// 	double q;
// 	string a;
// 	while(fin >> a){
// 		fin >> q;
// 		x.push_back(q);
// 		propx.push_back(q);
// 		fin >> q;
// 		dx.push_back(q);
// 		fin >> q;
// 		xmin.push_back(q);
// 		fin >> q;
// 		xmax.push_back(q);
// 	}
// 	fin.close();
// 	spec(x);}

