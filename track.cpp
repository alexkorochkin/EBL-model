#include <cmath>
#include <ctime>
// #include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;
struct Coef{
	Coef();
	double p[98][5];
	double q[58][5];
	double func_a(int n, double met);
	double func_b(int n, double met);
	double a(int n, double met);
	double b(int n, double met);

	double t_hook(double M, double met);
	double t_MS(double M, double met);
	double t_BGB(double M, double met);
	double t_HeI(double M, double met);
	double t_HeMS(double M);
	double t_He(double M, double met);
	double tau_bl(double M, double met);
	double t_DU(double M, double met);
	double t_AGB(double M, double met);

	double M_cBGB(double M, double met);
	double M_cHeI(double M, double met);
	double M_cBAGB(double M, double met);
	double M_cGB(double t, double M, double met);
	double M_c(double t, double M, double met);

	double M_HeF(double met);
	double M_FGB(double met);
	double L_BGB(double M, double met);
	
	double L_ZAMS(double M, double met);
	double R_ZAMS(double M, double met);
	
	double L_TMS(double M, double met);
	double R_TMS(double M, double met);
	
	double L_EHG(double M, double met);
	double R_EHG(double M, double met);
	
	double L_ZHe(double M);
	double R_ZHe(double M);

	double L_minHe(double M, double met);
	double R_mHe(double M, double met);

	double L_HeI(double M, double met);	
	double R_HeI(double M, double met);

	double L_ZAHB(double M, double met);
	double R_ZAHB(double M, double met);

	double L_BAGB(double M, double met);
	double R_BAGB(double M, double met);

	double R_GBadd(double L, double M, double met);
	double R_AGBadd(double L, double M, double met);

	double L_MS(double t, double M, double met);
	double R_MS(double t, double M, double met);
	double L_HG(double t, double M, double met);
	double R_HG(double t, double M, double met);
	double L_GB(double t, double M, double met);
	double R_GB(double t, double M, double met);
	double L_CHeB(double t, double M, double met);
	double R_CHeB(double t, double M, double met);
	double L_EAGB(double t, double M, double met);
	double R_EAGB(double t, double M, double met);
	double L_TPAGB(double t, double M, double met);
	double R_TPAGB(double t, double M, double met);

	double Lumt(double t, double M, double met);
	double Radt(double t, double M, double met);
	double Tt(double t, double M, double met);
	double lt(double M, double met);};
Coef::Coef(){
	double x;
	ifstream fin("in/a.txt");
	for(int i = 1; i < 98; i++){
		for(int j = 0; j < 5; j++){
			fin >> x;
			p[i][j] = x;
		}
	}
	fin.close();
	fin.open("in/b.txt");
	for(int i = 1; i < 58; i++){
		for(int j = 0; j < 5; j++){
			fin >> x;
			q[i][j] = x;
		}
	}
	fin.close();}
double Coef::func_a(int n, double met){
	double x, zeta;
	zeta = log10(met/0.02);
	x = p[n][0] + p[n][1]*zeta + p[n][2]*pow(zeta, 2) +
		p[n][3]*pow(zeta, 3) + p[n][4]*pow(zeta, 4);
	return x;}
double Coef::func_b(int n, double met){
	double x, zeta;
	zeta = log10(met/0.02);
	x = q[n][0] + q[n][1]*zeta + q[n][2]*pow(zeta, 2) +
		q[n][3]*pow(zeta, 3) + q[n][4]*pow(zeta, 4);
	return x;}
double Coef::a(int n, double met){
	double x, y, z, zeta, sigma;
	zeta = log10(met/0.02);
	sigma = log10(met);
	switch(n) {
	case 11:
		x = func_a(11, met);
		y = func_a(14, met);
		return x*y;
		break;
	case 12:
		x = func_a(12, met);
		y = func_a(14, met);
		return x*y;
		break;
	case 17:
		x = max(0.097 - 0.1072*(sigma + 3), max(0.0097, min(0.1461, 0.1461 +
			0.1237*(sigma + 2))));
		return pow(10, x);
		break;
	case 18:
		x = func_a(18, met);
		y = func_a(20, met);
		return x*y;
		break;
	case 19:
		x = func_a(19, met);
		y = func_a(20, met);
		return x*y;
		break;
	case 29:
		x = func_a(29, met);
		y = func_a(32, met);
		return pow(x, y);
	case 33:
		x = min(1.4, 1.5135 + 0.3769*zeta);
		x = max(0.6355 - 0.4192*zeta, max(1.25, x));
		return x;
		break; 
	case 42:
		x = func_a(42, met);
		x = min(1.25, max(1.1, x));
		return x;
		break;
	case 44:
		x = func_a(44, met);
		x = min(1.3, max(0.45, x));
		return x;
		break;
	case 49:
		x = func_a(49, met);
		x = max(x, 0.145);
		return x;
		break;	
	case 50:
		x = func_a(50, met);
		x = min(x, 0.306 + 0.053*zeta);
		return x;
		break;
	case 51:
		x = func_a(51, met);
		x = min(x, 0.3625 + 0.062*zeta);
		return x;
		break;
	case 52:
		x = func_a(52, met);
		if(met <= 0.01) {x = max(x, 0.9);}
		else {x = min(x, 1.0);}
		return x;
		break;
	case 53:
		x = func_a(53, met);
		if(met <= 0.01) {x = max(x, 1.0);}
		else {x = min(x, 1.1);}
		return x;
		break;
	case 57:
		x = func_a(57, met);
		x = min(1.4, x);
		x = max(0.6355 - 0.4192*zeta, max(1.25, x));
		return x;
		break;
	case 62:
		x = func_a(62, met);
		x = max(0.065, x);
		return x;
		break;
	case 63:
		x = func_a(63, met);
		if(met < 0.004) {x = min(0.055, x);}
		return x;
		break;
	case 64:
		x = func_a(64, met);
		y = func_a(66, met);
		y = max(y, min(1.6, -0.308 - 1.046*zeta));
		y = max(0.8, min(0.8 - 2.0*zeta, y));
		z = func_a(68, met);
		z = max(0.9, min(z, 1.0));
		if(z > y) {x = (a(58, met)*pow(y, a(60, met)))/
					   (a(59, met)*pow(y, a(61, met)));}
		else {x = max(0.091, min(0.121, x));}
		return x;
		break;
	case 66:
		x = func_a(66, met);
		x = max(x, min(1.6, -0.308 - 1.046*zeta));
		x = max(0.8, min(0.8 - 2.0*zeta, x));
		return x;
		break;
	case 68:
		y = func_a(66, met);
		y = max(y, min(1.6, -0.308 - 1.046*zeta));
		y = max(0.8, min(0.8 - 2.0*zeta, y));
		z = func_a(68, met);
		z = max(0.9, min(z, 1.0));
		x = min(y, z);
		return x;
		break;
	case 72:
		x = func_a(72, met);
		if(met > 0.01) {x = max(x, 0.95);}
		return x;
		break;
	case 74:
		x = func_a(74, met);
		x = max(1.4, min(x, 1.6));
		return x;
		break;
	case 75:
		x = func_a(75, met);
		x = max(1.0, min(x, 1.27));
		x = max(x, 0.6355 - 0.4192*zeta);
		return x;
		break;
	case 76:
		x = func_a(76, met);
		x = max(x, -0.1015564 - 0.2161264*zeta - 0.05182516*pow(zeta, 2));
		return x;
		break;
	case 77:
		x = func_a(77, met);
		x = max(-0.3868766 - 0.5457078*zeta - 0.1463472*pow(zeta, 2), min(0.0, x));
		return x;
		break;
	case 78:
		x = func_a(78, met);
		x = max(0.0, min(x, 7.454 + 9.046*zeta));
		return x;
		break;
	case 79:
		x = func_a(79, met);
		x = min(x, max(2.0, -13.3 - 18.6*zeta));
		return x;
		break;
	case 80:
		x = func_a(80, met);
		x = max(0.0585542, x);
		return x;
		break;
	case 81:
		x = func_a(81, met);
		x = min(1.5, max(0.4, x));
		return x;
		break;
	default:
		x = func_a(n, met);
		return x;
		break;
	}}
double Coef::b(int n, double met){
	double x, y, zeta, sigma, ro;
	zeta = log10(met/0.02);
	sigma = log10(met);
	ro = zeta + 1.0;
	switch(n) {
	case 1:
		x = func_b(1, met);
		x = min(0.54, x);
		return x;
		break;
	case 2:
		x = pow(10, -4.6739 - 0.9394*sigma);
		x = min(max(x, -0.04167 + 55.67*met), 0.4771 - 9329.21*pow(met, 2.94));
		return x;
		break;
	case 3:
		x = max(-0.1451, -2.2794 - 1.5175*sigma - 0.254*pow(sigma, 2));
		x = pow(10, x);
		if(met > 0.004) {x = max(x, 0.7307 + 14265.1*pow(met, 3.395));}
		return x;
		break;
	case 4:
		x = func_b(4, met);
		x = x + 0.1231572*pow(zeta, 5);
		return x;
		break;
	case 6:
		x = func_b(6, met);
		x = x + 0.01640687*pow(zeta, 5);
		return x;
		break;
	case 11:
		x = func_b(11, met);
		x = pow(x, 2);
		return x;
		break;
	case 13:
		x = func_b(13, met);
		x = pow(x, 2);
		return x;
		break;
	case 14:
		x = func_b(14, met);
		y = func_b(15, met);
		return pow(x, y);
		break;
	case 16:
		x = func_b(16, met);
		y = func_b(15, met);
		return pow(x, y);
		break;
	case 17:
		x = 1.0;
		if(zeta > -1.0) {x = 1.0 - 0.3880523*pow(zeta + 1.0, 2.862149);}
		return x;
		break;
	case 24:
		x = func_b(24, met);
		y = func_b(28, met);
		return pow(x, y);
		break;
	case 26:
		x = 5.0 - 0.09138012*pow(met, -0.3671407);
		return x;
		break;
	case 27:
		x = func_b(27, met);
		y = func_b(28, met);
		return pow(x, 2*y);
		break;
	case 31:
		x = func_b(31, met);
		y = func_b(33, met);
		return pow(x, y);
		break;
	case 34:
		x = func_b(34, met);
		y = func_b(33, met);
		return pow(x, y);
		break;
	case 36:
		x = func_b(36, met);
		return pow(x, 4);
		break;
	case 37:
		x = func_b(37, met);
		return 4.0*x;
		break;
	case 38:
		x = func_b(38, met);
		return pow(x, 4);
		break;
	case 40:
		x = func_b(40, met);
		x = max(x, 1.0);
		return x;
		break;
	case 41:
		x = func_b(41, met);
		y = func_b(42, met);
		return pow(x, y);
		break;
	case 44:
		x = func_b(44, met);
		return pow(x, 5);
		break;
	case 45:
		x = 1.0 - (2.47162*ro - 5.401682*pow(ro, 2) + 3.247361*pow(ro, 3));
		if(ro <= 0) {x = 1.0;}
		return x;
		break;
	case 46:
		double M_HeF, M_FGB;
		M_HeF = 1.995 + 0.25*zeta + 0.087*pow(zeta, 2);
		M_FGB = 13.048*pow(met/0.02, 0.06)/(1 + 0.0012*pow(0.02/met, 1.27));
		x = func_b(46, met);
		x = -1.0*x*log10(M_HeF/M_FGB);
	case 47:
		x = 1.127733*ro + 0.2344416*pow(ro, 2) - 0.3793726*pow(ro, 3);
		return x;
		break;
	case 51:
		x = func_b(51, met);
		x = x - 0.1343798*pow(zeta, 5);
		return x;
		break;
	case 53:
		x = func_b(53, met);
		x = x + 0.4426929*pow(zeta, 5);
		return x;
		break;
	case 55:
		x = func_b(55, met);
		x = min(0.99164 - 743.123*pow(met, 2.83), x);
		return x;
		break;
	case 56:
		x = func_b(56, met);
		x = x + 0.1140142*pow(zeta, 5);
		return x;
		break;
	case 57:
		x = func_b(57, met);
		x = x - 0.01308728*pow(zeta, 5);
		return x;
		break;
	default:
		x = func_b(n, met);
		return x;
		break;
	}}

double Coef::t_hook(double M, double met){
	double mu;
	mu = max(0.5, 1.0 - 0.01*max(a(6, met)/
			pow(M, a(7, met)), a(8, met) + a(9, met)/pow(M, a(10, met))));
	return mu*t_BGB(M, met);}
double Coef::t_MS(double M, double met){
	double x, zeta;
	zeta = log10(met/0.02);
	x = max(0.95, min(0.95 - 0.03*(zeta + 0.30103), 0.99));
	return max(t_hook(M, met), x*t_BGB(M, met));}
double Coef::t_BGB(double M, double met){
	double x, y;
	x = a(1, met) + a(2, met)*pow(M, 4) + a(3, met)*pow(M, 5.5) + pow(M, 7);
	y = a(4, met)*pow(M, 2) + a(5, met)*pow(M, 7);
	return x/y;}
double Coef::t_HeI(double M, double met){
	double L_x, x, y, tinf1, tinf2, tx, A_H, B, D, D0, p, q, zeta;

	A_H = max(-4.8, min(-5.7 + 0.8*M, -4.1 + 0.14*M));
	A_H = pow(10, A_H);

	if(M <= M_HeF(met)) {p = 6;}
	if(M >= 2.5) 		{p = 5;}
	if((M > M_HeF(met)) && (M < 2.5))
						{p = 6 - (M - M_HeF(met))/(2.5 - M_HeF(met));}

	if(M <= M_HeF(met)) {q = 3;}
	if(M >= 2.5) 		{q = 2;}
	if((M > M_HeF(met)) && (M < 2.5))
	 					{q = 3 - (M - M_HeF(met))/(2.5 - M_HeF(met));}
	
	zeta = log10(met/0.02);	
	B = max(3*1e4, 500 + 1.75*1e4*pow(M, 0.6));
	D0 = 5.37 + 0.135*zeta;
	if(M <= M_HeF(met))	{D = D0;}
	if(M >= 2.5)		{D = max(-1.0, max(0.975*D0 - 0.18*M, 0.5*D0 - 0.06*M));}
	if((M > M_HeF(met)) && (M < 2.5))
			{x = max(-1.0, max(0.975*D0 - 0.18*2.5, 0.5*D0 - 0.06*2.5)) - D0;
			 y = 2.5 - M_HeF(met);
			 D = D0 + (M - M_HeF(met))*x/y;}
	D = pow(10, D);

	L_x 	= B*pow(B/D, q/(p - q));
	tinf1 	= t_BGB(M, met) + pow(D/L_BGB(M, met), (p - 1)/p) / ((p - 1)*A_H*D);
	tx 		= tinf1 - (tinf1 - t_BGB(M, met))*pow(L_BGB(M, met)/L_x, (p - 1)/p);
	tinf2 	= tx + pow(B/L_x, (q - 1)/q) / ((q - 1)*A_H*B);

	if(L_HeI(M, met) <= L_x) {return tinf1 - 
			pow(D, (p - 1)/p)/((p - 1)*A_H*D*pow(L_HeI(M, met), (p - 1)/p));}
	if(L_HeI(M, met) > L_x) {return tinf2 - 
			pow(B, (q - 1)/q)/((q - 1)*A_H*B*pow(L_HeI(M, met), (q - 1)/q));}}
double Coef::t_HeMS(double M){
	double x, y;
	x = 0.4129 + 18.81*pow(M, 4) + 1.853*pow(M, 6);
	y = pow(M, 6.5);
	return x/y;}
double Coef::t_He(double M, double met){
	double alpha4, mu, t_He, x, y, M_c;
	M_c = M_cHeI(M, met);
	if(M >= M_HeF(met)) {
		x = b(41, met)*pow(M, b(42, met)) + b(43, met)*pow(M, 5);
		y = b(44, met) + pow(M, 5);
		t_He = t_BGB(M, met)*x/y;
		return t_He;
	}
	if(M < M_HeF(met)) {
		alpha4 = (Coef::t_He(M_HeF(met), met) - b(39, met))/b(39, met);
		mu = (M - M_c)/(M_HeF(met) - M_c);
		t_He = (b(39, met) + (t_HeMS(M_c) - b(39, met))*pow(1 - mu, b(40 ,met)))*
				(1 + alpha4*exp(15*(M - M_HeF(met))));
		return t_He;
	}}
double Coef::tau_bl(double M, double met){
	double alpha_bl, f_blM, f_blM_FGB, tau_bl;
	alpha_bl = 	(1 - b(45, met)*pow(M_HeF(met)/M_FGB(met), 0.414));

	if((1 - R_mHe(M, met)/R_AGBadd(L_HeI(M, met), M, met)) > 0){
		f_blM = pow(M, b(48, met))*
				pow(1 - R_mHe(M, met)/R_AGBadd(L_HeI(M, met), M, met), b(49, met));}
	else {f_blM = 0;}
	f_blM_FGB = pow(M_FGB(met), b(48, met))*
				pow(1 - R_mHe(M_FGB(met), met)/R_AGBadd(L_HeI(M_FGB(met), met), M, met), b(49, met));
	
	if(M < M_HeF(met)) {
		tau_bl = 1;
		return tau_bl;}
	
	if((M >= M_HeF(met)) && (M <= M_FGB(met)))  {
		tau_bl = b(45, met)*pow(M/M_FGB(met), 0.414) + alpha_bl*
			pow(log10(M/M_FGB(met))/log10(M_HeF(met)/M_FGB(met)), b(46, met));
		return tau_bl;}
	
	if(M > M_FGB(met)) {
		tau_bl = (1 - b(47, met))*f_blM/f_blM_FGB;
		return tau_bl;}}
double Coef::t_DU(double M, double met){
	double L_x, L_DU, p, q, B, D, D0, zeta, x, y,
			 A_He, tinf1, tinf2, tx, t_DU, M_cDU, M_cSN, M_Ch;
	A_He = 7.66e-5;
	M_Ch = 1.44;

	if(M <= M_HeF(met)) p = 6.;
	if(M >= 2.5) 		p = 5.;
	if((M > M_HeF(met)) && (M < 2.5))
						p = 6. - (M - M_HeF(met))/(2.5 - M_HeF(met));

	if(M <= M_HeF(met)) q = 3.;
	if(M >= 2.5) 		q = 2.;
	if((M > M_HeF(met)) && (M < 2.5))
	 					q = 3. - (M - M_HeF(met))/(2.5 - M_HeF(met));
	
	zeta = log10(met/0.02);	
	B = max(3*1e4, 500 + 1.75*1e4*pow(M, 0.6));
	D0 = 5.37 + 0.135*zeta;
	if(M <= M_HeF(met))	D = D0;
	if(M >= 2.5)		D = max(-1.0, max(0.975*D0 - 0.18*M, 0.5*D0 - 0.06*M));
	if((M > M_HeF(met)) && (M < 2.5))
			{x = max(-1.0, max(0.975*D0 - 0.18*2.5, 0.5*D0 - 0.06*2.5)) - D0;
			 y = 2.5 - M_HeF(met);
			 D = D0 + (M - M_HeF(met))*x/y;}
	D = pow(10, D);	
	M_cSN = max(M_Ch, 0.773*M_cBAGB(M, met) - 0.35);

	if(M_cBAGB(M, met) <= 0.8) 								{M_cDU = M_cBAGB(M, met);}
	if((0.8 < M_cBAGB(M, met)) && (M_cBAGB(M, met)) < 2.25) {M_cDU = 0.44*M_cBAGB(M, met) + 0.448;}
	if(2.25 <= M_cBAGB(M, met))								{M_cDU = M_cSN;}

	L_DU 	= min(B*pow(M_cDU, q), D*pow(M_cDU, p));
	L_x 	= B*pow(B/D, q/(p - q));
	tinf1 	= t_HeI(M, met) + t_He(M, met) + pow(D/L_BAGB(M, met), (p - 1.)/p) / ((p - 1.)*A_He*D);
	tx 		= tinf1 - (tinf1 - t_HeI(M, met) - t_He(M, met))*pow(L_BAGB(M, met)/L_x, (p - 1.)/p);
	tinf2 	= tx + pow(B/L_x, (q - 1.)/q) / ((q - 1.)*A_He*B);
	if(L_DU <= L_x) 
		{t_DU = tinf1 - pow(D/L_DU, (p - 1.)/p)/((p - 1.)*A_He*D);}
	if(L_DU > L_x) 
		{t_DU = tinf2 - pow(B/L_DU, (q - 1.)/q)/((q - 1.)*A_He*B);}
	return t_DU;}
double Coef::t_AGB(double M, double met){
	double L_x, L_DU, L_SN, p, q, B, D, D0, zeta, x, y, lambda,
			 A_HHe, tinf1, tinf2, tx, t_AGB, M_cDU, M_cSN, M_cSN_real, M_Ch;
	if(M_cBAGB(M, met) >= 2.25) return t_DU(M, met);
	A_HHe = 1.27e-5;
	M_Ch = 1.44;
	lambda = min(0.9, 0.3 + 0.001*pow(M, 5));

	if(M <= M_HeF(met)) p = 6.;
	if(M >= 2.5) 		p = 5.;
	if((M > M_HeF(met)) && (M < 2.5))
						p = 6. - (M - M_HeF(met))/(2.5 - M_HeF(met));

	if(M <= M_HeF(met)) q = 3.;
	if(M >= 2.5) 		q = 2.;
	if((M > M_HeF(met)) && (M < 2.5))
	 					q = 3. - (M - M_HeF(met))/(2.5 - M_HeF(met));
	
	zeta = log10(met/0.02);	
	B = max(3*1e4, 500 + 1.75*1e4*pow(M, 0.6));
	D0 = 5.37 + 0.135*zeta;
	if(M <= M_HeF(met))	D = D0;
	if(M >= 2.5)		D = max(-1.0, max(0.975*D0 - 0.18*M, 0.5*D0 - 0.06*M));
	if((M > M_HeF(met)) && (M < 2.5))
			{x = max(-1.0, max(0.975*D0 - 0.18*2.5, 0.5*D0 - 0.06*2.5)) - D0;
			 y = 2.5 - M_HeF(met);
			 D = D0 + (M - M_HeF(met))*x/y;}
	D = pow(10, D);	
	M_cSN = max(M_Ch, 0.773*M_cBAGB(M, met) - 0.35);

	if(M_cBAGB(M, met) <= 0.8) 								{M_cDU = M_cBAGB(M, met);}
	if((0.8 < M_cBAGB(M, met)) && (M_cBAGB(M, met)) < 2.25) {M_cDU = 0.44*M_cBAGB(M, met) + 0.448;}
	M_cSN_real = M_cDU + (M_cSN - M_cDU)/(1. - lambda);

	L_SN	= min(B*pow(M_cSN_real, q), D*pow(M_cSN_real, p));
	L_DU 	= min(B*pow(M_cDU, q), D*pow(M_cDU, p));
	L_x 	= B*pow(B/D, q/(p - q));
	if(L_DU >= L_x){
		tinf2 = t_DU(M, met) + pow(B/L_DU, (q - 1.)/q);
		t_AGB = tinf2 - pow(B/L_SN, (q - 1.)/q)/((q - 1.)*A_HHe*B);
	}
	else{
		tinf1 	= t_DU(M, met) + pow(D/L_DU, (p - 1.)/p) / ((p - 1.)*A_HHe*D);
		tx 		= tinf1 - (tinf1 - t_DU(M, met))*pow(L_DU/L_x, (p - 1.)/p);
		tinf2 	= tx + pow(B/L_x, (q - 1.)/q) / ((q - 1.)*A_HHe*B);
		if(L_SN <= L_x) 
			{t_AGB = tinf1 - pow(D/L_SN, (p - 1.)/p)/((p - 1.)*A_HHe*D);}
		if(L_SN > L_x) 
			{t_AGB = tinf2 - pow(B/L_SN, (q - 1.)/q)/((q - 1.)*A_HHe*B);}
	}
	return t_AGB;}

double Coef::M_cBGB(double M, double met){
	double L_x, M_cHeI, p, q, B, D, zeta, x, y, A_H, tinf1, tinf2, tx;
	double c1, c2, C;

	if(M <= M_HeF(met)) {	
		A_H = max(-4.8, min(-5.7 + 0.8*M, -4.1 + 0.14*M));
		A_H = pow(10, A_H);
		p = 6;
		q = 3;
		zeta = log10(met/0.02);	
		B = max(3*1e4, 500 + 1.75*1e4*pow(M, 0.6));
		D = 5.37 + 0.135*zeta;
		D = pow(10, D);
		L_x 	= B*pow(B/D, q/(p - q));
		tinf1 	= t_BGB(M, met) + pow(D/L_BGB(M, met), (p - 1)/p) / ((p - 1)*A_H*D);
		tx 		= tinf1 - (tinf1 - t_BGB(M, met))*pow(L_BGB(M, met)/L_x, (p - 1)/p);
		tinf2 	= tx + pow(B/L_x, (q - 1)/q) / ((q - 1)*A_H*B);
		if(L_HeI(M, met) <= L_x) {
			M_cHeI = pow((p - 1)*A_H*D*(tinf1 - t_BGB(M, met)), 1/(1 - p));}
	}

	if(M > M_HeF(met)) {
		c1 = 9.20925e-5;
		c2 = 5.402216;
		C = pow(Coef::M_cBGB(M_HeF(met), met), 4) - c1*pow(M_HeF(met), c2);
		M_cHeI = min(0.95*M_cBAGB(M, met), pow(C + c1*pow(M, c2), 0.25));
	} 

	return M_cHeI;}
double Coef::M_cHeI(double M, double met){
	double L_x, M_cHeI, p, q, B, D, zeta, x, y, A_H, tinf1, tinf2, tx;
	double c1, c2, C;

	if(M <= M_HeF(met)) {	
		A_H = max(-4.8, min(-5.7 + 0.8*M, -4.1 + 0.14*M));
		A_H = pow(10, A_H);
		p = 6;
		q = 3;
		zeta = log10(met/0.02);	
		B = max(3*1e4, 500 + 1.75*1e4*pow(M, 0.6));
		D = 5.37 + 0.135*zeta;
		D = pow(10, D);
		L_x 	= B*pow(B/D, q/(p - q));
		tinf1 	= t_BGB(M, met) + pow(D/L_BGB(M, met), (p - 1)/p) / ((p - 1)*A_H*D);
		tx 		= tinf1 - (tinf1 - t_BGB(M, met))*pow(L_BGB(M, met)/L_x, (p - 1)/p);
		tinf2 	= tx + pow(B/L_x, (q - 1)/q) / ((q - 1)*A_H*B);
		if(L_HeI(M, met) <= L_x) {
			M_cHeI = pow((p - 1)*A_H*D*(tinf1 - t_HeI(M, met)), 1/(1 - p));}
		if(L_HeI(M, met) > L_x) {
			M_cHeI = pow((q - 1)*A_H*B*(tinf2 - t_HeI(M, met)), 1/(1 - q));}
	}

	if(M > M_HeF(met)) {
		c1 = 9.20925e-5;
		c2 = 5.402216;
		C = pow(Coef::M_cHeI(M_HeF(met), met), 4) - c1*pow(M_HeF(met), c2);
		M_cHeI = min(0.95*M_cBAGB(M, met), pow(C + c1*pow(M, c2), 0.25));
	} 

	return M_cHeI;}
double Coef::M_cBAGB(double M, double met){
	double x;
	x = pow(b(36, met)*pow(M, b(37, met)) + b(38, met), 0.25);
	return x;}
double Coef::M_cGB(double t, double M, double met){
	double L_x, M_cGB, p, q, B, D, zeta, x, y, A_H, tinf1, tinf2, tx;
	double c1, c2, C, M_cBAGB, tau;

	if(M <= M_HeF(met)) {	
		A_H = max(-4.8, min(-5.7 + 0.8*M, -4.1 + 0.14*M));
		A_H = pow(10, A_H);
		p = 6;
		q = 3;
		zeta = log10(met/0.02);	
		B = max(3*1e4, 500 + 1.75*1e4*pow(M, 0.6));
		D = 5.37 + 0.135*zeta;
		D = pow(10, D);

		L_x 	= B*pow(B/D, q/(p - q));
		tinf1 	= t_BGB(M, met) + pow(D/L_BGB(M, met), (p - 1)/p) / ((p - 1)*A_H*D);
		tx 		= tinf1 - (tinf1 - t_BGB(M, met))*pow(L_BGB(M, met)/L_x, (p - 1)/p);
		tinf2 	= tx + pow(B/L_x, (q - 1)/q) / ((q - 1)*A_H*B);
		if(L_HeI(M, met) <= L_x) {
			M_cGB = pow((p - 1)*A_H*D*(tinf1 - t), 1/(1 - p));}
		if(L_HeI(M, met) > L_x) {
			M_cGB = pow((q - 1)*A_H*B*(tinf2 - t), 1/(1 - q));}
	}

	if(M > M_HeF(met)) {
		tau = (t - t_BGB(M, met))/(t_HeI(M, met) - t_BGB(M, met));
		M_cGB = M_cBGB(M, met) + (M_cHeI(M, met) - M_cBGB(M, met))*tau;
	}
	return M_cGB;}
double Coef::M_c(double t, double M, double met){
	double tau, M_c;
	tau = (t - t_HeI(M, met))/t_He(M, met);
	M_c = (1 - tau)*M_cHeI(M, met) + tau*M_cBAGB(M, met);
	return M_c;}

double Coef::M_HeF(double met){
	double x;
	x = 1.995 + 0.25*log10(met/0.02) + 0.087*pow(log10(met/0.02), 2);
	return x;}
double Coef::M_FGB(double met){
	double x, y;
	x = 13.048*pow(met/0.02, 0.06);
	y = 1 + 0.0012*pow(0.02/met, 1.27);
	return x/y;}
double Coef::L_BGB(double M, double met){
	double c2, c3, x, y;
	c2 = 9.301992;
	c3 = 4.637345;
	x = a(27, met)*pow(M, a(31, met)) + a(28, met)*pow(M, c2);
	y = a(29, met) + a(30, met)*pow(M, c3) + pow(M, a(32, met));
	return x/y;}

double Coef::L_ZAMS(double M, double met){
	double x, y;
	x = a(82, met)*pow(M, 5.5) + a(83, met)*pow(M, 11);
	y = a(84, met) + pow(M, 3) + a(85, met)*pow(M, 5) +
		a(86, met)*pow(M, 7) + a(87, met)*pow(M, 8) +
		a(88, met)*pow(M, 9.5);
	return x/y;}
double Coef::R_ZAMS(double M, double met){
	double x, y;
	x = a(89, met)*pow(M, 2.5) + a(90, met)*pow(M, 6.5) +
		a(91, met)*pow(M, 11) + a(92, met)*pow(M, 19) +
		a(93, met)*pow(M, 19.5);
	y = a(94, met) + a(95, met)*pow(M, 2) + 
		a(96, met)*pow(M, 8.5) + pow(M, 18.5) +
		a(97, met)*pow(M, 19.5);
	return x/y;}

double Coef::L_TMS(double M, double met){
	double x, y;
	x = a(11, met)*pow(M, 3) + a(12, met)*pow(M, 4) + 
		a(13, met)*pow(M, a(16, met) + 1.8);
	y = a(14, met) + a(15, met)*pow(M, 5) + pow(M, a(16, met));
	return x/y;}
double Coef::R_TMS(double M, double met){
	double x, y, R_TMS;
	if(M <= a(17, met)){
		x = a(18, met) + a(19, met)*pow(M, a(21, met));
		y = a(20, met) + pow(M, a(22, met));
		R_TMS = x/y;}
	if(M >= (a(17, met) + 0.1)){
		x = (-8.672073e-2)*pow(M, 3) + a(23, met)*pow(M, a(26, met)) + 
			a(24, met)*pow(M, a(26, met) + 1.5);
		y = a(25, met) + pow(M, 5);
		R_TMS = x/y;}
	if((M > a(17, met)) && (M < (a(17, met) + 0.1))){
		x = Coef::R_TMS(a(17, met), met);
		y = Coef::R_TMS(a(17, met) + 0.1, met);
		R_TMS = x + (y - x)/0.1*(M - a(17, met));}
	if(M < 0.5) {R_TMS = max(R_TMS, 1.5*R_ZAMS(M, met));}
	return R_TMS;}

double Coef::L_EHG(double M, double met){
	if(M < M_FGB(met)) 	{return L_BGB(M, met);}
	else 				{return L_HeI(M, met);}}
double Coef::R_EHG(double M, double met){
	if(M < M_FGB(met)) 	return R_GBadd(L_BGB(M, met), M, met);
	else				return R_HeI(M, met);}

double Coef::L_ZHe(double M){
	double x, y;
	x = 15262*pow(M, 10.25);
	y = pow(M, 9) + 29.54*pow(M, 7.5) + 31.18*pow(M, 6) + 0.0469;
	return x/y;}
double Coef::R_ZHe(double M){
	double x, y;
	x = 0.2391*pow(M, 4.6);
	y = pow(M, 4) + 0.162*pow(M, 3) + 0.0065;
	return x/y;}

double Coef::L_minHe(double M, double met){
	double x, y, c;
	c = b(17, met)/pow(M_FGB(met), 0.1) + 
		(b(16, met)*b(17, met) - b(14, met))/pow(M_FGB(met), b(15, met) + 0.1);
	x = L_HeI(M, met)*(b(14, met) + c*pow(M, b(15, met) + 0.1));
	y = b(16, met) + pow(M, b(15, met));
	return x/y;}
double Coef::R_mHe(double M, double met){
	double x, y, mu, M_c, R_mHe;
	if(M >= M_HeF(met)) {
		x = b(24, met)*M + pow(b(25, met)*M, b(26, met))*pow(M, b(28, met));
		y = b(27,  met) + pow(M, b(28, met));
		R_mHe = x/y;
		//следующие три строки моя личная вставка
		if(R_mHe > R_AGBadd(L_HeI(M, met), M, met)){
			R_mHe = R_AGBadd(L_HeI(M, met), M, met);
		}
		return R_mHe;}
	if(M < M_HeF(met)) {
		M_c = M_cHeI(M, met);
		mu = (M - M_c)/(M_HeF(met) - M_c);
		R_mHe = R_GBadd(L_ZAHB(M, met), M, met)*
				pow(Coef::R_mHe(M_HeF(met), met)/
				R_GBadd(L_ZAHB(M_HeF(met), met), M_HeF(met), met), mu);
		return R_mHe;}}		

double Coef::L_HeI(double M, double met){
	double L_HeI, L_HeIM, alpha1;
	L_HeIM = (b(11, met) + b(12, met)*pow(M_HeF(met), 3.8))/
			 (b(13, met) + pow(M_HeF(met), 2));
	alpha1 = (b(9, met)*pow(M_HeF(met), b(10, met)) - L_HeIM)/L_HeIM;
	
	if(M >= M_HeF(met)) {L_HeI = 
		(b(11, met) + b(12, met)*pow(M, 3.8))/(b(13, met) + pow(M, 2));}
	else 				{L_HeI = 
		b(9, met)*pow(M, b(10, met))/(1 + alpha1*pow(15,(M - M_HeF(met))));}
	return L_HeI;}
double Coef::R_HeI(double M, double met){
	double mu;
	if(M < M_FGB(met)) 					{return R_GB(L_HeI(M, met), M, met);}
	if(M >= max(M_FGB(met), 12.0)) 		{return R_mHe(M, met);}
	if((M > M_FGB(met)) && (M < 12.0))	{
		mu = log10(M/12.0)/log10(M_FGB(met)/12.0);
		return R_mHe(M, met)*pow(R_GB(L_HeI(M, met), M, met)/R_mHe(M, met), mu);}}

double Coef::L_ZAHB(double M, double met){
	double mu, alpha2, M_c, L_ZAHB, x, y;
	M_c = M_cHeI(M, met);
	mu = (M - M_c)/(M_HeF(met) - M_c);
	alpha2 = (b(18, met) + L_ZHe(M_c) - L_minHe(M_HeF(met), met))/
			 (L_minHe(M_HeF(met), met) - L_ZHe(M_c));
	x = (1 + b(20, met))/(1 + b(20, met)*pow(mu, 1.6479));
	y = b(18, met)*pow(mu, b(19, met))/(1 + alpha2*exp(15*(M - M_HeF(met))));
	L_ZAHB = L_ZHe(M_c) + x*y;
	return L_ZAHB;}
double Coef::R_ZAHB(double M, double met){
	double f, R_ZAHB, M_c, mu;
	M_c = M_cHeI(M, met);
	mu = (M - M_c)/(M_HeF(met) - M_c);
	f = (1 + b(21, met))*pow(mu, b(22, met))/(1 + b(21, met)*pow(mu, b(23, met)));
	R_ZAHB = (1 - f)*R_ZHe(M_c) + f*R_GBadd(L_ZAHB(M, met), M, met);
	return R_ZAHB;}

double Coef::L_BAGB(double M, double met){
	double L_BAGB, x, y, alpha3;
	
	if(M >= M_HeF(met)) {
		x = b(31, met) + b(32, met)*pow(M, b(33, met) + 1.8);
		y = b(34, met) + pow(M, b(33, met));
		L_BAGB = x/y;
		return L_BAGB;}
	
	if(M < M_HeF(met)) {
		alpha3 = (b(29, met)*pow(M_HeF(met), b(30, met)) - Coef::L_BAGB(M_HeF(met), met))/
		 			Coef::L_BAGB(M_HeF(met), met);
		x = b(29, met)*pow(M, b(30, met));
		y = 1 + alpha3*exp(15*(M - M_HeF(met)));
		L_BAGB = x/y; 
		return L_BAGB;}}
double Coef::R_BAGB(double M, double met){
	double x;
	x = R_AGBadd(L_BAGB(M, met), M, met);
	return x;}

double Coef::R_GBadd(double L, double M, double met){
	double A, x;
	A = min(b(4, met)*pow(M, -b(5, met)), b(6, met)*pow(M, -b(7, met)));
	x = A*(pow(L, b(1, met)) + b(2, met)*pow(L, b(3, met)));
	return x;}
double Coef::R_AGBadd(double L, double M, double met){
	double A, R_AGBadd;

	if(M >= M_HeF(met)) {
		A = min(b(51, met)*pow(M, -b(52, met)), b(53, met)*pow(M, -b(54, met)));
		R_AGBadd = A*(pow(L, b(1, met)) + b(2, met)*pow(L, b(55, met)*b(3, met)));
		return R_AGBadd;}
	if(M <= (M_HeF(met) - 0.2)) {
		A = b(56, met) + b(57, met)*M;
		R_AGBadd = A*(pow(L, b(1, met)) + b(2, met)*pow(L, b(3, met)));
		return R_AGBadd;}
	if((M < M_HeF(met)) && (M > (M_HeF(met) - 0.2))) {
		R_AGBadd = Coef::R_AGBadd(L, M_HeF(met), met) + 
		(Coef::R_AGBadd(L, M_HeF(met), met) - Coef::R_AGBadd(L, M_HeF(met) - 0.2, met))/
		0.2*(M - M_HeF(met) + 0.2);
		return R_AGBadd;}}

double Coef::L_MS(double t, double M, double met){
	double tau, tau1, tau2;
	double zeta, B, alpha_L, beta_L, eps, eta; 
	double L_MS, dL, M_hook;
	zeta = log10(met/0.02);
	eps = 0.01;
	tau = t/t_MS(M, met);
	tau1 = min(1.0, t/t_hook(M, met));
	tau2 = max(0.0, min(1.0, (t - (1.0 - eps)*t_hook(M, met))/
			(eps*t_hook(M, met))));
	M_hook = 1.0185 + 0.16015*zeta + 0.0892*pow(zeta, 2);

	B = (a(45, met) + a(46, met)*pow(2.0, a(48, met)))/
		(pow(2.0, 0.4) + a(47, met)*pow(2.0, 1.9));
	if(M >= 2.0) 								{alpha_L = 
		(a(45, met) + a(46, met)*pow(M, a(48, met)))/
		(pow(M, 0.4) + a(47, met)*pow(M, 1.9));}
	if(M < 0.5) 								{alpha_L = 
		a(49, met);}	
	if((M >= 0.5) && (M < 0.7)) 				{alpha_L =
	 	a(49, met) + 5*(0.3 - a(49, met))*(M - 0.5);}	
	if((M >= 0.7) && (M < a(52, met))) 			{alpha_L =
		 0.3 + (a(50, met) - 0.3)*(M - 0.7)/(a(52, met) - 0.7);}	
	if((M >= a(52, met)) && (M < a(53, met))) 	{alpha_L =
	 	a(50, met) + (a(51, met) - a(50, met))*
	 	(M - a(52, met))/(a(53, met) - a(52, met));}	
	if((M >= a(53, met)) && (M < 2.0)) 			{alpha_L =
	 	a(51, met) + (B - a(51, met))*
	 	(M - a(53, met))/(2 - a(53, met));}

	B = max(0.0, a(54, met) - a(55, met)*pow(a(57, met), a(56, met))); 
	beta_L = max(0.0, a(54, met) - a(55, met)*pow(M, a(56, met)));
	if((M > a(57, met)) && (beta_L > 0.0)) {beta_L =
		max(0.0, B - 10.0*(M - a(57, met))*B);}

	B = min(a(34, met)/pow(a(33, met), a(35, met)),
			a(36, met)/pow(a(33, met), a(37, met)));
	if(M <= M_hook) 						{dL =
		0;}
	if((M_hook < M) && (M < a(33, met))) 	{dL =
		B*pow((M - M_hook)/(a(33, met) - M_hook), 0.4);}
	if(M >= a(33, met)) 					{dL = 
		min(a(34, met)/pow(M, a(35, met)), a(36, met)/pow(M, a(37, met)));}
	
	if(met > 0.0009) {eta = 10;}
	else {
		if(M <= 1.0) {eta = 10;}
		if(M >= 1.1) {eta = 20;}
		if((M > 1.0) && (M < 1.1)) {eta = 10 + 100*(M - 1.0);}}

	L_MS = alpha_L*tau + beta_L*pow(tau, eta) + 
			(log10(L_TMS(M, met)/L_ZAMS(M, met)) - alpha_L - beta_L)*
			pow(tau, 2) - dL*(pow(tau1, 2) - pow(tau2, 2));	

	// cout 	<< "alphaL = " << alpha_L << '\n'
	// 		<< "betaL = " << beta_L << '\n'
	// 		<< endl;
	// cin.get();		
	return L_MS + log10(L_ZAMS(M, met));}
double Coef::R_MS(double t, double M, double met){
	double tau, tau1, tau2;
	double zeta, B, C, alpha_R, beta_R, gamma, eps; 
	double R_MS, dR, M_hook;
	zeta = log10(met/0.02);
	eps = 0.01;
	tau = t/t_MS(M, met);
	tau1 = min(1.0, t/t_hook(M, met));
	tau2 = max(0.0, min(1.0, (t - (1.0 - eps)*t_hook(M, met))/
			(eps*t_hook(M, met)) ));
	M_hook = 1.0185 + 0.16015*zeta + 0.0892*pow(zeta, 2);

	B = (a(58, met)*pow(a(66, met), a(60, met)))/
		(a(59, met)*pow(a(66, met), a(61, met)));
	C = (a(58, met)*pow(a(67, met), a(60, met)))/
		(a(59, met)*pow(a(67, met), a(61, met)));
	if(M < 0.5) 								{alpha_R = 
		a(62, met);}
	if((M >= 0.5) && (M < 0.65)) 				{alpha_R = 
		a(62, met) + (a(63, met) - a(62, met))*(M - 0.5)/0.15;}	
	if((M >= 0.65) && (M < a(68, met))) 		{alpha_R = 
		a(63, met) + (a(64, met) - a(63, met))*(M - 0.65)/(a(68, met) - 0.65);}
	if((M >= a(68, met)) && (M < a(66, met))) 	{alpha_R =
		a(64, met) + (B - a(64, met))*(M - a(68, met))/(a(66, met) - a(68, met));}	
	if((M >= a(66, met)) && (M <= a(67, met))) 	{alpha_R = 
		(a(58, met)*pow(M, a(60, met)))/(a(59, met)*pow(M, a(61, met)));}	
	if(M > a(67, met)) 							{alpha_R = 
		C + a(65, met)*(M - a(67, met));}


	B = a(69, met)*pow(2, 3.5)/(a(70, met) + pow(2, a(71, met)));
	C = a(69, met)*pow(16.0, 3.5)/(a(70, met) + pow(16.0, a(71, met)));
	if(M <= 1.0) 						{beta_R = 
		1.06;}	
	if((M > 1.0) && (M < a(74, met))) 	{beta_R =
		1.06 + (a(72, met) - 1.06)*(M - 1.0)/(a(74, met) - 1);}	
	if((M >= a(74, met)) && (M < 2.0)) 	{beta_R =
		a(72, met) + (B - a(72, met))*(M - a(74, met))/(2.0 - a(74, met));}	
	if((M >= 2.0) && (M <= 16.0)) 		{beta_R =
		a(69, met)*pow(M, 3.5)/(a(70, met) + pow(M, a(71, met)));}	
	if(M > 16.0) 						{beta_R = 
		C + a(73, met)*(M - 16.0);}
	beta_R = beta_R - 1;

	
	B = a(76, met) + a(77, met)*pow(1.0 - a(78, met), a(79, met));
	if(a(75, met) > 1.0) {C = a(80, met);}
	else {C = B;}	
	
	if(M <= 1.0) 									{gamma =
		a(76, met) + a(77, met)*pow(M - a(78, met), a(79, met));}
	if((M > 1.0) && (M <= a(75, met))) 				{gamma =
		B + (a(80, met) - B)*pow((M - 1.0)/(a(75, met) - 1.0), a(81, met));}
	if((M > a(75, met)) && (M < (a(75, met) + 0.1))){gamma =
		C - 10.0*(M - a(75, met))*C;}
	if(M >= (a(75, met) + 0.1)) 					{gamma =
		0.0;}
	if(gamma < 0) {gamma = 0;}
	

	B = (a(38, met) + a(39, met)*pow(2.0, 3.5))/
		(a(40, met)*pow(2.0, 3) + pow(2.0, a(41, met))) - 1;
	if(M <= M_hook) 						{dR = 
		0.0;}
	if((M > M_hook) && (M <= a(42, met))) 	{dR =
		a(43, met)*pow((M - M_hook)/(a(42, met) - M_hook), 0.5);}
	if((M > a(42, met)) && (M < 2.0)) 		{dR =
		a(43, met) + (B - a(43, met))*
		pow((M - a(42, met))/(2.0 - a(42, met)), a(44, met));}
	if(M >= 2.0) 							{dR = 
		(a(38, met) + a(39, met)*pow(M, 3.5))/
		(a(40, met)*pow(M, 3) + pow(M, a(41, met))) - 1;}

	R_MS = alpha_R*tau + beta_R*pow(tau, 10) + gamma*pow(tau, 40) + 
			(log10(R_TMS(M, met)/R_ZAMS(M, met)) - alpha_R - beta_R - gamma)*
			pow(tau, 3) - dR*(pow(tau1, 3) - pow(tau2, 3));			
	return R_MS + log10(R_ZAMS(M, met));}
double Coef::L_HG(double t, double M, double met){
	double tau, L_HG;
	tau = (t - t_MS(M, met))/(t_BGB(M, met) - t_MS(M, met));
	L_HG = L_TMS(M, met)*pow(L_EHG(M, met)/L_TMS(M, met), tau);
	return log10(L_HG);}
double Coef::R_HG(double t, double M, double met){
	double tau, R_HG;
	tau = (t - t_MS(M, met))/(t_BGB(M, met) - t_MS(M, met));
	R_HG = R_TMS(M, met)*pow(R_EHG(M, met)/R_TMS(M, met), tau);
	return log10(R_HG);}
double Coef::L_GB(double t, double M, double met){
	double L_GB, L_x, M_x, M_cGB, p, q, B, D, D0, zeta, x, y, A_H, tinf1, tinf2, tx;
	ofstream fout;
	A_H = max(-4.8, min(-5.7 + 0.8*M, -4.1 + 0.14*M));
	A_H = pow(10, A_H);

	if(M <= M_HeF(met)) {p = 6;}
	if(M >= 2.5) 		{p = 5;}
	if((M > M_HeF(met)) && (M < 2.5))
						{p = 6 - (M - M_HeF(met))/(2.5 - M_HeF(met));}

	if(M <= M_HeF(met)) {q = 3;}
	if(M >= 2.5) 		{q = 2;}
	if((M > M_HeF(met)) && (M < 2.5))
	 					{q = 3 - (M - M_HeF(met))/(2.5 - M_HeF(met));}
	
	zeta = log10(met/0.02);	
	B = max(3*1e4, 500 + 1.75*1e4*pow(M, 0.6));
	D0 = 5.37 + 0.135*zeta;
	if(M <= M_HeF(met))	{D = D0;}
	if(M >= 2.5)		{D = max(-1.0, max(0.975*D0 - 0.18*M, 0.5*D0 - 0.06*M));}
	if((M > M_HeF(met)) && (M < 2.5))
			{x = max(-1.0, max(0.975*D0 - 0.18*2.5, 0.5*D0 - 0.06*2.5)) - D0;
			 y = 2.5 - M_HeF(met);
			 D = D0 + (M - M_HeF(met))*x/y;}
	D = pow(10, D);

	L_x 	= B*pow(B/D, q/(p - q));
	M_x 	= pow(B/D, 1/(p - q));
	tinf1 	= t_BGB(M, met) + pow(D/L_BGB(M, met), (p - 1)/p) / ((p - 1)*A_H*D);
	tx 		= tinf1 - pow(D/L_x, (p - 1)/p) / ((p - 1)*A_H*D);
	tinf2 	= tx + pow(B/L_x, (q - 1)/q) / ((q - 1)*A_H*B);


	if((t <= tx) && (t >= t_BGB(M, met))) 
		{L_GB = D*pow((p - 1)*A_H*D*(tinf1 - t), p/(1 - p));}
	if((t > tx) && (t <= t_HeI(M, met))) 
		{L_GB = B*pow((q - 1)*A_H*B*(tinf2 - t), q/(1 - q));}	

	return L_GB;}
double Coef::R_GB(double t, double M, double met){
	double A, L, x;
	L = L_GB(t, M, met);
	A = min(b(4, met)*pow(M, -b(5, met)), b(6, met)*pow(M, -b(7, met)));
	x = A*(pow(L, b(1, met)) + b(2, met)*pow(L, b(3, met)));
	return x;}
double Coef::L_CHeB(double t, double M, double met){
	double L_x, R_x, tau, tau_x, lambda, lambda_sh, xi, L;

	if(M < M_HeF(met)) {L_x = L_ZAHB(M, met);}
	if((M >= M_HeF(met)) && M < M_FGB(met)) {L_x = L_minHe(M, met);}
	if(M >= M_FGB(met)) {L_x = L_HeI(M, met);}

	if(M < M_HeF(met)) {R_x = R_ZAHB(M, met);}
	if((M >= M_HeF(met)) && (M < M_FGB(met))) {R_x = R_GBadd(L_minHe(M, met), M, met);}
	if(M >= M_FGB(met)) {R_x = R_HeI(M, met);}

	tau = (t - t_HeI(M, met))/t_He(M, met);
	if((M > M_HeF(met)) && (M < M_FGB(met))) {tau_x = 1 - tau_bl(M, met);}
	if((M <= M_HeF(met)) || (M >= M_FGB(met))) {tau_x = 0;}
	
	xi = min(2.5, max(0.4, R_mHe(M, met)/R_x));
	lambda = pow((tau - tau_x)/(1 - tau_x), xi);
	lambda_sh = pow((tau_x - tau)/tau_x, 3);

	if(tau > 1) {tau = 1;}
	if(tau < 0) {tau = 0;}
	if((tau_x <= tau) && (tau <= 1)) {L = L_x*pow(L_BAGB(M, met)/L_x, lambda);}
	if((tau_x > tau) && (tau >= 0)) {L = L_x*pow(L_HeI(M, met)/L_x, lambda_sh);}

	return L;}
double Coef::R_CHeB(double t, double M, double met){
	double R, R_x, R_min, R_y, L_y, tau, tau_x, tau_y, ro;

	if(M < M_HeF(met)) {R_x = R_ZAHB(M, met);}
	if((M >= M_HeF(met)) && (M < M_FGB(met))) {R_x = R_GBadd(L_minHe(M, met), M, met);}
	if(M >= M_FGB(met)) {R_x = R_HeI(M, met);}

	R_min = min(R_mHe(M, met), R_x);
	tau = (t - t_HeI(M, met))/t_He(M, met);
	
	if((M > M_HeF(met)) && (M < M_FGB(met))) 	{tau_x = 1 - tau_bl(M, met);}
	if((M <= M_HeF(met)) || (M >= M_FGB(met))) 	{tau_x = 0;}
	if(M < M_FGB(met)) 	{tau_y = 1;}
	if(M >= M_FGB(met)) {tau_y = tau_bl(M, met);}
	
	if(M <= M_FGB(met)) {L_y = L_BAGB(M, met);}
	if(M > M_FGB(met)) 	{L_y = L_CHeB(t_HeI(M, met) + t_He(M, met)*tau_y, M, met);} 	 

	R_y = R_AGBadd(L_y, M, met);
	if(tau > 1) {tau = 1;}
	if(tau < 0) {tau = 0;}

	if((0 <= tau) && (tau < tau_x)) {R = R_GBadd(L_CHeB(t, M, met), M, met);}
	if((tau_y <= tau) && (tau <= 1)) {R = R_AGBadd(L_CHeB(t, M, met), M, met);}
	if((tau_x <= tau) && (tau < tau_y)) {
		ro = pow(log(R_y/R_min), 1./3.)*(tau - tau_x)/(tau_y - tau_x) -
		 	 pow(log(R_x/R_min), 1./3.)*(tau_y - tau)/(tau_y - tau_x);
		R = R_min*exp(pow(fabs(ro), 3));}

	return R;}
double Coef::L_EAGB(double t, double M, double met){
	double L_x, L_DU, L_EAGB, p, q, B, D, D0, zeta, x, y,
			 A_He, tinf1, tinf2, tx, t_DU, M_cDU, M_cSN, M_Ch;
	A_He = 7.66e-5;
	M_Ch = 1.44;

	if(M <= M_HeF(met)) {p = 6;}
	if(M >= 2.5) 		{p = 5;}
	if((M > M_HeF(met)) && (M < 2.5))
						{p = 6 - (M - M_HeF(met))/(2.5 - M_HeF(met));}

	if(M <= M_HeF(met)) {q = 3;}
	if(M >= 2.5) 		{q = 2;}
	if((M > M_HeF(met)) && (M < 2.5))
	 					{q = 3 - (M - M_HeF(met))/(2.5 - M_HeF(met));}
	
	zeta = log10(met/0.02);	
	B = max(3*1e4, 500 + 1.75*1e4*pow(M, 0.6));
	D0 = 5.37 + 0.135*zeta;
	if(M <= M_HeF(met))	{D = D0;}
	if(M >= 2.5)		{D = max(-1.0, max(0.975*D0 - 0.18*M, 0.5*D0 - 0.06*M));}
	if((M > M_HeF(met)) && (M < 2.5))
			{x = max(-1.0, max(0.975*D0 - 0.18*2.5, 0.5*D0 - 0.06*2.5)) - D0;
			 y = 2.5 - M_HeF(met);
			 D = D0 + (M - M_HeF(met))*x/y;}
	D = pow(10, D);	
	M_cSN = max(M_Ch, 0.773*M_cBAGB(M, met) - 0.35);

	if(M_cBAGB(M, met) <= 0.8) 								{M_cDU = M_cBAGB(M, met);}
	if((0.8 < M_cBAGB(M, met)) && (M_cBAGB(M, met)) < 2.25) {M_cDU = 0.44*M_cBAGB(M, met) + 0.448;}
	if(2.25 <= M_cBAGB(M, met))								{M_cDU = M_cSN;}

	L_DU 	= min(B*pow(M_cDU, q), D*pow(M_cDU, p));
	L_x 	= B*pow(B/D, q/(p - q));
	tinf1 	= t_HeI(M, met) + t_He(M, met) + pow(D/L_BAGB(M, met), (p - 1)/p) / ((p - 1)*A_He*D);
	tx 		= tinf1 - (tinf1 - t_HeI(M, met) - t_He(M, met))*pow(L_BAGB(M, met)/L_x, (p - 1)/p);
	tinf2 	= tx + pow(B/L_x, (q - 1)/q) / ((q - 1)*A_He*B);
	if(t <= tx) 
		L_EAGB = D*pow((p - 1)*A_He*D*(tinf1 - t), p/(1 - p));
	if(t > tx) 
		L_EAGB = B*pow((q - 1)*A_He*B*(tinf2 - t), q/(1 - q));
	return L_EAGB;}
double Coef::R_EAGB(double t, double M, double met){
	double A, R_AGB, L;
	L = L_EAGB(t, M, met);

	if(M >= M_HeF(met)) {
		A = min(b(51, met)*pow(M, -b(52, met)), b(53, met)*pow(M, -b(54, met)));
		R_AGB = A*(pow(L, b(1, met)) + b(2, met)*pow(L, b(55, met)*b(3, met)));
		return R_AGB;}
	if(M <= (M_HeF(met) - 0.2)) {
		A = b(56, met) + b(57, met)*M;
		R_AGB = A*(pow(L, b(1, met)) + b(2, met)*pow(L, b(3, met)));
		return R_AGB;}
	if((M < M_HeF(met)) && (M > (M_HeF(met) - 0.2))) {
		R_AGB = Coef::R_AGBadd(L, M_HeF(met), met) + 
		(Coef::R_AGBadd(L, M_HeF(met), met) - Coef::R_AGBadd(L, M_HeF(met) - 0.2, met))/
		0.2*(M - M_HeF(met) + 0.2);
		return R_AGB;}}
double Coef::L_TPAGB(double t, double M, double met){
	double L_x, L_DU, L_SN, L_AGB, p, q, B, D, D0, zeta, x, y, lambda,
			 A_HHe, tinf1, tinf2, tx, M_cDU, M_cSN, M_cSN_real, M_Ch;
	A_HHe = 1.27e-5;
	M_Ch = 1.44;
	lambda = min(0.9, 0.3 + 0.001*pow(M, 5));

	if(M <= M_HeF(met)) p = 6.;
	if(M >= 2.5) 		p = 5.;
	if((M > M_HeF(met)) && (M < 2.5))
						p = 6. - (M - M_HeF(met))/(2.5 - M_HeF(met));

	if(M <= M_HeF(met)) q = 3.;
	if(M >= 2.5) 		q = 2.;
	if((M > M_HeF(met)) && (M < 2.5))
	 					q = 3. - (M - M_HeF(met))/(2.5 - M_HeF(met));
	
	zeta = log10(met/0.02);	
	B = max(3*1e4, 500 + 1.75*1e4*pow(M, 0.6));
	D0 = 5.37 + 0.135*zeta;
	if(M <= M_HeF(met))	D = D0;
	if(M >= 2.5)		D = max(-1.0, max(0.975*D0 - 0.18*M, 0.5*D0 - 0.06*M));
	if((M > M_HeF(met)) && (M < 2.5))
			{x = max(-1.0, max(0.975*D0 - 0.18*2.5, 0.5*D0 - 0.06*2.5)) - D0;
			 y = 2.5 - M_HeF(met);
			 D = D0 + (M - M_HeF(met))*x/y;}
	D = pow(10, D);	
	M_cSN = max(M_Ch, 0.773*M_cBAGB(M, met) - 0.35);

	if(M_cBAGB(M, met) <= 0.8) 								{M_cDU = M_cBAGB(M, met);}
	if((0.8 < M_cBAGB(M, met)) && (M_cBAGB(M, met)) < 2.25) {M_cDU = 0.44*M_cBAGB(M, met) + 0.448;}
	M_cSN_real = M_cDU + (M_cSN - M_cDU)/(1. - lambda);

	L_SN	= min(B*pow(M_cSN_real, q), D*pow(M_cSN_real, p));
	L_DU 	= min(B*pow(M_cDU, q), D*pow(M_cDU, p));
	L_x 	= B*pow(B/D, q/(p - q));
	if(L_DU >= L_x){
		tinf2 = t_DU(M, met) + pow(B/L_DU, (q - 1.)/q);
		L_AGB = B*pow((q - 1.)*A_HHe*B*(tinf2 - t), q/(1. - q));
	}
	else{
		tinf1 	= t_DU(M, met) + pow(D/L_DU, (p - 1.)/p) / ((p - 1.)*A_HHe*D);
		tx 		= tinf1 - (tinf1 - t_DU(M, met))*pow(L_DU/L_x, (p - 1.)/p);
		tinf2 	= tx + pow(B/L_x, (q - 1.)/q) / ((q - 1.)*A_HHe*B);
		if(t <= tx) 
			L_AGB = D*pow((p - 1.)*A_HHe*B*(tinf1 - t), p/(1. - p));
		if(t > tx) 
			L_AGB = B*pow((q - 1.)*A_HHe*B*(tinf2 - t), q/(1. - q));
	}
	return L_AGB;}
double Coef::R_TPAGB(double t, double M, double met){
	double A, R_AGB, L;
	L = L_TPAGB(t, M, met);

	if(M >= M_HeF(met)) {
		A = min(b(51, met)*pow(M, -b(52, met)), b(53, met)*pow(M, -b(54, met)));
		R_AGB = A*(pow(L, b(1, met)) + b(2, met)*pow(L, b(55, met)*b(3, met)));
		return R_AGB;}
	if(M <= (M_HeF(met) - 0.2)) {
		A = b(56, met) + b(57, met)*M;
		R_AGB = A*(pow(L, b(1, met)) + b(2, met)*pow(L, b(3, met)));
		return R_AGB;}
	if((M < M_HeF(met)) && (M > (M_HeF(met) - 0.2))) {
		R_AGB = Coef::R_AGBadd(L, M_HeF(met), met) + 
		(Coef::R_AGBadd(L, M_HeF(met), met) - Coef::R_AGBadd(L, M_HeF(met) - 0.2, met))/
		0.2*(M - M_HeF(met) + 0.2);
		return R_AGB;}}

double Coef::Lumt(double t, double M, double met){
	double L = 1;
	if((t >= 0) && (t <= t_MS(M, met))) 								L = L_MS(t, M, met);
	if((t > t_MS(M, met)) && (t <= t_HeI(M, met)))						L = L_HG(t, M, met);
	if(M < M_FGB(met))
		if((t >= t_BGB(M, met)) && (t <= t_HeI(M, met))) 				L = log10(L_GB(t, M, met));
	if((t > t_HeI(M, met)) && (t <= (t_HeI(M, met) + t_He(M, met)))) 	L = log10(L_CHeB(t, M, met));
	if((t > t_HeI(M, met) + t_He(M, met)) && (t <= t_DU(M, met))) 		L = log10(L_EAGB(t, M, met));
	if(M_cBAGB(M, met) < 2.25)
		if((t > t_DU(M, met)) && (t <= t_AGB(M, met)))					L = log10(L_TPAGB(t, M, met));
	return L;}
double Coef::Radt(double t, double M, double met){
	double R = 1;
	if((t >= 0) && (t <= t_MS(M, met))) 								R = R_MS(t, M, met);
	if((t > t_MS(M, met)) && (t <= t_HeI(M, met)))						R = R_HG(t, M, met);
	if(M < M_FGB(met))
		if((t >= t_BGB(M, met)) && (t <= t_HeI(M, met))) 				R = log10(R_GB(t, M, met));
	if((t > t_HeI(M, met)) && (t <= (t_HeI(M, met) + t_He(M, met)))) 	R = log10(R_CHeB(t, M, met));
	if((t > t_HeI(M, met) + t_He(M, met)) && (t <= t_DU(M, met))) 		R = log10(R_EAGB(t, M, met));
	if(M_cBAGB(M, met) < 2.25)
		if((t > t_DU(M, met)) && (t <= t_AGB(M, met)))					R = log10(R_TPAGB(t, M, met));
	return R;}
double Coef::Tt(double t, double M, double met){
	double x;
	x = pow(pow(10, Lumt(t, M, met)), 0.25)*
		pow(pow(10, Radt(t, M, met)), -0.5)*5778;
	return x;}
double Coef::lt(double M, double met){return 1e6*t_DU(M, met);} //1e6*t_DU(M, met);}

// int main(){
// 	cout.precision(10);
// 	cout << fixed;
// 	Coef a;
// 	double M, met, tmax, tmax1, tmin, imin, imax, jmax, imax1, t, T;
// 	system("rm -rf out/*.txt");
// 	double mas[10] = {0.64, 1, 1.6, 2.5, 4, 6.35, 10, 16, 25, 40};
// 	string name[10] = 	{"out/out0.txt", "out/out1.txt", "out/out2.txt", 
// 						"out/out3.txt", "out/out4.txt", "out/out5.txt",
// 						"out/out6.txt", "out/out7.txt", "out/out8.txt",
// 						"out/out9.txt"};				

// 	for(int j = 0; j < 10; j++){
// 		ofstream fout(name[j].c_str());
// 		M = mas[j];
// 		met = 0.02;
// 		imax = 50;
// 		// cout << M << "   " << a.M_cBAGB(M, met) << "    "  << a.t_HeI(M, met) + a.t_He(M, met) 
// 		// 	 << "    " << a.t_DU(M, met) << "    " << a.t_AGB(M, met) << endl;
// 		for(int i = 0; i < imax; i++){
// 			t = a.t_MS(M, met)*i/(imax - 1.);
// 			T = pow(pow(10, a.Lumt(t, M, met)), 0.25)*
// 				pow(pow(10, a.Radt(t, M, met)), -0.5)*5778;
// 			fout << t << "   " << log10(T) << "   " 
// 					<< a.Lumt(t, M, met) << "  " << a.Radt(t, M, met)
// 					<< "   " << i << endl;
// 		}
// 		for(int i = 0; i < imax; i++){
// 			t = a.t_MS(M, met) + (a.t_BGB(M, met) - a.t_MS(M, met))*i/(imax - 1.);
// 			T = pow(pow(10, a.Lumt(t, M, met)), 0.25)*
// 				pow(pow(10, a.Radt(t, M, met)), -0.5)*5778;
// 			fout << t << "   " << log10(T) << "   " 
// 					<< a.Lumt(t, M, met) << "  " << a.Radt(t, M, met)
// 					<< "   " << i << endl;
// 		}
// 		for(int i = 0; i < imax; i++){
// 			t = a.t_BGB(M, met) + (a.t_HeI(M, met) - a.t_BGB(M, met))*i/(imax - 1.);
// 			T = pow(pow(10, a.Lumt(t, M, met)), 0.25)*
// 				pow(pow(10, a.Radt(t, M, met)), -0.5)*5778;
// 			fout << t << "   " << log10(T) << "   " 
// 					<< a.Lumt(t, M, met) << "  " << a.Radt(t, M, met) 
// 					<< "   " << i << endl;
// 		}
// 		for(int i = 0; i < imax; i++){
// 			t = a.t_HeI(M, met) + a.t_He(M, met)*i/(imax - 1.);
// 			T = pow(pow(10, a.Lumt(t, M, met)), 0.25)*
// 				pow(pow(10, a.Radt(t, M, met)), -0.5)*5778;
// 			fout << t << "   " << log10(T) << "   " 
// 					<< a.Lumt(t, M, met) << "  " << a.Radt(t, M, met) 
// 					<< "   " << i << endl;
// 		}
// 		for(int i = 0; i < imax; i++){
// 			t = a.t_HeI(M, met) + a.t_He(M, met) + 
// 				(a.t_DU(M, met) - a.t_HeI(M, met) - a.t_He(M, met))*i/(imax - 1.);
// 			T = pow(pow(10, a.Lumt(t, M, met)), 0.25)*
// 				pow(pow(10, a.Radt(t, M, met)), -0.5)*5778;
// 			fout << t << "   " << log10(T) << "   " 
// 					<< a.Lumt(t, M, met) << "  " << a.Radt(t, M, met) 
// 					<< "   " << i << endl;
// 		}
// 		if(a.M_cBAGB(M, met) < 2.25){
// 			for(int i = 0; i < imax; i++){
// 				t = a.t_DU(M, met) + (a.t_AGB(M, met) - a.t_DU(M, met))*i/(imax - 1.);
// 				T = pow(pow(10, a.Lumt(t, M, met)), 0.25)*
// 					pow(pow(10, a.Radt(t, M, met)), -0.5)*5778;
// 				fout << t << "   " << log10(T) << "   " 
// 						<< a.Lumt(t, M, met) << "  " << a.Radt(t, M, met) 
// 						<< "   " << i << endl;
// 			}
// 		}
// 	}
// 	return 0;}