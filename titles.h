#include <vector>
#include <cmath>
using namespace std;


//track
struct Coef{
	Coef();
	double p[98][5];
	double q[58][5];
	double func_a(int, double);
	double func_b(int, double);
	double a(int, double);
	double b(int, double);

	double t_hook(double, double);
	double t_MS(double, double);
	double t_BGB(double, double);
	double t_HeI(double, double);
	double t_HeMS(double);
	double t_He(double, double);
	double tau_bl(double, double);
	double t_DU(double, double);
	double t_AGB(double, double);

	double M_cBGB(double, double met);
	double M_cHeI(double, double met);
	double M_cBAGB(double, double met);
	double M_cGB(double, double M, double met);
	double M_c(double, double M, double met);

	double M_HeF(double);
	double M_FGB(double);
	double L_BGB(double, double);
	
	double L_ZAMS(double, double);
	double R_ZAMS(double, double);
	
	double L_TMS(double, double);
	double R_TMS(double, double);
	
	double L_EHG(double, double);
	double R_EHG(double, double);
	
	double L_ZHe(double);
	double R_ZHe(double);

	double L_minHe(double, double);
	double R_mHe(double, double);

	double L_HeI(double, double);	
	double R_HeI(double, double);

	double L_ZAHB(double, double);
	double R_ZAHB(double, double);

	double L_BAGB(double, double);
	double R_BAGB(double, double);

	double R_GBadd(double, double, double);
	double R_AGBadd(double, double, double);

	double L_MS(double, double, double);
	double R_MS(double, double, double);
	double L_HG(double, double, double);
	double R_HG(double, double, double);
	double L_GB(double, double, double);
	double R_GB(double, double, double);
	double L_CHeB(double, double, double);
	double R_CHeB(double, double, double);
	double L_EAGB(double, double, double);
	double R_EAGB(double, double, double);
	double L_TPAGB(double, double, double);
	double R_TPAGB(double, double, double);

	double Lumt(double, double, double);
	double Radt(double, double, double);
	double Tt(double, double, double);
	double lt(double, double);};



//spec
struct Supp;
double zt(const double);
double tz(const double);
double rz(const double);
double zr(const double);
double imf(const double, const Supp *);
double sfr(const double, const Supp *);

void print(const double *, const Supp *);
void starspec(const double, const double, const double, double *, const Supp *);
void basicspec(const double, const double, double *, const Supp *);
void basicspecdust(const double, const double, double *, const Supp *);
double dusttemp(const int, const int, const double *, const Supp *);
void spec(const vector<double> &);





//proposal
struct Ran{
	unsigned long long int u, v, w;
	Ran(unsigned long long int);
	unsigned long long int int64();
	double doub();
	double normaldev(double, double);};
struct Parameters: Ran{
	double chi2, prob, pchi2, pprob;
	vector<double> x, dx, xmin, xmax, propx;
	vector<double> lower_limits_x, lower_limits_y, lower_limits_upsigma, lower_limits_downsigma,
				   upper_limits_x, upper_limits_y, upper_limits_upsigma, upper_limits_downsigma,
				   direct_x, direct_y, direct_upsigma, direct_downsigma;
	Parameters(unsigned long long int);
	void propose();
	void chsq(double &, double &);
	void print();};





//mcmc
void mcmcstep(Parameters &, double &);
