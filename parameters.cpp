#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "titles.h"
#include "spline.h"
using namespace std;


Ran::Ran(unsigned long long int j){
	v = 4101842887655102017LL;
	w = 1;
	u = j^v;
	int64();
	v = u;
	int64();
	w = v;
	int64();}
unsigned long long int Ran::int64(){
	u = u* 2862933555777941757LL + 7046029254386353087LL;
	v ^= v >> 17;
	v ^= v << 31;
	v ^= v >> 8;
	w = 4294957665U*(w & 0xffffffff) + (w >> 32);
	unsigned long long int x = u ^ (u << 21);
	x ^= x >> 35;
	x ^= x << 4;
	return (x + v) ^ w;}
double Ran::doub(){return 5.42101086242752217E-20*int64();}
double Ran::normaldev(double mu, double sig){
	double a, b, x, y, q;
	do {
		a = doub();
		b = 1.7156* (doub() - 0.5);
		x = a - 0.449871;
		y = fabs(b) + 0.386595;
		q = pow(x, 2) + y*(0.19600*y - 0.25472*x);
	} while (q > 0.27597 && (q > 0.27846 || pow(b, 2) > -4.*log(a)*pow(a, 2)));
	return mu + sig*b/a;}

Parameters::Parameters(unsigned long long int ranseed): Ran(ranseed){
	double q;
	string a;

	ifstream fin("in/param.txt");
	x.clear();
	propx.clear();
	xmax.clear();
	xmin.clear();
	while(fin >> a){ 
		fin >> q;
		x.push_back(q);
		propx.push_back(q);
		fin >> q;
		xmin.push_back(q);
		fin >> q;
		xmax.push_back(q);
		fin >> q;
		dx.push_back(q);
	}
	fin.close();
	fin.open("in/lower_limits.txt");
	while(fin >> q){
		lower_limits_x.push_back(q);
		fin >> q;
		lower_limits_y.push_back(q);
		fin >> q;
		lower_limits_upsigma.push_back(q);
		fin >> q;
		lower_limits_downsigma.push_back(q);
	}
	fin.close();
	fin.open("in/upper_limits.txt");
	while(fin >> q){
		upper_limits_x.push_back(q);
		fin >> q;
		upper_limits_y.push_back(q);
		fin >> q;
		upper_limits_upsigma.push_back(q);
		fin >> q;
		upper_limits_downsigma.push_back(q);		
	}
	fin.close();
	fin.open("in/direct.txt");
	while(fin >> q){
		direct_x.push_back(q);
		fin >> q;
		direct_y.push_back(q);
		fin >> q;
		direct_upsigma.push_back(q);		
		fin >> q;
		direct_downsigma.push_back(q);		
	}
	fin.close();


	spec(x);
	chsq(chi2, prob);
	print();}
void Parameters::propose(){
	int point, flag = 0;
	double prop_point;
	for(int i = 0; i < x.size(); i++) propx[i] = x[i];
	while(flag != 1){
		point = int(doub()*(x.size() + 1)); 
		if (dx[point] != 0.){
			prop_point = x[point] + normaldev(0, dx[point]);
			if ((prop_point > xmin[point]) && (prop_point < xmax[point])){
				flag = 1;
				propx[point] = prop_point;
			}
		}
	}
	spec(propx);
	chsq(pchi2, pprob);}
void Parameters::print(){
	ofstream fout("out/mcmc.txt", iostream::app);
	fout << 1 << "  " << prob << "  ";
	for(int i = 0; i < x.size(); i++) fout << x[i] << "  ";
	fout << endl;
	fout.close();}
void Parameters::chsq(double &chi2, double &prob){
	vector<double> xx, yy;
	tk::spline spln;
	double temp;
	chi2 = 0.0;

	ifstream fin("out/spec.txt");
	fin >> temp;
	while(!fin.eof()){
		xx.push_back(temp);
		fin >> temp;
		yy.push_back(temp);
		fin >> temp;
	}		
	fin.close();
	spln.set_points(xx, yy);

	for(int j = 0; j < lower_limits_x.size(); j++)
		if(spln(lower_limits_x[j]) < lower_limits_y[j]){
			temp = spln(lower_limits_x[j]) - lower_limits_y[j];
			chi2 += temp*temp/(2*lower_limits_downsigma[j]);
		}
	for(int j = 0; j < upper_limits_x.size(); j++)
		if(spln(upper_limits_x[j]) > upper_limits_y[j]){
			temp = spln(upper_limits_x[j]) - upper_limits_y[j];
			chi2 += temp*temp/(2*upper_limits_upsigma[j]);
		}
	for(int j = 0; j < direct_x.size(); j++){
		temp = spln(direct_x[j]) - direct_y[j];
		chi2 += temp*temp/(2*direct_upsigma[j]);
	}
	prob = exp(-chi2);}
