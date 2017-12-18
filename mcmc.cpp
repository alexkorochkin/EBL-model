#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>
#include "titles.h"
#include "gsl/gsl_spline.h"
using namespace std;




int main(){
	double rand, accept;
	time_t sec = time(NULL);
	Parameters param(sec);
	
	while(1 > 0){
		param.propose();
		rand = param.doub();
		accept = min(1., param.pprob/param.prob);
		if(rand < accept){
			param.chi2 = param.pchi2;
			param.prob = param.pprob;
			for(int i = 0; i < param.x.size(); i++) param.x[i] = param.propx[i];
			param.print();
		}
		else param.print();
	}
	return 0;
}
