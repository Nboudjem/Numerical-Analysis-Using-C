#include <stdio.h>
#include <math.h>



double f1(double x){
	/* f1 = x^5 + 3 x^2 + 1 */
	return x*x*(x*x*x+3.) +1. ;
}


double f1_prime(double x){
	/* f1_prime est la dérivée de f1.
		f1_prime = 4 x^4 + 6 x  */
	return x*(4.*x*x*x+6.)  ;
}



double f2(double x){
	/* f2 = x log(1+x^4) - e^x sin(x)  */
	return x*log(x*x*x*x+1.) - exp(x)*sin(x) ;
}

double f2_prime(double x){
	/* f2_prime est la dérivée de f2.
	   f2_prime =  log(1+x^4) + x (4 x^3) / (1+x^4) - e^x sin(x) - e^x cos(x)  */
	double x4 ;
	x4 = x*x*x*x ;
	return log(1.+x4) + 4. * x4 / (1.+x4) - exp(x)*( cos(x) + sin(x) ) ;
}



double f3(double x){
	/* f3 = \arcsin(x/10) - 3 x \sqrt(1+x^2) + \pi \arctan(x)  */
	double pi;
	pi = 4. * atan(1.) ;
	return asin(x/10.) - 3. * x * sqrt(1.+x*x) + pi * atan(x) ;
}

double f3_prime(double x){
	/* f3_prime est la dérivée de f3.
	   f3_prime = 1 / 10 / \sqrt( 1 - (x/10)^2 ) - 3 \sqrt(1+x^2) - 3 x^2 / \sqrt(1+x^2) 
		+ \pi / ( 1 + x^2)  */
	double pi, x2;
	pi = 4. * atan(1.) ;
	x2 = x*x ;
	return 	0.1/ sqrt( 1. - 0.01*x2 ) 
		- 3. * sqrt(1.+x2) - 3. * x2 / sqrt(1.+x2) 
		+ pi / (1.+x2) ;
}
