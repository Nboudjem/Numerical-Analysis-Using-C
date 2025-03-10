/*
	This code uses netwton method to find f(x) = 0
	where the used initial value x_0 = 2.
*/

#include <stdio.h>
#include <math.h>


double f( double x){
	return 3.*x + exp(x) ;
}
double f_prime( double x){
	return 3. + exp(x) ;
}



void main(){
	long n=0;
	double x,y, Delta, Delta_relative , eps=1.e-15;

	x = 2. ; /* vinitial value x_0 */

	printf("\n n              x                   y                     Delta               Delta_r\n" );

	do {/* use of loop do {...} while(...) */
		n = n + 1 ;
		y = x - f(x)/f_prime(x) ;
		Delta = y-x ;
		Delta_relative = Delta/y  ;
		printf("\n%2ld  %19.16f  %19.16f  %19.1e  %19.1e " , n , x , y , Delta, Delta_relative);

		x = y ;
	}
	while( fabs(Delta_relative) > eps );

	printf("\n\n" );

	return;
}
