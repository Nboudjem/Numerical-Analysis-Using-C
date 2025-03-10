#include <stdio.h>
//extern int printf (__const char *__restrict __format, ...);

#include <math.h>

int main(){
	long n;
	double eps=0.e-6, x,y,A;
	/*double eps=0., x,y,A;*/

		n=0;
		A=17.3;
	/*	x=4. ;*/
		x=4000.16 ;

	n=1;
	y = 0.5*(x+A/x) ;

	while( fabs(x-y) > eps ){
		x=y;
		n=n+1;
		y = 0.5*(x+A/x) ;
	}
	/* Nos affichons les r�sultats pour x*/
	printf("\nRacine carr�e donn�e par sqrt(%f) = %.16f" , A , sqrt(A));
	printf("\nNombre d'iterations n = %ld" , n );
	printf("\nRacine carr�e (x) approch�e de %f = %.16f" , A , x);
	printf("\nRacine carr�e (y) approch�e de %f = %.16f" , A , y);
	printf("\n\nErreur absolue  \t= %1.e" , fabs(y-x) );
	printf(  "\nErreur relative \t= %1.e\n\n" , fabs(y-x)/y );

	return 0;
}
