#include <stdio.h>
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
	/* Nos affichons les résultats pour x*/
	printf("\nRacine carrée donnée par sqrt(%f) = %.16f" , A , sqrt(A));
	printf("\nNombre d'iterations n = %ld" , n );
	printf("\nRacine carrée (x) approchée de %f = %.16f" , A , x);
	printf("\nRacine carrée (y) approchée de %f = %.16f" , A , y);
	printf("\n\nErreur absolue  \t= %1.e" , fabs(y-x) );
	printf(  "\nErreur relative \t= %1.e\n\n" , fabs(y-x)/y );
	return 0;
}
