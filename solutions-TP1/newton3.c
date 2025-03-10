#include <stdio.h>
#include <math.h>

/* 
Programme qui lie les nombres reels A, eps et x0
et qui calcule sqrt(A) avec une precision de eps,
et une valeure initiiale x0.
*/

int main(){
	long n;
	double eps, x,y,A;

	printf("Calculer la racine carr�e de A = ");
		scanf("%lf",&A);
	printf("avec une pr�cision de  eps = ");
		scanf("%lf",&eps);
	printf("et une valeur initiale = ");
	scanf("%lf",&x);

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
