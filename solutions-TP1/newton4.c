#include <stdio.h>
#include <math.h>

/* 
Programme qui define une fonction f(A,eps,x0).
Cette fonction est utilis�e dans le programme
principal "main" pour calculer sqrt(A)
avec l'algorithme de Newton.
	precision = eps
	valeur initial =  x0
*/

double f(double A, double eps, double x){
	double y;

	y = 0.5*(x+A/x) ;
	while( fabs(x-y) > eps ){
		x=y;
		y = 0.5*(x+A/x) ;
	}
return y;
}


int main(){
	long n;
	double eps, x,y,A;

	printf("Calculer la racine carr�e de A = ");
		scanf("%lf",&A);
	printf("avec une pr�cision de  eps = ");
		scanf("%lf",&eps);
	printf("et une valeur initiale = ");
	scanf("%lf",&x);

	y = f(A,eps,x);
	/* Nos affichons les r�sultats pour y*/
	printf("\n\n");
	printf("\n\nRacine carr�e donn�e par sqrt(%f) = %.16f" , A , sqrt(A));
	printf("\nRacine carr�e (y) approch�e de %f = %.16f" , A , y);
	printf("\n\n");

	return 0;
}
