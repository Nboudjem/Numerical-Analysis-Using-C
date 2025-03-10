#include <stdio.h>
#include <math.h>

/* 
Programme qui define une fonction f(A,eps,x0).
Cette fonction est utilis�e dans le programme
principal "main" pour calculer sqrt(1), sqrt(2) ... sqrt(1000.)
avec l'algorithme de Newton.
	precision = eps = 10^(-6)
	valeur initial =  ???????????? Debrouiller-vous.
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
	long i;
	double eps, x,y,A;

   eps = 1.e-6 ;

   for(i=1;i<= 1000 ; i++){
	A = (double)i;
	x = A ;

	y = f(A,eps,x);
	/* Nos affichons les r�sultats pour y*/
	printf("\n\n");
	printf("\n\nRacine carr�e donn�e par sqrt(%f) = %.16f" , A , sqrt(A));
	printf("\nRacine carr�e (y) approch�e de %f = %.16f" , A , y);
	printf("\n\n");
  }
  return 0;
}
