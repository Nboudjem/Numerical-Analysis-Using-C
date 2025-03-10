/*
	Ce programme utilise la methode de bissection 
	(ou de la position fausse) pour
	trouver une racine de f(x)=0 dans l'intervalle
	[a,b], verifiant  f(a)f(b)  < 0
*/

#include <stdio.h>
#include <math.h>

double f( double x){ 
	return 3.*x + exp(x); 
}

void main(){
	long n=0;
	double a,b,c, fa,fb,fc, Delta, Delta_relative , eps=1.e-15;

	a = -1. ;
	b = 0. ;

	fa = f(a) ;
	fb = f(b) ;

	if ( fa*fb > 0. ) {
		printf("\n\n f(a) f(b) > 0. Choisissez un autre intervalle verifiant f(a)f(b)<0 \n\n");
		return;
	}

	Delta = 0.5*(b-a) ;
	c = 0.5*(a+b) ;	/* Pour la position fausse on utilise un autre "c" : 
				c = (a*fb - b*fa)/(fb-fa) ;
			Autrement le programme est partout identique */

	printf("\n n              a                   b                     c                     Delta \n" );

	while( fabs(Delta/c) > eps ){
		n = n + 1 ;
		/* nous affichons les resultats pour chaque iteration */
		printf("\n%3ld  %19.15f  %19.15f  %19.15f  %19.1e " , n , a , b , c , Delta);
	
		fc = f(c) ;
		if (fc == 0.) {
			printf("\n-------------------------------------");
			printf("\nWe have an exact solution x = %.16g", c);
			printf("\n-------------------------------------");
			break; /* break veut dire: nous sortons de la boucle while */
		}
		else if (fa*fc < 0){
			b = c ;
		}
		else{
			a = c ;
		}

		Delta = 0.5*(b-a) ;
		c = 0.5*(a+b) ;	/* Pour la position fausse on utilise un autre "c" : 
					c = (a*fb - b*fa)/(fb-fa) ;
				Autrement le programme est partout identique */
	}
	printf("\n\n" );

	return;
}
