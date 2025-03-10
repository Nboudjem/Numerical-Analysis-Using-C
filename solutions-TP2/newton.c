/*
	Ce programme utilise la methode de newton pour
	trouver les racines de f2(x)=0 
	en commen�ant par une valeur initiale x_0 obtenue
	de la question 3.
	Pour la fonction f2. Une solution triviale x =0 et 2 autres solutions appartenant 
	aux intervalles [ -0.9 , -0.8 ]  et [ 2.2 , 2.3 ] 
*/

#include <stdio.h>
#include <math.h>

#include "mes_fonctions.h"


void main(){
	long n;
	double x,y, Delta, Delta_relative , eps=1.e-15; /* ici eps=erreur relative */


	/* Cherchons d'abord pour f2 la solution dans [ -0.9 , -0.8 ] */
	x = -0.8 ; /* valeur initiale x_0 */
	n = 0 ;
	do {
		n = n + 1 ;
		y = x - f2(x)/f2_prime(x) ; /* x = x_{i}  et  y = x_{i+1} */
		Delta = fabs(y-x) ;	/* erreur absolue */
		Delta_relative = fabs (Delta/y)  ; /* erreur relative */
		x = y ;
	}
	while( fabs(Delta_relative) > eps );

	printf("\n\nPour f2, nous avons une solution approch�e x = %.16g", x ) ;
	printf("\n\tv�rifiant f2(x) = %.16g\n", f2(x) );
	printf("\n\tL'erreur relative calcul�e est : %.1g\n\n", Delta_relative );







	/* Cherchons maintenant la solution dans [ 2.2 , 2.3 ]  */
	x = 2.3 ; /* valeur initiale x_0 */
	n = 0 ;
	do {
		n = n + 1 ;
		y = x - f2(x)/f2_prime(x) ; /* x = x_{i}  et  y = x_{i+1} */
		Delta = fabs(y-x) ;	/* erreur absolue */
		Delta_relative = fabs (Delta/y)  ; /* erreur relative */
		x = y ;
	}
	while( fabs(Delta_relative) > eps );

	printf("\n\nPour f2, nous avons une solution approch�e x = %.16g", x ) ;
	printf("\n\tv�rifiant f2(x) = %.16g\n", f2(x) );
	printf("\n\tL'erreur relative calcul�e est : %.1g\n\n", Delta_relative );

	return;
}
