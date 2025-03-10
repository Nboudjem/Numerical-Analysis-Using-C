#include "mes_fonctions.h"
#include <stdio.h>
#include <math.h>

void main(){
	/*
	Ce programme affiche les intervalles [x_i , x_{i+1}] qui
	v�rifient  f_k( x_i ) * f_k( x_{i+1} )  < 0  ,
		pour k=1,2,3 et i=0,1,...,99.
	Les fonctions f1, f2, et f3 sont d�finies dans le fichier "mes_fonctions.h"
	*/

	int i;
	double xi, yi; /* Nous utiliserons yi pour representer x_{i+1} */
	double fxi, fyi; /* valeurs des fonctions f1, f2, f3 au point x_i et  yi = x_{i+1} */



	/* **************************************************************** */
	/* On commence par la fonction f1 */
	/* **************************************************************** */
	xi = -5. ; /* c'est le premier xi */
	fxi = f1(xi) ;
	for (i=0 ; i<100 ; i++){
		yi = - 5. + 0.1 * (double)(i+1) ;
		fyi = f1(yi) ;
		if ( fxi == 0. ){/* Dans ce cas, nous avons une solution que nous affichons*/
			printf("f1 a une solution exacte pour i = %d et x_i = %.16g   avec f(xi) = %.16g\n" , 
				i , xi , fxi );
		}
		if ( fyi == 0. ){/* Dans ce cas, nous avons une solution que nous affichons*/
			printf("f1 a une solution exacte pour i = %d et x_i = %.16g   avec f(xi) = %.16g\n" , 
				i+1 , yi , fyi );
		}
		if ( fxi*fyi < 0. ){/* Dans ce cas, nous avons une solution dans l'intervalle xi,yi */
			printf("f1 a une solution dans l'intervalle [ x[%d] , x[%d] ] = [ %.6g , %.6g ]\n" , 
				i , i+1 , xi , yi );
			printf("\t\t f( %.6g ) = %.6g \t\t f( %.6g ) = %.6g\n" , xi , fxi , yi , fyi );
		}

		xi = yi ;
		fxi = fyi ;		
	}


	/* **************************************************************** */
	/* On passe � la fonction f2 */
	/* **************************************************************** */

	printf("\n\n") ; /* On saute 2 lignes avant d'afficher les r�sultats de f2 */ 

	xi = -5. ; /* c'est le premier xi */
	fxi = f2(xi) ;
	for (i=0 ; i<100 ; i++){
		yi = - 5. + 0.1 * (double)(i+1) ;
		fyi = f2(yi) ;
		if ( fxi == 0. ){/* Dans ce cas, nous avons une solution que nous affichons*/
			printf("f2 a une solution exacte pour i = %d et x_i = %.16g   avec f(xi) = %.16g\n" , 
				i , xi , fxi );
		}
		if ( fyi == 0. ){/* Dans ce cas, nous avons une solution que nous affichons*/
			printf("f2 a une solution exacte pour i = %d et x_i = %.16g   avec f(xi) = %.16g\n" , 
				i+1 , yi , fyi );
		}
		if ( fxi*fyi < 0. ){/* Dans ce cas, nous avons une solution dans l'intervalle xi,yi */
			printf("f2 a une solution dans l'intervalle [ x[%d] , x[%d] ] = [ %.6g , %.6g ]\n" , 
				i , i+1 , xi , yi );
			printf("\t\t f( %.6g ) = %.6g \t\t f( %.6g ) = %.6g\n" , xi , fxi , yi , fyi );
		}

		xi = yi ;
		fxi = fyi ;		
	}



	/* **************************************************************** */
	/* On passe � la fonction f3 */
	/* **************************************************************** */

	printf("\n\n") ; /* On saute 2 lignes avant d'afficher les r�sultats de f3 */ 

	xi = -5. ; /* c'est le premier xi */
	fxi = f3(xi) ;
	for (i=0 ; i<100 ; i++){
		yi = - 5. + 0.1 * (double)(i+1) ;
		fyi = f3(yi) ;
		if ( fxi == 0. ){/* Dans ce cas, nous avons une solution que nous affichons*/
			printf("f3 a une solution exacte pour i = %d et x_i = %.16g   avec f(xi) = %.16g\n" , 
				i , xi , fxi );
		}
		if ( fyi == 0. ){/* Dans ce cas, nous avons une solution que nous affichons*/
			printf("f3 a une solution exacte pour i = %d et x_i = %.16g   avec f(xi) = %.16g\n" , 
				i+1 , yi , fyi );
		}
		if ( fxi*fyi < 0. ){/* Dans ce cas, nous avons une solution dans l'intervalle xi,yi */
			printf("f3 a une solution dans l'intervalle [ x[%d] , x[%d] ] = [ %.6g , %.6g ]\n" , 
				i , i+1 , xi , yi );
			printf("\t\t f( %.6g ) = %.6g \t\t f( %.6g ) = %.6g\n" , xi , fxi , yi , fyi );
		}

		xi = yi ;
		fxi = fyi ;		
	}
return;
}
