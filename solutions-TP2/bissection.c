#include <stdio.h>
#include <math.h>

/*
	Ce programme utilise la bissection pour trouver les z�ros
	des fonctions f1, f2, f3. Ces fonctions sont d�finies dans le
	fichier "mes_fonctions.h". Nous incluons ce fichier.
	Les intervalles [a,b] contenant les z�ros peuvent etre trouv�s dans
	le fichier "zeros.output" (ce fichier contient le resultat de la question 3).
*/

#include "mes_fonctions.h"

int main(){
	long n ; 
	int exact_solution ;
	double a,b,c, fa,fb,fc, Delta, eps=1.e-15 ;




	/* On commence avec la fonction f1. Une solution appartenant � [ -1.6 , -1.5 ] */
	a = -1.6 ;		b = -1.5 ;
	fa = f1(a) ;		fb = f1(b) ;
	exact_solution = 0 ;
	n = 0 ;
	do {
		n = n + 1 ; /* n=nombre d'it�rations */
		c = 0.5*(a+b) ;	/* Pour la position fausse on utilise un autre "c" : 
				c = (a*fb - b*fa)/(fb-fa) ; */

		fc = f1(c) ;
		if (fc == 0.) {
			exact_solution = 1 ;
			printf("\n\nPour f1, nous avons une solution exacte x = %.16g\n", c);
			break; /* break veut dire: nous sortons de la boucle do */
		}
		else if (fa*fc < 0) { b = c ; }
		else		    { a = c ; }

		Delta = 0.5*(b-a) ; /* erreur absolue */
	}
	while( fabs(Delta) > eps ) ;

	if (exact_solution == 0 ) { /* 0 veut dire que la solution n'est pas exacte */
		printf("\n\nPour f1, nous avons une solution approch�e x = %.16g", c ) ;
		printf("\n\tv�rifiant f1(x) = %.16g\n", f1(c) );
	}







	/* Pour la fonction f2. Une solution triviale x =0 et 2 autres solutions appartenant 
	   aux intervalles [ -0.9 , -0.8 ]  et [ 2.2 , 2.3 ] */

	/* Cherchons d'abord pour f2 la solution dans [ -0.9 , -0.8 ] */
	a = -0.9 ;		b = -0.8 ;
	fa = f2(a) ;		fb = f2(b) ;
	exact_solution = 0 ;
	n = 0 ;
	do {
		n = n + 1 ; /* n=nombre d'it�rations */
		c = 0.5*(a+b) ;	/* Pour la position fausse on utilise un autre "c" : 
				c = (a*fb - b*fa)/(fb-fa) ; */

		fc = f2(c) ;
		if (fc == 0.) {
			exact_solution = 1 ;
			printf("\n\nPour f2, nous avons une solution exacte x = %.16g\n", c);
			break; /* break veut dire: nous sortons de la boucle do */
		}
		else if (fa*fc < 0) { b = c ; }
		else		    { a = c ; }

		Delta = 0.5*(b-a) ; /* erreur absolue */
	}
	while( fabs(Delta) > eps ) ;

	if (exact_solution == 0 ) { /* 0 veut dire que la solution n'est pas exacte */
		printf("\n\nPour f2, nous avons une solution approch�e x = %.16g", c ) ;
		printf("\n\tv�rifiant f2(x) = %.16g\n", f2(c) );
	}


	/* Cherchons la deuxi�me solution pour f2  [ 2.2 , 2.3 ] */
	a = 2.2 ;		b = 2.3 ;
	fa = f2(a) ;		fb = f2(b) ;
	exact_solution = 0 ;
	n = 0 ;
	do {
		n = n + 1 ; /* n=nombre d'it�rations */
		c = 0.5*(a+b) ;	/* Pour la position fausse on utilise un autre "c" : 
				c = (a*fb - b*fa)/(fb-fa) ; */

		fc = f2(c) ;
		if (fc == 0.) {
			exact_solution = 1 ;
			printf("\n\nPour f2, nous avons une solution exacte x = %.16g\n", c);
			break; /* break veut dire: nous sortons de la boucle do */
		}
		else if (fa*fc < 0) { b = c ; }
		else		    { a = c ; }

		Delta = 0.5*(b-a) ; /* erreur absolue */
	}
	while( fabs(Delta) > eps ) ;

	if (exact_solution == 0 ) { /* 0 veut dire que la solution n'est pas exacte */
		printf("\n\nPour f2, nous avons une solution approch�e x = %.16g", c ) ;
		printf("\n\tv�rifiant f2(x) = %.16g\n", f2(c) );
	}









	/* Pour la fonction f3. Une solution triviale x =0 et 2 autres solutions appartenant 
	   aux intervalles  [ -0.4 , -0.3 ] et [ 0.3 , 0.4 ]. 
	   Comme f3 est impaire, donc si x est solution, alors -x est aussi solution
	   On en calcule une seule appartenant � [ 0.3 , 0.4 ] */

	a = 0.3 ;		b = 0.4 ;
	fa = f3(a) ;		fb = f3(b) ;
	exact_solution = 0 ;
	n = 0 ;
	do {
		n = n + 1 ; /* n=nombre d'it�rations */
		c = 0.5*(a+b) ;	/* Pour la position fausse on utilise un autre "c" : 
				c = (a*fb - b*fa)/(fb-fa) ; */

		fc = f3(c) ;
		if (fc == 0.) {
			exact_solution = 1 ;
			printf("\n\nPour f3, nous avons une solution exacte x = %.16g avec f(x)=%.16g\n", c, f3(c));
			printf("\t -x est aussi solution = %.16g \n", -c ) ;
			break; /* break veut dire: nous sortons de la boucle do */
		}
		else if (fa*fc < 0) { b = c ; }
		else		    { a = c ; }

		Delta = 0.5*(b-a) ; /* erreur absolue */
	}
	while( fabs(Delta) > eps ) ;

	if (exact_solution == 0 ) { /* 0 veut dire que la solution n'est pas exacte */
		printf("\n\nPour f3, nous avons une solution approch�e x = %.16g", c ) ;
		printf("\n\tv�rifiant f3(x) = %.16g", f2(c) );
		printf("\n\tL'autre solution est -x = %.16g\n", -c );
	}


return 0;
}
