#include "mes_fonctions.h"
#include <stdio.h>
#include <math.h>

void main(){
	/*
	Ce programme calcule f_k(x_i) pour k=1,2,3 et i=0,1,...,100.
	Les fonctions f1, f2, et f3 sont d�finies dans le fichier "mes_fonctions.h"
	Les r�sultats sont affich�s comme un tableau � 4 colonnes, chaque ligne
	ressemblant � :
	x_i	f1(x_i)		f2(x_i)		f3(x_i)
	*/

	int i;
	double xi;

	for (i=0 ; i<=100 ; i++){
		xi = - 5. + 0.1 * (double)i ;
		printf("%11.6g %11.6g %11.6g %11.6g\n" , xi , f1(xi) , f2(xi) , f3(xi) );
	}
return;
}
