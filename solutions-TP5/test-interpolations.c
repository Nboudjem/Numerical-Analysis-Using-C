#include <stdio.h>
#include <stdlib.h>
#include "differences_divisees.h"

#define N 5
void main(){
	int n ;
	double X[N+1], Y[N+1], B[(N+1)*(N+1)]; /* X et Y sont les points (X_i,Y_i), suppos�s donn�s */
	/*
		La matrice B contiendra la table des diff. divis�es que l'on calculera :
		B_{ij} == [Y_i, Y_{i+1}, ... , Y_{j-1}, Y_j] 
		Relation de r�currence pour les B_{ij}:
			B(i,j) = ( B(i+1,j) - B(i,j-1) )   /   ( X(j) - X(i) )
			avec, par d�finition, B(i,i) = Y(i)
			et avec i,j=0,...,N
	*/
	double x,h,x_milieu;
	int i,j,l;

	n = Read_X_Y("table5.txt", X, Y , N); /* lire au maximum N+1 points du fichier
						"table5.txt" et mettre les nombres dans X[i], Y[Yi]
						i=0,1,...,n */
	differences_divisees(n, X, Y, B);


	/* On calcule 501 points pour P_n(x) pour x= (X[0]+X[500)/2 + i*1.05*h ;  avec h = ( X[N] - X[0] )/500.
		et i=-250 , ... 250  */
	h = ( X[n] - X[0] )/500. ;
	x_milieu = (X[0]+X[n])/2. ;

	printf("\n# R�sultats obtenus avec P_n(x) pour 501 points\n");
	for (i=-250 ; i <= 250 ; i++){
		x = x_milieu + (double)i*1.05*h ;
		printf("\n%20.16g %20.16g" , x , Polynome_Newton2(x, X, n, B ) );
	}
	printf("\n");

 


	/* ***************************************************************************
		C-idessous, nous G�n�rons une table de 501 valeurs en utilisant
		le polynome de Legendre , et nous mettons nos r�sultats dans le fichier
		"table5_legendre.txt".
	***************************************************************************/
	{	
	FILE *in=fopen("table5_legendre.txt","w");

	/* On calcule 501 points pour P_n(x) pour x= (X[0]+X[500)/2 + i*1.05*h ;  avec h = ( X[N] - X[0] )/500.
		et i=-250 , ... 250  */
	h = ( X[n] - X[0] )/500. ;
	x_milieu = (X[0]+X[n])/2. ;

	fprintf(in,"\n# R�sultats obtenus avec P_n(x) (polyn�me de Legendre) pour 501 points\n");
	for (i=-250 ; i <= 250 ; i++){
		x = x_milieu + (double)i*1.05*h ;
		fprintf(in,"\n%20.16g %20.16g" , x , Polynome_Legendre(x, X, Y, n) );
	}
	fprintf(in,"\n");
	
	fclose(in);
	}	

	return;
}

