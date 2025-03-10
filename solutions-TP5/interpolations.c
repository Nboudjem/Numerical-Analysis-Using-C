/* 
	Ce fichier contient des fonctions qui permettent de calculer
	les polynomes d'interpolations P_n(x) qui passent
	par un ensemble de (n+1) points donn�s { (X_i, Y_i)   ,  i=0,...,n 
		1. Calcul de la table des diff�rences divis�es
		2. Obtention de la premi�re formule d'interpolation de Newton
		3. Obtention de la deuxieme formule d'interpolation de Newton
		4. M�thode directe pour obtenir P_n(x) en calculant l'inverse d'une matrice
			(de type Vondermonde)
		5. Comparaison avec les r�sultats obtenus avec le polynome de Lagrange.
*/

#include <stdio.h>
#define N 5
#define mat(B,n,i,j) 	B[ (i)*(n) + (j) ] /* element B_ij , B matrice (n+1)x(n+1) , i,j = 0,...,n */

double Polynome_Newton1(double x, double X[], int n, double B[] ){
	/*
	Cette fonction calcule P_n(x), 1�re formule d'interpolation polynomiale de Newton.
	La matrice B, (n+1)x(n+1), est la table des diff�rences divis�es, suppos�ee donn�.

	RAPPEL: la table des diff. divis�es est calcul�e comme ceci :
		B(i,j) == [Y_i, Y_{i+1}, ... , Y_{j-1}, Y_j] 
		Relation de r�currence pour les B(i,j):
			B(i,j) = ( B(i+1,j) - B(i,j-1) )   /   ( X(j) - X(i) )  avec i<= j
			avec, par d�finition, B(i,i) = Y(i)
			et avec i<=j  qui peuvent prendre les valeurs 0,...,N
	
	P_n(x) = a0     + a1 * (x - X[0])
			+ a2 * (x - X[0]) * (x - X[1])
			+ ...
			+ an * (x - X[0]) * (x - X[1]) * ... * (x-X[n-1])

	Les coefficients ai sont les differences divis�es, �gales �: a_i = B_{0,i} pour la
	1�re formule d'interpolation de Newton.
	*/
	double y=1. , sum;
	int i;
		
	sum = mat(B,n,0,0) ; /* + a0 */
	for (i=1 ; i<= n ; i++){
		y *= x-X[i-1] ;
		sum += mat(B,n,0,i) * y ;
	}	

	return sum;
}



void main(){

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

	/* exemple de points pour N=5 */
	X[0] =	1.1	;	Y[0]=	12.1	;
	X[1] =	1.15	;	Y[1]=	1.5	;
	X[2] =	1.2	;	Y[2]=	-2.5	;
	X[3] =	1.3	;	Y[3]=	-16.35	;
	X[4] =	1.45	;	Y[4]=	0.13	;
	X[5] =	1.65	;	Y[5]=	18.34	;


	/* ***************************************************************************
	Calcul de la table des diff�rences divis�es :
	1. On calcule d'abord B(i,i+1) pour i=0,...,n-1
	2. On calcule ensuite B(i,i+2) pour i=0,...,n-2
	3. ........
	4. On calcule finalement B(i,n). pour i=0,...,0 (une seule valeur)
   	***************************************************************************/

	for(i=0 ; i<= N ; i++)	mat(B,N,i,i) = Y[i] ; /* initialisation de B(i,i) = Y_i */

	for(l=1 ; l<= N ; l++){
		for(i=0 ; i<= N-l ; i++){
			j = i+l ;
			mat(B,N,i,j) = ( mat(B,N,i+1,j) - mat(B,N,i,j-1) )/ ( X[j] - X[i] ) ;
		}
	}

	/* Affichons la table B_ij (uniquement les �lements calcul�s: i<= j ) */
	printf("\n\nB = \n");
	for(i=0 ; i<= N ; i++){
	printf("\n\t");
		for(j=0 ; j<= N ; j++) {
			if (i>j)	printf("%16.0f", 0. );
			else 		printf("%16.8f",mat(B,N,i,j) );
		}
	}
	printf("\n\n\n\n");


	/* On calcule 501 points pour P_n(x) pour x= (X[0]+X[500)/2 + i*1.05*h ;  avec h = ( X[N] - X[0] )/500.
		et i=-250 , ... 250  */
	h = ( X[N] - X[0] )/500. ;
	x_milieu = (X[0]+X[N])/2. ;

	printf("\n R�sultats obtenus avec P_n(x) pour 501 points\n");
	for (i=-250 ; i <= 250 ; i++){
		x = x_milieu + (double)i*1.05*h ;
		printf("\n%20.16g %20.16g" , x , Polynome_Newton1(x, X, N, B ) );
	}
	printf("\n");

 
	/* ***************************************************************************/
		/* R�ponse � la question 1 */
	/* ***************************************************************************/

	return;
}

