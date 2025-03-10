#include <stdlib.h>
#include "algebre.c"

double P( double A[], int n, double lambda){
	/*
	Cette fonction calcule det(A-lambda*I)
	*/
	int i,j;
	double Ap[n*n];
	for (i=1; i<=n ; i++){/* on calcule Ap = A - lambda I */
		for (j=1; j<=n ; j++){
			if (i==j) 	mat(Ap,n,i,i) = mat(A,n,i,i) - lambda ;
			else 		mat(Ap,n,i,j) = mat(A,n,i,j) ;
		}
	}
	return determinant(Ap,n);
}

double f( double x){
	/* $$ f(x) = \sqrt{ 1 - e^{-x^2} } { \log ( 1+x^6 )\over 1+\tan^6( (1+x^2)/5)} $$ */
	double x2,x6,y,ttt;
	x2 = x*x ;
	x6 = x2*x2*x2 ;
	ttt = tan((1.+x2)/5.) ;
	ttt = ttt*ttt*ttt ;
	ttt = ttt*ttt ;
	y = sqrt( 1. - exp(-x2) ) * log( 1. + x6 ) / ( 1. + ttt ) ;
	return y;
}

#define N 100
main(){
/*
	Nous voulons calculer 
		B = \sum_{k=1}^N   A^k / k!
		C = \sum_{k=1}^N    (-1)^k  A^k / k!
*/
	double A[] = { /* matrice 4 x 4 */
		 1,  2,  3, 4,
		 0,  1, -2, 3,
		-2,  1,  5, 2,
		 3, -7,  6, 1
	} ;

	double a[] ={ /* vecteur à 4 composantes */
		 1, 
		 2, 
		 3,
		 4
	};

	double B[16], C[16], U[16], x[4], b[4], D[16], detA, detB, detC, trA, trB, trC;
	int ij,j,k;


	/* Calcul de B  */
	/* Initialisation B = I et U = I */
	matrice_identite(B,4);
	matrice_identite(U,4);

	for(k=1 ; k<=N ; k++){
		produit_matrice_matrice(U,A,D,4);
		produit_matrice_scalaire(D,1./(double)k,U,4);
		somme_matrice_matrice(B,U,B,4);
	}

	/* Calcul de C  */
	/* Initialisation C = I et U = I */
	matrice_identite(C,4);
	matrice_identite(U,4);

	for(k=1 ; k<=N ; k++){
		produit_matrice_matrice(U,A,D,4);
		produit_matrice_scalaire(D,-1./(double)k,U,4);
		somme_matrice_matrice(C,U,C,4);
	}

	printf("B = \n");
	afficher_m(B,4,6);	

	printf("C = \n");
	afficher_m(C,4,6);	

	printf("det B = %.8g\n", determinant(B,4));
	printf("det C = %.8g\n", determinant(C,4));

	trA = mat(A,4,1,1) + mat(A,4,2,2) + mat(A,4,3,3) + mat(A,4,4,4) ;
	printf("e^Tr A = %.8g\n", exp(trA) );
	printf("e^-Tr A = %.8g\n", exp(-trA) );

	/* Question 5 */
    {	int i;
	double lambda;
	double A[] = { /* matrice 4 x 4 */
		 1,   2,   3,  4 ,
		 2,   1,  -2, -7 ,
		 3,  -2,   5 , 6 ,
		 4,  -7,   6 , 1
	} ;
	double xi, xi_plus1, xi_minus1, fxi, fxi_plus1, fxi_minus1, Delta, eps=1.e-12, v_propres[5];

	printf("   i          lambda_i         P(lambda_i)\n");
	for(i=0 ; i<= 100 ; i++){/* Ici on calcule et affiche P(lambda_i) */
		lambda = -20. + 0.4*i ;
		printf(" %d %12.8g  %12.8g\n" , i , lambda , P(A,4,lambda) );
		}
	/* RESULTATS: 4 solutions dans les intervalles suivants:
		[ -8.4 , - 8.0 ]
		[ -1.2 , - 0.8 ]
		[  3.2 ,   3.6 ]
		[ 13.2 ,  13.6 ]
	*/

	/* On commence avec lambda appartenant à  [ -8.4 , - 8.0 ] */
	xi_minus1 	=  - 8.4  ;
	xi 		=  - 8.0  ;
	fxi_minus1	= P(A,4,xi_minus1)  ;
	fxi		= P(A,4,xi)  ;
	do {
		xi_plus1 =  ( xi_minus1 * fxi - xi * fxi_minus1) / (fxi-fxi_minus1) ;
		fxi_plus1	= P(A,4,xi_plus1)  ;

		Delta = xi_plus1 - xi ;

		xi_minus1 = xi ;		fxi_minus1 	= fxi;
		xi	  = xi_plus1 ;		fxi 		= fxi_plus1;
	}
	while( fabs(Delta) > eps );
	v_propres[1] = xi;
	printf("\nValeur propre 1 = %18.12g \t Erreur absolue calculée = %.12g\n" , xi, Delta);

	/* On passe à lambda appartenant à  [ -1.2 , - 0.8 ] */
	xi_minus1 	=  - 1.2  ;
	xi 		=  - 0.8  ;
	fxi_minus1	= P(A,4,xi_minus1)  ;
	fxi		= P(A,4,xi)  ;
	do {
		xi_plus1 =  ( xi_minus1 * fxi - xi * fxi_minus1) / (fxi-fxi_minus1) ;
		fxi_plus1	= P(A,4,xi_plus1)  ;

		Delta = xi_plus1 - xi ;

		xi_minus1 = xi ;		fxi_minus1 	= fxi;
		xi	  = xi_plus1 ;		fxi 		= fxi_plus1;
	}
	while( fabs(Delta) > eps );
	v_propres[2] = xi;
	printf("\nValeur propre 2 = %18.12g \t Erreur absolue calculée = %.12g\n" , xi, Delta);


	/* On passe à lambda appartenant à  [  3.2 ,   3.6 ] */
	xi_minus1 	=  3.2  ;
	xi 		=  3.6  ;
	fxi_minus1	= P(A,4,xi_minus1)  ;
	fxi		= P(A,4,xi)  ;
	do {
		xi_plus1 =  ( xi_minus1 * fxi - xi * fxi_minus1) / (fxi-fxi_minus1) ;
		fxi_plus1	= P(A,4,xi_plus1)  ;

		Delta = xi_plus1 - xi ;

		xi_minus1 = xi ;		fxi_minus1 	= fxi;
		xi	  = xi_plus1 ;		fxi 		= fxi_plus1;
	}
	while( fabs(Delta) > eps );
	v_propres[3] = xi;
	printf("\nValeur propre 3 = %18.12g \t Erreur absolue calculée = %.12g\n" , xi, Delta);


	/* On passe à lambda appartenant à  [ 13.2 ,  13.6 ] */
	xi_minus1 	=  13.2  ;
	xi 		=  13.6  ;
	fxi_minus1	= P(A,4,xi_minus1)  ;
	fxi		= P(A,4,xi)  ;
	do {
		xi_plus1 =  ( xi_minus1 * fxi - xi * fxi_minus1) / (fxi-fxi_minus1) ;
		fxi_plus1	= P(A,4,xi_plus1)  ;

		Delta = xi_plus1 - xi ;

		xi_minus1 = xi ;		fxi_minus1 	= fxi;
		xi	  = xi_plus1 ;		fxi 		= fxi_plus1;
	}
	while( fabs(Delta) > eps );
	v_propres[4] = xi;
	printf("\nValeur propre 4 = %18.12g \t Erreur absolue calculée = %.12g\n" , xi, Delta);

	/* Nous vérifions que det(A) = produit des valeurs propres */
	printf("\nVERIFICATION:\n");
	printf("\tdet(A)                      = %.16g\n" , determinant(A,4));
	printf("\tproduit des valeurs propres = %.16g\n\n" , 
		v_propres[1] * v_propres[2] * v_propres[3] * v_propres[4]
		);


    }



/* *************************************************** */
	/* Probleme 2 */
    {
	int n,i;
	double xi,fxi, xi_p1, fxi_p1, h, sum200=0. , sum300=0. , Delta;
	FILE *outfile=fopen("emd.data","w");



/* ******************************************************************/
	n=200 ;
	h = 1./(double)n ; /* on découpe en "n" l'intervalle [1,2] ;  xi = 1 + h*i */
	for(i=1; i<=n ; i++){
		xi = 1. + h*(double)i ;
		fxi = f(xi) ;
	/*	printf("xi = %.16g\t\t fxi = %.16g\n" , xi , fxi); */
		sum200 += fxi ;
	}		
	sum200 *= h ;




/* ******************************************************************/
	n=300 ;
	h = 1./(double)n ; /* on découpe en "n" l'intervalle [1,2] ;  xi = 1 + h*i */
	for(i=1; i<=n ; i++){
		xi = 1. + h*(double)i ;
		fxi = f(xi) ;
		sum300 += fxi ;
	}		
	sum300 *= h ;

	/* Affichons les résultats: */
	printf("\nL'intégrale est égale à : %.12g +/- %.1g\n" , sum300 , fabs(sum300-sum200) );
	printf("\nL'intégrale est égale à : %.12g +/- %.1g\n" , sum200 , fabs(sum300-sum200) );


/* ******************************************************************/
	/* Voici une autre maniere de calculer l'intégrale et l'erreur commise */
	/* I = sum h * (f_{i+1}+f_{i})/2 +/- E avec erreur E = sum  h* fabs(f_{i+1}-f_{i})/2 */
	/* Donc I = h * sum_{1}^{n-1} f(xi) + 0.5*h* (f(x0) +f(xn)) +/- Delta */
	n=1000 ;
	h = 1./(double)n ; /* on découpe en "n" l'intervalle [1,2] ;  xi = 1 + h*i */
	Delta =  0. ;
	sum300 = 0. ;

	xi = 1. ;
	fxi = f(xi) ;
	sum300 += fxi/2. ;

	for(i=1; i<n ; i++){
		xi_p1 = 1. + h*(double)i ;
		fxi_p1 = f(xi_p1) ;
		sum300 += fxi_p1 ;
		if (fxi > fxi_p1) 	Delta += fxi - fxi_p1 ;
		else		 	Delta -= fxi - fxi_p1 ;
		fxi = fxi_p1 ;
	}		

	xi_p1 = 1. + h*(double)i ;
	fxi_p1 = f(xi_p1) ;
	sum300 += fxi_p1/2. ;
	if (fxi > fxi_p1) 	Delta += fxi - fxi_p1 ;
	else		 	Delta -= fxi - fxi_p1 ;

	sum300 *= h ;
	Delta *= 0.5*h ;
	/* Affichons les résultats: */
	printf("\nL'intégrale (trapèzes) est égale à : %.12g +/- %.1g\n" , sum300 , Delta );



/* ******************************************************************/
	/*On veut faire un graphe de f(x) */
/*
	n=1000 ;
	h = 1./(double)n ; 
	for(i=1; i<=n ; i++){
		xi = 1. + h*(double)i ;
		fxi = f(xi) ;
		fprintf(outfile,"%8g\t%12.6g\n",xi,fxi);
	}		
*/
/* ******************************************************************/





	n=1000000 ;
	sum300 = 0. ;
	h = 1./(double)n ; /* on découpe en "n" l'intervalle [1,2] ;  xi = 1 + h*i */
	for(i=1; i<=n ; i++){
		xi = 1. + h*(double)i ;
		fxi = f(xi) ;
		sum300 += fxi ;
	}		
	sum300 *= h ;

	/* Affichons les résultats: */
	printf("\nL'intégrale est égale à : %.12g +/- %.1g\n" , sum300 , fabs(sum300-sum200) );


	n=200000 ;
	h = 1./(double)n ; /* on découpe en "n" l'intervalle [1,2] ;  xi = 1 + h*i */
	for(i=1; i<n ; i++){
		xi = 1. + h*(double)i ;
		fxi = f(xi) ;
	/*	printf("xi = %.16g\t\t fxi = %.16g\n" , xi , fxi); */
		sum200 += fxi ;
	}		
	sum200 += (f(1.) + f(2.))/2. ;
	sum200 *= h ;

	n=300000 ;
	h = 1./(double)n ; /* on découpe en "n" l'intervalle [1,2] ;  xi = 1 + h*i */
	for(i=1; i<n ; i++){
		xi = 1. + h*(double)i ;
		fxi = f(xi) ;
		sum300 += fxi ;
	}		
	sum300 += (f(1.) + f(2.))/2. ;
	sum300 *= h ;

	printf("\nL'intégrale''''''' est égale à : %.12g +/- %.1g\n" , sum300 , fabs(sum300-sum200) );


    }
/* *************************************************** */


return;
}

