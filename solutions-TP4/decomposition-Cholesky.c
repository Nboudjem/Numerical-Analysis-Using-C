/* ***************************************************************************/
/* ***************************************************************************/

#include "algebre.c"

void decomposition_Cholesky( double A[], int n, double L[])
{
/* 
   Cette fonction décompose une matrice $A$ carrée, réelle, symétrique,
   définie positive, en un produit $L.L^t$, avec $L$ triangulaire inférieure.
   Les équations sont:
	${\rm\ pour\ } i=1,\dots,n :$
	1.	$L_{ij} = 0 {\rm\ si\ } j>i $
	2.	$L_{ij} = { 1 \over L_{jj} } \ ( A_{ij} - \sum_{k=1}^{j-1} L_{ik} L_{jk} ), {\rm pour\ }j=1,\dots,i-1$	
	3.	$L_{ii} = \sqrt{ A_{ii} - \sum_{k=1}^{i-1} L_{ik}^2 }$	
*/

	int i,j,k;
	for(i=1; i<=n ; i++)
		for (j=i+1; j<=n ; j++)
			mat(L,n,i,j) = 0. ;/* des zéros dans la partie supérieure de L */
	for(i=1; i<=n ; i++){
		for (j=1; j<i ; j++){
			/* Calcul de L_ij */
			mat(L,n,i,j) = mat(A,n,i,j) ;
			for (k=1; k<j ; k++)	mat(L,n,i,j) -= mat(L,n,i,k)*mat(L,n,j,k) ;
			mat(L,n,i,j) /=  mat(L,n,j,j) ;
		}

		/* Calcul de L_ii */
		mat(L,n,i,i) = mat(A,n,i,i) ;
		for (k=1; k<i ; k++)	mat(L,n,i,i) -= mat(L,n,i,k)*mat(L,n,i,k) ;
		mat(L,n,i,i) = sqrt( mat(L,n,i,i) ) ;
	}
return;
}


main(){
	double A[] = { /* matrice 3 x 3 */
		 8,  -2,  -2 , 
		-2,   4,  -2 , 
		-2,  -2,  13
	} ;

	double C[] = { /* matrice 8 x 8 */
		 2,   3,  -1 ,    4,   1,  -2,   5,  32, 
		 1,  30,   5 ,    3,  -2,   3,   0,  10, 
		 6,  -3,  40 ,   -1,   1,   2,  -5, -11, 
		-2,   3,   0 ,   50,  +3,   4,  -2, -37, 
		 1,  -1,   7 ,    0,  60,   1,   2,  27, 
		 2,   2,  13 ,    7, -31,   70,  -6,   17, 
		 3,  -3,  11 ,   18,  43,   7,   80,  -14, 
		 4,  -4, -17,   -23,  37,  -9,    3,   740
	} ;

	double Q[9],R[9],Ap[9],B[64],L[64],Lt[64];

	long i,j,k=0;



/*	TEST matrice C, 8x8. Ok.... */
	for(i=1; i<=8 ; i++)
		for (j=i+1; j<=8 ; j++)
			mat(C,8,i,j) = mat(C,8,j,i) ;/* On symétrise C */
	afficher_matrice(C,8);
	decomposition_Cholesky(C,8,L);
	afficher_matrice(L,8);

	transpose(L,Lt,8);
	afficher_matrice(Lt,8);

	produit_matrice_matrice(L,Lt,B,8);
	afficher_matrice(B,8);


/*	TEST matrice A, 3x3. Ok....
	afficher_matrice(A,3);
	decomposition_Cholesky(A,3,L);
	afficher_matrice(L,3);

	transpose(L,Lt,3);
	afficher_matrice(Lt,3);

	produit_matrice_matrice(L,Lt,B,3);
	afficher_matrice(B,3);
*/
}

