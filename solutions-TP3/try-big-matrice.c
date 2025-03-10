#include "algebre.c"
#include <math.h>

#define N 400
int main(int argc, char *argv[]){
	int i,j;
	FILE *ifp ;
	double A[N*N], A_inverse[N*N] ;
	double I[N*N], B[N*N] ;

	ifp = fopen("/home/kassa/mozilla.ps" , "r" ) ;
	fread( A , sizeof(double) , N*N , ifp ) ;
	close(ifp);

	for (i=1;i<=N;i++)
		for (j=1;j<=N;j++)
			mat(A,N,i,j) = log ( fabs(mat(A,N,i,j)) ) ;

/*	afficher_matrice( A , N ) ; */

/* for(i=1;i<=10;i++){ */
	inverse_matrice(A,A_inverse,N);
/*	afficher_matrice( A_inverse , N ) ; */

	/* checking results */
	produit_matrice_matrice(A,A_inverse,B,N);
/* } */

/*
	for (i=1;i<=N;i++)
		for (j=1;j<=N;j++)
			if ( mat(B,N,i,j) <= 1.e-10 ) mat(B,N,i,j) = 0. ;
	afficher_matrice( B , N ) ;
*/

}
