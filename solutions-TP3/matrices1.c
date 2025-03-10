#include <stdio.h>
main(){ 
    int i,j,k;

    double x[6]  = { 1. , 3. , 6.75 , 0. , -2. , -1.e-2 } ;
    double B[18] = { 	 1,  -2,  3,  4, 5,   0,
			-1, 2.3,  0, -6, 9,   5,
			10,  20, 31, 42, 55, 89	};
    double u[3], R[9], S[36]; /* u=B.x ; R = B.Bt ; S = Bt.B */

#define B(i,j)	   B[ ((i)-1)*6 + ((j)-1)]
#define x(i) 	   x[ (i)-1 ]
#define R(i,j) 	   R[ ((i)-1)*3 + ((j)-1)]
#define u(i) 	   u[ (i)-1 ]

    /* Calcul et affichage du vecteur u */
    printf("u = " );
    for(i=1 ; i<=3 ; i++) {
	u(i) = 0. ;	for(j=1 ; j<=6 ; j++)	u(i) += B(i,j) * x(j) ;
	printf(" %10g" , u(i) );
    }
    /* Calcul et affichage de la matrice R */
    printf("\n\nR = \n" );
    for(i=1 ; i<=3 ; i++)
      {
	for(j=1 ; j<=3 ; j++) 
	   {
		R(i,j) = 0. ;	for(k=1 ; k<=3 ; k++)	R(i,j)  += B(i,k) * B(j,k)  ;
		printf(" %10g" , R(i,j) ) ;
	   }
	printf("\n");
      }
}
