/* ***************************************************************************/
/* ***************************************************************************/

#include "algebre.c"
#include <math.h>

void decomposition_QR( double A[], int n, double Q[], double R[])
{
/* 
   Cette fonction décompose une matrice A carrée, réelle, inversible, en Q.R
   avec Q orthogonale, et R triangulaire supérieure.
	A = [\vec a_1  \vec a_2 ...]  ===>  A_ij = (\vec a_j)_i
	Q = [\vec q_1  \vec q_2 ...]  ===>  Q_ij = (\vec q_j)_i
   Les équations sont:
	pour i=1,...,n
		1. Calculer $\vec a'_i = \vec a_i 
					- \sum_{k=1}^{i-1} (\vec a_i  \cdot \vec q_k) q_k $	
		2. $ R_{ii} = | \vec a'_i | $ (norme)
		3. $ \vec q_i = \vec a'_i / R_{ii} $
		4. $ R_{ij} = \vec q_i  \cdot  \vec a_j $  , pour $j=i+1,...,n$

  L'algorithme est donc:
     Remplir de 0, la partie inférieure de R,
     puis itérer 1,2,3,4 ci-dessous (i=1,...,n):
	1. Calculer $\vec a'_i = \vec a_i 
					- \sum_{k=1}^{i-1} R_{ki} \vec q_k $
		et mettre ce résultat intermédiaire dans $\vec q_i$, c'est à dire
		utiliser l'équation: Q_ji = A_ji - sum_{k=1}^{i-1} R_ki Q_jk pour j=1,...,n
	2. Calculer la norme de a'_i  ===>  R_ii = SQRT( sum_{j=1}^{n} Q_ij^2 )
	3. normaliser q_i avec Q_ji /= R_ii pour j=1,...,n
	4. Calculer R_ij = sum_{k=1}^{n} Q_ki A_kj   pour j=1,...,i-1
*/

	int i,j,k,m;
	/* Etape 0: R_ij = 0 si i>j */
	for(i=1; i<=n ; i++)
		for (j=1; j<i ; j++)
			mat(R,n,i,j) = 0. ;

	for(i=1; i<=n ; i++){
		/* Etape 1: Calcul de a'_i . On met ce vecteurs temporairement dans q_i */
		/* Utilisons: Q_ji = A_ji - sum_{k=1}^{i-1} R_ki Q_jk pour j=1,...,n */
		for (j=1; j<=n ; j++){
			mat(Q,n,j,i) = mat(A,n,j,i) ;
			for (k=1; k<i ; k++)	mat(Q,n,j,i) -= mat(R,n,k,i)*mat(Q,n,j,k) ;
		}

		/* Etape 2: Calcul de la norme de a'_i = R_ii = sum_{j=1}^{n} Q_ij^2 */
		mat(R,n,i,i) = 0. ;
		for (j=1; j<=n ; j++)	mat(R,n,i,i) += mat(Q,n,j,i)*mat(Q,n,j,i) ;
		mat(R,n,i,i) = sqrt( mat(R,n,i,i) ) ;

		/* Etape 3: On normalise q_i avec avec Q_ji /= R_ii */
		for (j=1; j<=n ; j++)	mat(Q,n,j,i) /= mat(R,n,i,i) ;

		/* Etape 4: Calculer R_ij = sum_{k=1}^{n} Q_ki A_kj   pour j=1,...,i-1 */
		for (j=i+1; j<=n ; j++){
			mat(R,n,i,j) = 0. ;
			for (k=1; k<=n ; k++)	mat(R,n,i,j) += mat(Q,n,k,i)*mat(A,n,k,j) ;
		}
	}
return;
}

void methode_QR( double A[], int n, double B[], int ITERATIONS)
{
/* 
	A = Q.R
	B = R.Q ===> B =  Q^t . A . Q
	Transforme A en B. Au bout de INTERATIONS fois, B est quasi
	triangulaire supérieure.
*/
	int i;
	double Q[n*n],R[n*n];

	copy_matrice(A,B,n);
	for(i=1; i<= ITERATIONS ; i++){
		decomposition_QR(B,n,Q,R);
		produit_matrice_matrice(R,Q,B,n);
	}
return;
}


main(){
	double A[] = { /* matrice 3 x 3 */
		 8,  -2,  -2 , 
		-2,   4,  -2 , 
		-2,  -2,  13
	} ;
	double B[9];
	double C1[] = { /* matrice 8 x 8 */
		 2,   3,  -1 ,    4,   1,  -2,   5,  32, 
		 1,  30,   5 ,    3,  -2,   3,   0,  10, 
		 6,  -3,  40 ,   -1,   1,   2,  -5, -11, 
		-2,   3,   0 ,   50,  +3,   4,  -2, -37, 
		 1,  -1,   7 ,    0,  60,   1,   2,  27, 
		 2,   2,  13 ,    7, -31,   70,  -6,   17, 
		 3,  -3,  11 ,   18,  43,   7,   80,  -14, 
		 4,  -4, -17,   -23,  37,  -9,    3,   740
	} ;
	double C[] = { /* matrice 8 x 8 */
		 2,   3,  -1 ,    4,   1,  -2,   5,  32, 
		 1,  -2,   5 ,    3,  -2,   3,   0,  10, 
		 6,  -3,  15 ,   -1,   1,   2,  -5, -11, 
		-2,   3,   0 ,    1,  +3,   4,  -2, -37, 
		 1,  -1,   7 ,    0,  13,   1,   2,  27, 
		 2,   2,  13 ,    7, -31,   4,  -6,   17, 
		 3,  -3,  11 ,   18,  43,   7,  -8,  -14, 
		 4,  -4, -17,   -23,  37,  -9,   3,   15
	} ;
	long i,j,k=0;
	double CCC[64];

	double D[] = { /* matrice 4 x 4 */
		2,   3,  -1 ,   4, 
		1,  -2,   5 ,   3, 
		6,  -3,  15 ,  -1, 
		-2,  3,   0 ,   1
	} ;
	double DDD[16];

	methode_QR(A,3,B,100);
	afficher_matrice(B,3);


	printf("\n det C = %g\n" , determinant(C,8) );
	methode_QR(C,8,CCC,100);
	afficher_matrice(CCC,8);

	printf("\n det D = %g\n" , determinant(D,4) );
	methode_QR(D,4,DDD,400);
	afficher_matrice(DDD,4);
}

