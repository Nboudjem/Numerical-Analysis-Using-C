#include "algebre.c"
#include "matrices.txt"

main(){

	double a_prime[3]; /* Nous déclarons les vecteurs a' , b' , et c' */
	double b_prime[4];
	double c_prime[8];


/* ***************************************************************************/
	/* Réponse à la question 1 */
/* ***************************************************************************/
	printf("\ndouble A[] =\n");	afficher_matrice(A,3) ;
	printf("\ndouble B[] =\n");	afficher_matrice(B,4) ;
	printf("\ndouble C[] =\n");	afficher_matrice(C,8) ;


/* ***************************************************************************/
	/* Réponse à la question 2 */
/* ***************************************************************************/
	printf("\ndouble a[] =\n");	afficher_vecteur(a,3) ;
	printf("\ndouble b[] =\n");	afficher_vecteur(b,4) ;
	printf("\ndouble c[] =\n");	afficher_vecteur(c,8) ;


/* ***************************************************************************/
	/* Réponse à la question 3 */
/* ***************************************************************************/
	produit_matrice_vecteur(A,a,a_prime,3) ;
	printf("\n a' = A . a =\n");		afficher_vecteur(a_prime,3) ;

	produit_matrice_vecteur(B,b,b_prime,4) ;
	printf("\n b' = B . b =\n");		afficher_vecteur(b_prime,4) ;

	produit_matrice_vecteur(C,c,c_prime,8) ;
	printf("\n c' = C . c =\n");		afficher_vecteur(c_prime,8) ;


/* ***************************************************************************/
	/* Réponse à la question 4. */
/* ***************************************************************************/
	/* Nous ne calculons ici que A.At et At.A.
	   Les autres produits, c'est pareil
	*/
	{
		double At[9], A_At[9]; /* Nous déclarons les matrices: At est la matrice transposée */
		int i,j;
		for(i=1;i<=3;i++){
			for(j=1;j<=3;j++) { mat(At,3,i,j) = mat(A,3,j,i); } /* on definit la matrice At */
		}

		transpose(A,At,3);	
		produit_matrice_matrice( A , At , A_At , 3 ) ;
		printf("\n A . At =\n");	afficher_matrice(A_At,3) ;

		produit_matrice_matrice( At , A , A_At , 3 ) ; /* ici A_At = At.A */
		printf("\n At . A =\n");	afficher_matrice(A_At,3) ;
	}

/* ***************************************************************************/
	/* réponse à la question 5 */
/* ***************************************************************************/
	/* Calcul du determinant de A, de B et de C. */
	{
		double z ;

		z = determinant(A,3);	
		printf("\n det A = %10g \n" , z );

		z = determinant(B,4);	
		printf("\n det B = %10g \n" , z );

		z = determinant(C,8);	
		printf("\n det C = %10g \n" , z );
	}

/* ***************************************************************************/
	/* réponse à la question 6 */
/* ***************************************************************************/
	/* Nous résolvons uniquement " A.x = a ". Le même code peut être
	   utilisé pour résoudre  B.x = b   et  C.x = c
	*/
	{
	double x[3];
	resoudre_avec_gauss(A,a,x,3);
	printf("\ndouble x[] =\n");	afficher_vecteur(x,3) ;
	/* Verification */
	produit_matrice_vecteur(A,x,a_prime,3) ;
	printf("\n a'  =\n");		afficher_vecteur(a_prime,3) ;
	}


/* ***************************************************************************/
	/* réponses à la question 7 */
/* ***************************************************************************/
	{ 
	double CCC[64], III[64];
	matrice_identite(CCC,8);
	gauss_multi(C,8,CCC,8); /* CCC contient maintenant l'inverse de C */
	afficher_matrice(CCC,8);

	/* VERIFIONS que C . CCC = CCC . C = I (8x8) */
	produit_matrice_matrice(C,CCC,III,8);
	afficher_matrice(III,8); /* resultat: III très proche de l'identité
				à cause des erreurs d'arrondis. */

		{ 	/*Nous rectifions maintenant: */
			int i,j;
			for(i=1;i<=8;i++)
				for(j=1;j<=8;j++)
					if ( mat(III,8,i,j) <= 1.e-12 ) 
						mat(III,8,i,j) =0. ;
			afficher_matrice(III,8);
		}

	produit_matrice_matrice(CCC,C,III,8);
	afficher_matrice(III,8);
		{ 	/*Nous rectifions maintenant: */
			int i,j;
			for(i=1;i<=8;i++)
				for(j=1;j<=8;j++)
					if ( mat(III,8,i,j) <= 1.e-12 ) 
						mat(III,8,i,j) =0. ;
			afficher_matrice(III,8);
		}

	}
/* ***************************************************************************/
	/* réponses à la question 7 . VERSION 2. */
/* ***************************************************************************/
	{ 
	double CCC[64], III[64];
	inverse_matrice(C,CCC,8);
	afficher_matrice(CCC,8);
	/* VERIFIONS que C . CCC = I (8x8) */
	produit_matrice_matrice(C,CCC,III,8);
	afficher_matrice(III,8);
	}


}

