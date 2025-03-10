#include <stdlib.h>
#include "algebre.c"

void methode_de_la_puissance_directe( double A[] , int n , double x0[] , double eps , double *lambda)
{
	/*
	Cette fonction utilise la m�thode de la puissance directe pour calculer
	la valeur et le vecteur propres de la matrice carr�e A, nxn. Cette valeur propre
	est celle qui a le module le plus grand, et est suppos�e distincte des autres.
	x0 est le vecteur initial avec lequel on d�marre l'it�ration.
	Le r�sultat pour la valeur propre est mis dans (*lambda), et le vecteur 
	propre, normalis� � 1, est remis dans x0[].
	On arr�te l'it�ration lorsque    | lambda^{n+1} - lambda^n | <= eps
	*/

	double lambda_p;
	double x1[n] ;
	int i,k=1; /* k=1: on choisit au d�but la composante x0[1] �gale � 1. On changera plus tard si n�cessaire */
	*lambda = 100000. ;

	for( i=1; i<=n; i++)	vec(x0,i) = mat(A,n,i,1) ;
	

	while ( vec(x0,k) == 0. &&  k<=n ) k++; /* on cherche x[k] != 0. */
	/* scale x0 d'abord, pour avoir x[k]=1 */
	for( i=1; i<=n; i++)	vec(x0,i) /= vec(x0,k) ;

	for( ; ; )
	{
		produit_matrice_vecteur(A,x0,x1,n) ;
		lambda_p = vec(x1,k) ;
		for( i=1; i<=n; i++)	vec(x0,i) = vec(x1,i)/lambda_p ;
		if( fabs(lambda_p- (*lambda) ) < eps ) { *lambda = lambda_p ; break;}
		*lambda = lambda_p;
	}

	/* scale x0 pour avoir une norme �gale � 1 */
	lambda_p=0. ; /* calcul de la norme=lambda_p */
	for(i=1;i<=n;i++)	lambda_p += vec(x0,i)*vec(x0,i) ;
	lambda_p = sqrt( lambda_p );
	for(i=1;i<=n;i++)	vec(x0,i) /= lambda_p ;
	return;
}

void resoudre_avec_PLU(int o[], double L[], double U[], double b[], int n, double x[]){
	/*
	La matrice A est carr�e,  (n x n): b est donn�. on r�soud A . x = b avec la m�thode
	P.A = L.U  . La matrice de permutation P est repr�sent�e par le vecteur o[n]
	(d'entiers).
	*/
	int i,j;

	/* On r�soud L.y = P^-1 . b */
	for (i=1; i<=n ; i++) {
		vec(x,vec(o,i)) = vec(b,vec(o,i)) ;
		for (j=1; j<i ; j++) 	vec(x,vec(o,i)) -= mat(L,n,i,j) * vec(x,vec(o,j)) ;
		/* Ici on ne divise pas par L_ii, car L_ii = 1 pour la d�composition Doolittle */
	}		

	/* ensuite U.x = y */
	for (i=n; i>=1 ; i--) {
		for (j=i+1; j<=n ; j++)	vec(x,vec(o,i)) -= mat(U,n,vec(o,i),j) * vec(x,vec(o,j)) ;
		vec(x,vec(o,i)) /= mat(U,n,vec(o,i),i) ;
	}
	return;
}

void methode_de_la_puissance_inverse( double A[] , int n , double x0[] , double eps , double *lambda)
{
	/*
	Cette fonction utilise la m�thode de la puissance inverse pour calculer
	la valeur et le vecteur propres de la matrice carr�e A, nxn. Cette valeur propre
	est celle qui a le module le plus petit, et est suppos�e distincte des autres.
	x0 est le vecteur initial avec lequel on d�marre l'it�ration.
	Le r�sultat pour la valeur propre est mis dans (*lambda), et le vecteur 
	propre, normalis� � 1, est remis dans x0[].
	On arr�te l'it�ration lorsque    | lambda^{n+1} - lambda^n | <= eps
	*/

	double lambda_p,L[n*n],U[n*n];
	double x1[n] ;
	int i,j,k=1; /* k=1: on choisit au d�but la composante x0[1] �gale � 1. On changera plus tard si n�cessaire */
	int o[n] ;
	*lambda = 100000. ;

	for( i=1; i<=n; i++)	vec(x0,i) = mat(A,n,i,1) ;

	while ( vec(x0,k) == 0. &&  k<=n ) {
		k++; /* on cherche x[k] != 0. */
	}
	/* scale x0 d'abord, pour avoir x[k]=1 */
	for( i=1; i<=n; i++)	vec(x0,i) /= vec(x0,k) ;

	gauss_PLU(A,n,o,L,U); /* on decompose A en PLU */


	for( ; ; )
	{
	/* Au lieu de calculer A^-1, puis calculer A^-1 . x0 = x1 ,
	on r�soud avec PLU: A.x1 = x0, avec x0 donn� et x1 le vecteur inconnu.
	*/
		resoudre_avec_PLU(o, L, U, x0, n, x1);
		lambda_p = vec(x1,k) ;
		for( i=1; i<=n; i++)	vec(x0,i) = vec(x1,i)/lambda_p ;
		if( fabs(lambda_p- (*lambda) ) < eps ) { *lambda = lambda_p ; break;}
		*lambda = lambda_p;
	}

	*lambda = 1./(*lambda) ; /* lambda est la valeur propre de A */

	/* scale x0 pour avoir une norme �gale � 1 */
	lambda_p=0. ; /* calcul de la norme=lambda_p */
	for(i=1;i<=n;i++)	lambda_p += vec(x0,i)*vec(x0,i) ;
	lambda_p = sqrt( lambda_p );
	for(i=1;i<=n;i++)	vec(x0,i) /= lambda_p ;
	return;
}

void methode_de_la_puissance_shifted( double A[] , int n , double x0[] , double eps , double *lambda)
{
	/*
	Cette fonction utilise la m�thode de la puissance shifted, combin�e avec
	la m�thode inverse pour calculer la valeur et le vecteur propres de 
	la matrice carr�e A, nxn. Cette valeur propre est souvent une valeur interm�diaire.
	si B = A - lambda_s I , alors A et B ont les m�mes vecteurs propres. et 
	lambda_A =  lambda_B + lambda_s. 
	*/

	double B[n*n], lambda_s;
	double x1[n] ;
	int i,j,k=1; /* k=1: on choisit au d�but la composante x0[1] �gale � 1. On changera plus tard si n�cessaire */
	int o[n] ;
	char sss[20];

	for( i=1; i<=n; i++)	vec(x0,i) = mat(A,n,i,1) ;

	copy_matrice(A,B,n);

	sprintf(sss,"%2g", *lambda);
	sscanf(sss,"%lf", &lambda_s);

	/* lambda_s = (*lambda) */; /* lambda_s est initialement dans (*lambda): BAD idea! */
	for (i=1 ; i<=n ; i++)	mat(B,n,i,i) -= lambda_s ; /* B = A -lambda_s I */

	methode_de_la_puissance_inverse( B, n , x0, eps , lambda) ;
	*lambda += lambda_s ; /* *lambda est maintenant �gale � la valeur propre de A */
	return;
}

void decomposition_QR( double A[], int n, double Q[], double R[])
{
	/* 
	Cette fonction d�compose une matrice A carr�e, r�elle, inversible, en Q.R
	avec Q orthogonale, et R triangulaire sup�rieure.
		A = [\vec a_1  \vec a_2 ...]  ===>  A_ij = (\vec a_j)_i
		Q = [\vec q_1  \vec q_2 ...]  ===>  Q_ij = (\vec q_j)_i
	Les �quations sont:
		pour i=1,...,n
			1. Calculer $\vec a'_i = \vec a_i 
						- \sum_{k=1}^{i-1} (\vec a_i  \cdot \vec q_k) q_k $	
			2. $ R_{ii} = | \vec a'_i | $ (norme)
			3. $ \vec q_i = \vec a'_i / R_{ii} $
			4. $ R_{ij} = \vec q_i  \cdot  \vec a_j $  , pour $j=i+1,...,n$

	L'algorithme est donc:
		Remplir de 0, la partie inf�rieure de R,
		puis it�rer 1,2,3,4 ci-dessous (i=1,...,n):
		1. Calculer $\vec a'_i = \vec a_i 
						- \sum_{k=1}^{i-1} R_{ki} \vec q_k $
			et mettre ce r�sultat interm�diaire dans $\vec q_i$, c'est � dire
			utiliser l'�quation: Q_ji = A_ji - sum_{k=1}^{i-1} R_ki Q_jk pour j=1,...,n
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
		triangulaire sup�rieure.
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

void v_propres( double A[], int n, int ITERATIONS, double eps, double lambda[], double S[]){
	double B[n*n], x[n];
	int i,j;
	methode_QR(A,n,B,ITERATIONS); /* On suppose que les valeurs propres approch�es sont les B_ii */

	for (i=1; i<=n ; i++){
		lambda[i-1] = mat(B,n,i,i);
		methode_de_la_puissance_shifted( A, n , x, eps , &lambda[i-1] ) ;
		for (j=1; j<= n; j++)	mat(S,n,j,i) = vec(x,j) ;
	}
	return;
}

void verifier_v_propres( double A[], int n, double lambda[], double S[], double B[]){
	int i,j;
	produit_matrice_matrice(A,S,B,n);
	for (j=1; j<=n ; j++){/* pour chaque vecteur propre (chaque colonne dans S)
				, on divise par la valeur propre */
		for (i=1; i<=n ; i++)	mat(B,n,i,j) /= vec(lambda,j) ;
	}
	return;
}

void main(){

	double A[] = { /* matrice 3 x 3 */
		 8,  -2,  -2 , 
		-2,   4,  -2 , 
		-2,  -2,  13
	} ;

	double aaa[] ={ /* vecteur � 3 composantes */
		 1, 
		 1, 
		 1
	};

	double a[] ={ /* vecteur � 3 composantes */
		 0, 
		 0, 
		 1
	};

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
	double c[] ={ /* vecteur � 8 composantes */
		 1, 
		 3, 
		 5, 
		-1, 
		23, 
		13, 
		-67, 
		99
	};


	double eps=1.e-10;
	int ITERATIONS=30 ;


	{/* Valeurs et vecteurs propres pour C, 8x8 */
	double lambda[8], S[64], B[64];
	v_propres( C, 8 , ITERATIONS , eps, lambda,  S);

	afficher_v(lambda,8,13);
	afficher_m(S,8,13);

	verifier_v_propres( C, 8, lambda, S, B);
	printf("\nV�rification:\n");
	afficher_m(B,8,13);
	}

	{/* Valeurs et vecteurs propres pour A, 3x3 */
	double lambda[3], S[9], B[9];	
	v_propres( A, 3 , ITERATIONS , eps, lambda,  S);

	afficher_v(lambda,3,10);
	afficher_m(S,3,10);

	verifier_v_propres( A, 3, lambda, S, B);
	printf("\nV�rification:\n");
	afficher_m(B,3,10);

	/* RESULTATS:
	       13.87058512,
	       8.620434078,
	       2.508980799
		};
		{
	      0.2774924805,      0.8690961095,      0.4094751221,
	      0.1364648562,     -0.4575483834,      0.8786505676,
	     -0.9509864701,      0.1879399618,      0.2455673114
		};
	*/
	}

	{/* Valeurs et vecteurs propres pour A, 4x4 pour EMD 3 eme ann�e */
	double lambda[4], S[4], B[4];	
	double A[] = { /* matrice 4 x 4 */
		 1,   2,   3,  4 ,
		 2,   1,  -2, -7 ,
		 3,  -2,   5 , 6 ,
		 4,  -7,   6 , 1
	} ;

	v_propres( A, 4 , ITERATIONS , eps, lambda,  S);

	afficher_v(lambda,4,10);
	afficher_m(S,4,10);

	verifier_v_propres( A, 4, lambda, S, B);
	printf("\nV�rification:\n");
	afficher_m(B,4,10);

	/* RESULTATS:
	{
       13.50806282,
       -8.23613497,
       3.597545738,
     -0.8694735865
	};
	{
      0.2817097916,      0.3837243353,      0.6371032896,      0.6062133504,
     -0.3969137097,      -0.583653571,      0.6879116209,     -0.1690728183,
       0.623856426,      0.1427122831,      0.3080210096,     -0.7039598154,
      0.6114754779,     -0.7012398647,     -0.1612453341,      0.3291813091
	};
	*/
	}



	{/* Valeurs et vecteurs propres pour C, 8x8 */
	double lambda[8], S[64], B[64];
	v_propres( C, 8 , ITERATIONS , eps, lambda,  S);

	afficher_v(lambda,8,10);
	afficher_m(S,8,10);

	verifier_v_propres( C, 8, lambda, S, B);
	printf("\nV�rification:\n");
	afficher_m(B,8,10);
	}


	/* Valeurs et vecteurs propres pour C, en utilisant Householder */
	{
	double lambda[8], S[64], B[64], E[64];
	int i,j;

	v_propres( E, 8 , ITERATIONS , eps, lambda,  S);

	afficher_v(lambda,8,10);
	afficher_m(S,8,10);

	verifier_v_propres( E, 8, lambda, S, B);
	printf("\nV�rification Hessenberg:\n");
	afficher_m(B,8,10);
	}


	{/* Valeurs et vecteurs propres pour C, que l'on aura sym�triser avant */
	double lambda[8], S[64], B[64],sum;
	int i,j,k;
	for (i=1; i<=8 ; i++)
		for (j=i+1; j<=8 ; j++)	mat(C,8,i,j) = mat(C,8,j,i) ;


	v_propres( C, 8 , ITERATIONS , eps, lambda,  S);

	afficher_v(lambda,8,10);
	afficher_m(S,8,10);

	verifier_v_propres( C, 8, lambda, S, B);
	printf("\nV�rification:\n");
	afficher_m(B,8,10);


	printf("\ndet S = %g \n", determinant(C,8));
	for (i=1; i<=8 ; i++){/* On calcule les produits scalaires entre les vecteurs propres x_i . x_j
				Nous devons trouver \delta_ij */
		for (j=i; j<=8 ; j++){
			sum=0. ;
			for (k=1; k<=8 ; k++)	sum += mat(S,8,k,i) * mat(S,8,k,j)  ;
			printf("x[%d] . x[%d] = %g\n" , i, j, sum );
		}
	}
	}



    { /* Valeurs propre pour A: 3x3 */
	double x0[3],x1[3], lambda;
	int i, k=0;

	copy_vecteur(a,x0,3);

	methode_de_la_puissance_inverse( A, 3 , x0, eps , &lambda) ;
	printf("\n lambda = %20.6g\n x0 = "  , lambda);

	afficher_vecteur(x0,3);

	/* V�rifions */
	produit_matrice_vecteur(A,x0,x1,3);
	for(i=1;i<=3;i++)	vec(x1,i) /= lambda ;
	afficher_vecteur(x1,3);

    }

    { /* Valeurs propre pour A: shifted by 13.87 */
	double x0[3],x1[3], lambda;
	int i, k=0;

	copy_vecteur(a,x0,3);
	lambda = 13.87 ;
	methode_de_la_puissance_shifted( A, 3 , x0, eps , &lambda) ;
	printf("\n lambda = %20.6g\n x0 = "  , lambda);

	afficher_vecteur(x0,3);

	/* V�rifions */
	produit_matrice_vecteur(A,x0,x1,3);
	for(i=1;i<=3;i++)	vec(x1,i) /= lambda ;
	afficher_vecteur(x1,3);

    }

    { /* Valeurs propre pour A: shifted by 8.6 */
	double x0[3],x1[3], lambda;
	int i, k=0;

	copy_vecteur(a,x0,3);
	lambda = 8.6 ;
	methode_de_la_puissance_shifted( A, 3 , x0, eps , &lambda) ;
	printf("\n lambda = %20.6g\n x0 = "  , lambda);

	afficher_vecteur(x0,3);

	/* V�rifions */
	produit_matrice_vecteur(A,x0,x1,3);
	for(i=1;i<=3;i++)	vec(x1,i) /= lambda ;
	afficher_vecteur(x1,3);

    }

    { /* Valeurs propre pour C: 8x8 */
	double x0[8],x1[8], lambda;
	int i, k=0;

	copy_vecteur(c,x0,8);

	methode_de_la_puissance_inverse( C, 8 , x0, eps , &lambda) ;
	printf("\n lambda = %15.10g\n x0 = "  , lambda);
	afficher_vecteur(x0,8);
	/* V�rifions */
	produit_matrice_vecteur(C,x0,x1,8);
	for(i=1;i<=8;i++)	vec(x1,i) /= lambda ;
	afficher_vecteur(x1,8);
    }

	return;
}

