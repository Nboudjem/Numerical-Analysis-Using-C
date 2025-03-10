#include "algebre.c"
double A[] = { /* matrice 3 x 3 */
	2,   3,  -1 , 
	1,  -2,   5 , 
	6,  -3,  15
} ;

double B[] = { /* matrice 4 x 4 */
	2,   3,  -1 ,   4, 
	1,  -2,   5 ,   3, 
	6,  -3,  15 ,  -1, 
	-2,  3,   0 ,   1
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

double a[] ={ /* vecteur � 3 composantes */
	 1, 
	 3, 
	 5
};

double b[] ={ /* vecteur � 4 composantes */
	 1, 
	 3, 
	 5, 
	-1
};

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
	double x1[n] ; /* le vecteur x1 correspond au vecteur y utilis� dans le cours et le TD*/
	int i,k; 

	for (k=1; k<= n; k++)
	   {	/* on choisit d'abord x0 = (1,0,...,0) , et si �a ne marche
		pas, on choisit x0=(0,1,0,...) , etc .... */
		for (i=1;i<=n;i++)
		   {
			if(i==k)	vec(x0,i) = 1; /* ce vecteur x[0] a uniquement
							une composante non nulle: x0[k]=1 */
			else		vec(x0,i) = 0;
		   }

		lambda_p = 100000. ; /* Notre notation est:	*lambda --> ancien
								lambda_p--> nouveau
					*/  
		do
		    {
			*lambda = lambda_p;
			produit_matrice_vecteur(A,x0,x1,n) ;
			lambda_p = vec(x1,k) ; /* On esp�re ici que lambda_p est diff�rent de 0.
					Sinon, nous avons un probl�me. Nous devons changer de vecteur
					initial x0, et tout recommencer. */
			if (lambda_p == 0. ){ break; /* on sort de la boucle "do" */}
			for( i=1; i<=n; i++)	vec(x0,i) = vec(x1,i)/lambda_p ;
		    }
		while( fabs(lambda_p- (*lambda) ) > eps ) ;


		if( fabs(lambda_p- (*lambda) ) < eps ) {/* nous avons une solution */
			*lambda = lambda_p;
			break; /* on sort de la boucle for(k=1;k<=n;k++) */
		}
		else { /* On n'a pas encore de solution. laissons k s'incr�m�nter de 1
			et recommen�ons avec un nouveau vecteur */
			continue;
		}
	   } /* fin de la boucle for(k=1;k<=n;k++)


	/* Arriv� ici, nous avons un vecteur propre x0 qu'on normalise x0 � 1 */
	lambda_p=0. ; /* calcul de la norme de x0 , appel�e ici "lambda_p". 
			On r�-utilise cette variable "lambda_p" 
			puisque l'on en a plus besoin */
	for(i=1;i<=n;i++)	lambda_p += vec(x0,i)*vec(x0,i) ;
	lambda_p = sqrt( lambda_p );
	for(i=1;i<=n;i++)	vec(x0,i) /= lambda_p ;
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

	/* Cette m�thode est �quivalente � la methode de la puisance directe
	pour la matrice inverse A^-1. Le vecteur propre obtenu est le m�me
	pour A et A^-1, mais la valeur propre de A est �gale a 1/(valeur propre de A^-1) */

	double A_i[n*n]; /* A_i == A^-1 */
	inverse_matrice(A,A_i,n,2) ; /* Si on augmente le dernier argument, on
					augmente la pr�cision dans le calcul de A^-1 */
	methode_de_la_puissance_directe(A_i,n,x0,eps,lambda);
	*lambda = 1./(*lambda);
	return;
}


/*void methode_de_la_puissance_avec_translation( double A[] , int n , double x0[] , double eps , double *lambda)*/
void methode_de_la_puissance_shifted( double A[] , int n , double lambda_s, double x0[] , double eps , double *lambda)
{
	/*
	Cette fonction utilise la m�thode de la puissance shifted, combin�e avec
	la m�thode inverse pour calculer la valeur et le vecteur propres de 
	la matrice carr�e A, nxn. Cette valeur propre est souvent une valeur interm�diaire.
	si B = A - lambda_s I , alors A et B ont les m�mes vecteurs propres. et 
	lambda_A =  lambda_B + lambda_s. 
	*/

	double B[n*n];
	int i,j;

	/* On calcule B = A - lamba_s I */
	for (i=1 ; i<=n ; i++){
		for (j=1 ; j<=n ; j++){
			if (i==j)	mat(B,n,i,j) = mat(A,n,i,j) - lambda_s ; 
			else 		mat(B,n,i,j) = mat(A,n,i,j) ; 
		}
	}

	methode_de_la_puissance_inverse( B, n , x0, eps , lambda) ;
	*lambda += lambda_s ;
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
						- \sum_{k=1}^{i-1} (\vec a_i  \cdot \vec q_k) \vec q_k $	
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
	/* Etape 0: R_ij = 0 si i>j car R triangulaire sup�rieure */
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


void main(){

	double eps=1.e-12 ;
	double z,lambda, At[3*3], As[3*3], Bt[4*4], Bs[4*4], Ct[8*8], Cs[8*8];
	int i,k,j;
	
	transpose(A,At,3);
	produit_matrice_matrice(A,At,As,3);
	z = determinant(As,3);
	printf("\ndet As = %.6g\n" , z);

	transpose(B,Bt,4);
	produit_matrice_matrice(B,Bt,Bs,4);
	z = determinant(Bs,4);
	printf("\ndet Bs = %.6g\n" , z);

	transpose(C,Ct,8);
	produit_matrice_matrice(C,Ct,Cs,8);
	z = determinant(Cs,8);
	printf("\ndet Cs = %.6g\n" , z);



	printf("\n-----------------------------------------------------------");
	printf("\nREPONSE A LA QUESTION 2");
	printf("\n-----------------------------------------------------------\n\n");


	/* **************************************************************************** */
    {
	printf("\n\n\n");
	printf("-----------------------------------------------------------");
	printf("\nValeur propre pour As: 3x3 avec methode directe \n");
	double x0[3], x1[3];

	methode_de_la_puissance_directe( As, 3 , x0, eps , &lambda) ;
	printf("\n lambda = %20.14g\n x0 = "  , lambda);

	afficher_v(x0,3,14);

	/* V�rifions */
	produit_matrice_vecteur(As,x0,x1,3);
	for(i=1;i<=3;i++)	vec(x1,i) /= lambda ;
	printf("\n Verification: a-t-on As . x0  = lambda x0?");
	printf("\n\tlambda = %.14g\n As . x0 / lambda = "  , lambda);
	afficher_v(x1,3,14);

    }



	/* **************************************************************************** */
    {
	printf("\n\n\n");
	printf("-----------------------------------------------------------");
	printf("\n Valeur propre pour Bs: 4x4 avec methode directe \n");
	double x0[4], x1[4];

	methode_de_la_puissance_directe( Bs, 4 , x0, eps , &lambda) ;
	printf("\n lambda = %20.14g\n x0 = "  , lambda);

	afficher_v(x0,4,14);

	/* V�rifions */
	produit_matrice_vecteur(Bs,x0,x1,4);
	for(i=1;i<=4;i++)	vec(x1,i) /= lambda ;
	printf("\n Verification: a-t-on Bs . x0  = lambda x0?");
	printf("\n\tlambda = %.14g\n Bs . x0 / lambda = "  , lambda);
	afficher_v(x1,4,14);

    }


	/* **************************************************************************** */
    {
	printf("\n\n\n");
	printf("-----------------------------------------------------------");
	printf("\n Valeur propre pour Cs: 8x8 avec methode directe \n");
	double x0[8], x1[8];

	methode_de_la_puissance_directe( Cs, 8 , x0, eps , &lambda) ;
	printf("\n lambda = %20.14g\n x0 = "  , lambda);

	afficher_v(x0,8,14);

	/* V�rifions */
	produit_matrice_vecteur(Cs,x0,x1,8);
	for(i=1;i<=8;i++)	vec(x1,i) /= lambda ;
	printf("\n Verification: a-t-on Cs . x0  = lambda x0?");
	printf("\n\tlambda = %.14g\n Cs . x0 / lambda =  "  , lambda);
	afficher_v(x1,8,14);

    }


	printf("\n-----------------------------------------------------------");
	printf("\nREPONSE A LA QUESTION 3");
	printf("\n-----------------------------------------------------------\n\n");


	/* **************************************************************************** */
    {
	printf("\n\n\n");
	printf("-----------------------------------------------------------");
	printf("\nValeur propre pour As: 3x3 avec methode inverse \n");
	double x0[3], x1[3];

	methode_de_la_puissance_inverse( As, 3 , x0, eps , &lambda) ;
	printf("\n lambda = %20.14g\n x0 = "  , lambda);

	afficher_v(x0,3,14);

	/* V�rifions */
	produit_matrice_vecteur(As,x0,x1,3);
	for(i=1;i<=3;i++)	vec(x1,i) /= lambda ;
	printf("\n Verification: a-t-on As . x0  = lambda x0?");
	printf("\n\tlambda = %.14g\n As . x0 / lambda = "  , lambda);
	afficher_v(x1,3,14);

    }

	/* **************************************************************************** */
    {
	printf("\n\n\n");
	printf("-----------------------------------------------------------");
	printf("\n Valeur propre pour Bs: 4x4 avec methode inverse \n");
	double x0[4], x1[4];

	methode_de_la_puissance_inverse( Bs, 4 , x0, eps , &lambda) ;
	printf("\n lambda = %20.14g\n x0 = "  , lambda);

	afficher_v(x0,4,14);

	/* V�rifions */
	produit_matrice_vecteur(Bs,x0,x1,4);
	for(i=1;i<=4;i++)	vec(x1,i) /= lambda ;
	printf("\n Verification: a-t-on Bs . x0  = lambda x0?");
	printf("\n\tlambda = %.14g\n Bs . x0 / lambda = "  , lambda);
	afficher_v(x1,4,14);

    }


	/* **************************************************************************** */
    {
	printf("\n\n\n");
	printf("-----------------------------------------------------------");
	printf("\n Valeur propre pour Cs: 8x8 avec methode inverse \n");
	double x0[8], x1[8];

	methode_de_la_puissance_inverse( Cs, 8 , x0, eps , &lambda) ;
	printf("\n lambda = %20.14g\n x0 = "  , lambda);

	afficher_v(x0,8,14);

	/* V�rifions */
	produit_matrice_vecteur(Cs,x0,x1,8);
	for(i=1;i<=8;i++)	vec(x1,i) /= lambda ;
	printf("\n Verification: a-t-on Cs . x0  = lambda x0?");
	printf("\n\tlambda = %.14g\n Cs . x0 / lambda =  "  , lambda);
	afficher_v(x1,8,14);

    }

	printf("\n-----------------------------------------------------------");
	printf("\nREPONSE A LA QUESTION 4. ON TESTE la methode avec translation");
	printf("\n-----------------------------------------------------------\n\n");


	/* **************************************************************************** */
    {
	printf("\n\n\n");
	printf("-----------------------------------------------------------");
	printf("\nValeur propre pour As: 3x3 avec methode avec translation. On utilise s = 13.8 \n");
	double x0[3], x1[3], s=13.8;

	methode_de_la_puissance_shifted( As, 3 , s , x0, eps , &lambda) ;
	printf("\n lambda = %20.14g\n x0 = "  , lambda);

	afficher_v(x0,3,14);

	/* V�rifions */
	produit_matrice_vecteur(As,x0,x1,3);
	for(i=1;i<=3;i++)	vec(x1,i) /= lambda ;
	printf("\n Verification: a-t-on As . x0  = lambda x0?");
	printf("\n\tlambda = %.14g\n As . x0 / lambda = "  , lambda);
	afficher_v(x1,3,14);

    }

	/* **************************************************************************** */
    {
	printf("\n\n\n");
	printf("-----------------------------------------------------------");
	printf("\n Valeur propre pour Bs: 4x4 avec methode avec translation. On utilise s = 13.8  \n");
	double x0[4], x1[4], s=13.8;

	methode_de_la_puissance_shifted( Bs, 4 , s, x0, eps , &lambda) ;
	printf("\n lambda = %20.14g\n x0 = "  , lambda);

	afficher_v(x0,4,14);

	/* V�rifions */
	produit_matrice_vecteur(Bs,x0,x1,4);
	for(i=1;i<=4;i++)	vec(x1,i) /= lambda ;
	printf("\n Verification: a-t-on Bs . x0  = lambda x0?");
	printf("\n\tlambda = %.14g\n Bs . x0 / lambda = "  , lambda);
	afficher_v(x1,4,14);

    }


	/* **************************************************************************** */
    {
	printf("\n\n\n");
	printf("-----------------------------------------------------------");
	printf("\n Valeur propre pour Cs: 8x8 avec methode avec translation. On utilise s = 13.8  \n");
	double x0[8], x1[8], s=13.8;

	methode_de_la_puissance_shifted( Cs, 8 , s, x0, eps , &lambda) ;
	printf("\n lambda = %20.14g\n x0 = "  , lambda);

	afficher_v(x0,8,14);

	/* V�rifions */
	produit_matrice_vecteur(Cs,x0,x1,8);
	for(i=1;i<=8;i++)	vec(x1,i) /= lambda ;
	printf("\n Verification: a-t-on Cs . x0  = lambda x0?");
	printf("\n\tlambda = %.14g\n Cs . x0 / lambda =  "  , lambda);
	afficher_v(x1,8,14);

    }




	printf("\n-----------------------------------------------------------");
	printf("\nREPONSE A LA QUESTION 5. On decompose A,B, et C en produit QR");
	printf("\n-----------------------------------------------------------\n\n");

    {
	double Q[3*3], R[3*3];
	printf("\n\n-----------------------------------------------------------");
	decomposition_QR(A,3,Q,R);
	printf("\nA = Q . R avec :");
	printf("\nQ = ");	afficher_m(Q,3,14);
	printf("\nR = ");	afficher_m(R,3,14);

	/* V�rifions les d�terminants */
	printf("\n det(A) = %.14g" , determinant(A,3) );
	printf("\n det(R) = %.14g" , determinant(R,3) );

    }

    {
	double Q[4*4], R[4*4];
	printf("\n\n-----------------------------------------------------------");
	decomposition_QR(B,4,Q,R);
	printf("\nB = Q . R avec :");
	printf("\nQ = ");	afficher_m(Q,4,14);
	printf("\nR = ");	afficher_m(R,4,14);

	/* V�rifions les d�terminants */
	printf("\n det(B) = %.14g" , determinant(B,4) );
	printf("\n det(R) = %.14g" , determinant(R,4) );

    }

    {
	double Q[8*8], R[8*8];
	printf("\n\n-----------------------------------------------------------");
	decomposition_QR(C,8,Q,R);
	printf("\nB = Q . R avec :");
	printf("\nQ = ");	afficher_m(Q,8,14);
	printf("\nR = ");	afficher_m(R,8,14);

	/* V�rifions les d�terminants */
	printf("\n det(C) = %.14g" , determinant(C,8) );
	printf("\n det(R) = %.14g" , determinant(R,8) );

    }




	printf("\n-----------------------------------------------------------");
	printf("\nREPONSE A LA QUESTION 6. Methode QR pour As, Bs, et Cs");
	printf("\n-----------------------------------------------------------\n\n");



	k = 60 ; /* nombre d'i�trations pour la m�thode QR */
	eps = 1.e-12 ;
    {
	double x0[3], x1[3], s;
	double Ad[3*3]; /* Ad sera presque triangulaire sup�rieure */
	char sss[50];
	printf("\n\n-----------------------------------------------------------");
	methode_QR(As,3,Ad,k);
	printf("\nApr�s %d it�rations, As est similaire � Ad = " ,  k );
	afficher_m(Ad,3,14);

	/* V�rifions les d�terminants */
	printf("\n det(As) = %.14g" , determinant(As,3) );
	printf("\n det(Ad) = %.14g" , determinant(Ad,3) );


	/* Cherchons maintenant toutes les valeurs et vecteurs propres de As :
	   nous utilisons s= Ad(i,i) comme valeur approch�e de la valeur propre,
	   et utilisons la m�thode de la puissance avec translation:
	*/
	printf("\n ******************************************************");
	printf("\n Pour As, nous avons les r�sultats suivants:");
	for( i=1 ; i <= 3 ; i++){
		sprintf(sss,"%.2g",mat(Ad,3,i,i));	sscanf(sss,"%lf",&s);
		/* s = mat(Ad,3,i,i) ; */
		methode_de_la_puissance_shifted(As,3,s, x0, eps, &lambda);
		printf("\n lambda(%d) = %24.14g , et son vecteur associ� est x%d = " , i, lambda,i);
		afficher_v(x0,3,16);
		/* V�rifions */
		produit_matrice_vecteur(As,x0,x1,3);
		for(j=1;j<=3;j++)	vec(x1,j) /= lambda ;
		printf("\n\tVerification: lambda = %.14g\n As . x0 / lambda = "  , lambda);
		afficher_v(x1,3,16);

	}
	printf("\n ******************************************************\n\n");

    }

    {
	double x0[4], x1[4], s;
	double Bd[4*4]; /* Bd presque triangulaire sup�rieure */
	char sss[50];
	printf("\n\n-----------------------------------------------------------");
	methode_QR(Bs,4,Bd,k);
	printf("\nApr�s %d it�rations, Bs est similaire � Bd = " ,  k );
	afficher_m(Bd,4,14);

	/* V�rifions les d�terminants */
	printf("\n det(Bs) = %.14g" , determinant(Bs,4) );
	printf("\n det(Bd) = %.14g" , determinant(Bd,4) );


	/* Cherchons maintenant toutes les valeurs et vecteurs propres de Bs :
	   nous utilisons s= Bd(i,i) comme valeur approch�e de la valeur propre,
	   et utilisons la m�thode de la puissance avec translation:
	*/
	printf("\n ******************************************************");
	printf("\n Pour Bs, nous avons les r�sultats suivants:");
	for( i=1 ; i <= 4 ; i++){
		sprintf(sss,"%.2g",mat(Bd,4,i,i));	sscanf(sss,"%lf",&s);
		methode_de_la_puissance_shifted(Bs,4,s, x0, eps, &lambda);
		printf("\n lambda(%d) = %24.14g , et son vecteur associ� est x%d = " , i, lambda,i);
		afficher_v(x0,4,16);
		/* V�rifions */
		produit_matrice_vecteur(Bs,x0,x1,4);
		for(j=1;j<=4;j++)	vec(x1,j) /= lambda ;
		printf("\n\tVerification: lambda = %.14g\n Bs . x0 / lambda = "  , lambda);
		afficher_v(x1,4,16);
	}
	printf("\n ******************************************************\n\n");

    }

    {
	double x0[8], x1[8], s;
	double Cd[8*8]; /* Cd presque triangulaire sup�rieure */
	char sss[50];
	printf("\n\n-----------------------------------------------------------");
	methode_QR(Cs,8,Cd,k);
	printf("\nApr�s %d it�rations, Cs est similaire � Cd = " ,  k );
	afficher_m(Cd,8,6);

	/* V�rifions les d�terminants */
	printf("\n det(Cs) = %.14g" , determinant(Cs,8) );
	printf("\n det(Cd) = %.14g" , determinant(Cd,8) );


	/* Cherchons maintenant toutes les valeurs et vecteurs propres de Cs :
	   nous utilisons s= Cd(i,i) comme valeur approch�e de la valeur propre,
	   et utilisons la m�thode de la puissance avec translation:
	*/
	printf("\n ******************************************************");
	printf("\n Pour Cs, nous avons les r�sultats suivants:");
	for( i=1 ; i <= 8 ; i++){
		sprintf(sss,"%.2g",mat(Cd,8,i,i));	sscanf(sss,"%lf",&s);
		methode_de_la_puissance_shifted(Cs,8,s, x0, eps, &lambda);
		printf("\n lambda(%d) = %24.14g , et son vecteur associ� est x%d = " , i, lambda,i);
		afficher_v(x0,8,16);
		/* V�rifions */
		produit_matrice_vecteur(Cs,x0,x1,8);
		for(j=1;j<=8;j++)	vec(x1,j) /= lambda ;
		printf("\n\tVerification: lambda = %.14g\n Cs . x0 / lambda = "  , lambda);
		afficher_v(x1,8,16);
	}
	printf("\n ******************************************************\n\n");

    }

	return;

}

