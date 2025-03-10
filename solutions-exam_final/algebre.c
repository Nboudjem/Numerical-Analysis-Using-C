/*
	Ce fichier contient des fonctions qui manipulent des matrices
	et des vecteurs (qui sont assimil�s a des matrices a 1 seule colonne).
	Une matrice a_ij (n x m) sera repr�sent�e par un tableau A[n*m],
	o� 	a[1][1] = A[0],		a[1][2] = A[1],	  ...,	 a[1][m] = A[m-1]
	 	a[2][1] = A[m],		a[2][2] = A[m+1],  ...,	 a[2][m] = A[2m-1]
		......................
		a[i][j] = A[ (i-1)*m + (j-1) ]
		......................

	Donc, un vecteur b_i qui a "n" composantes, sera assimil� a une
	matrice (n x 1)  et sera repr�sent� par la matrice B[n*1],
	avec 
		b[i] = B[ (i-1)*1 + (1-1) ] == B[ i-1 ]

   Nous d�finirons nos matrices comme ceci (ci dessous, une matrice A, 3x3):
	double A[] = 
		{
			2,  3, -1 ,
			1, -2,  5 ,
			6, -3, 15
		} ;

   Pour les vecteurs, nous definirons de m�me (ci-dessous, un vecteur � 
	3 composantes):
	double B[] = 
		{
			-1 ,
			 5 ,
			15
		} ;

*/

#include <stdio.h>
#include <math.h>

#define mat(A,m,i,j) 	A[ (i)*(m) - (m) + (j)-1 ] /* element A_ij , avec "m" colonnes */
#define vec(b,i) 	b[ (i)-1 ]		/* composante b_i */
/* Avec les d�finitions ci-dessus, on a �videmment: vec(B,i) == mat(B,1,i,1) */



void afficher_matrice( double A[], int n){
	/*
	Cette fonction affiche la matrice A (CARR�E) sous une forme qui peut
	�tre r�-utilis�e dans un programme C.
	les �l�ments sont A_ij avec i,j=1, ..., n 
	*/
	int i,j;
	printf("\t{"); /* cette chaine de caract�res sera affich�e telle quelle */
	for (i=1; i<n ; i++){/* on affiche les (n-1) premi�res lignes */
		printf("\n");
		for (j=1; j<=n ; j++)	printf("   %10.4g,", mat(A,n,i,j) );
	}

	/* pour i=n (derni�re ligne), on affiche diff�remment */
	printf("\n");
	for (j=1; j<n ; j++)	printf("   %10.4g,", mat(A,n,i,j) );
	printf("   %10.4g\n\t};\n", mat(A,n,i,j) ); /* pour j=n */
}

void afficher_m( double A[], int n, int digits){
	/*
	Cette fonction affiche la matrice A (CARR�E) sous une forme qui peut
	�tre r�-utilis�e dans un programme C.
	les �l�ments sont A_ij avec i,j=1, ..., n 
	*/
	int i,j;
	char sss[50];
	sprintf(sss,"%%%d.%dg," , digits+8, digits);

	printf("\t{"); /* cette chaine de caract�res sera affich�e telle quelle */
	for (i=1; i<n ; i++){/* on affiche les (n-1) premi�res lignes */
		printf("\n");
		for (j=1; j<=n ; j++)	printf(sss, mat(A,n,i,j) );
	}

	/* pour i=n (derni�re ligne), on affiche diff�remment */
	printf("\n");
	for (j=1; j<n ; j++)	printf(sss, mat(A,n,i,j) );
	sprintf(sss,"%%%d.%dg\n\t};\n" , digits+8, digits);
	printf(sss, mat(A,n,i,j) ); /* pour j=n */
}


void afficher_vecteur( double x[], int n){
	/* Cette fonction affiche le vecteur x qui a n �l�ments.	*/
	int i,j;
	printf("\t{"); /* cette chaine de caract�res sera affich�e telle quelle */
	for (i=1; i<n ; i++){ /* on affiche les (n-1) premi�res composantes */
		printf("\n   %10.4g,", vec(x,i) );
	}

	/* pour i=n (derni�re composante), on affiche diff�remment */
	printf("\n   %10.4g\n\t};\n", vec(x,i) );
}
void afficher_v( double x[], int n, int digits){
	/* Cette fonction affiche le vecteur x qui a n �l�ments, avec digits chiffres significatifs */
	int i,j;
	char sss[50];
	sprintf(sss,"\n%%%d.%dg," , digits+8, digits);
	printf("\t{"); /* cette chaine de caract�res sera affich�e telle quelle */
	for (i=1; i<n ; i++){ /* on affiche les (n-1) premi�res composantes */
		printf(sss, vec(x,i) );
	}

	/* pour i=n (derni�re composante), on affiche diff�remment */
	sprintf(sss,"\n%%%d.%dg\n\t};\n" , digits+8, digits);
	printf(sss, vec(x,i) );
}

void copy_matrice( double A[], double B[], int n){
	/* copy la matrice A et met le resultat dans la matrice B */
	int i,j;
	for(i=1;i<=n;i++){
		for(j=1;j<=n;j++) 	mat(B,n,i,j) = mat(A,n,i,j);
	}
}


void copy_vecteur( double x[], double y[], int n){
	/* copy le vecteur x et met le resultat dans le vecteur y */
	int i;
	for(i=1;i<=n;i++)		vec(y,i) = vec(x,i);
}


void matrice_identite( double A[], int n){
	/* poser A = I (nxn) */
	int i,j;
	for(i=1;i<=n;i++){
		for(j=1;j<i;j++)	mat(A,n,i,j) = 0. ; 
		mat(A,n,i,i) = 1.; 
		for(j=i+1;j<=n;j++)	mat(A,n,i,j) = 0. ; 
	}
}


void transpose( double A[], double B[], int n){
	/* calcule B = transpose de A */
	int i,j;
	for(i=1;i<=n;i++){
		for(j=1;j<=n;j++)	mat(B,n,i,j) = mat(A,n,j,i); 
	}
}


void produit_matrice_vecteur(double M[], double x[], double y[], int n){
	/* La matrice M est carr�e,  (n x n), les vecteurs x et y ont n �l�ments.
	   On calcule le produit y = M.x
	*/
   int i,j;

   for (i=1; i<=n ; i++){
	vec(y,i)  = 0. ;
	for(j=1;j<=n;j++) vec(y,i) += mat(M,n,i,j) * vec(x,j) ;
   }
return;
}


void produit_matrice_matrice(double M[], double N[], double R[], int n){
	/* Les matrice M et N sont carr�es,  (n x n): 
	   On calcule la matrice R �gale au produit M.N
	*/
   int i,j,k;

   for (i=1; i<=n ; i++){
	for (j=1; j<=n ; j++){
		mat(R,n,i,j)  = 0. ;
		for(k=1;k<=n;k++) 	mat(R,n,i,j) += mat(M,n,i,k) * mat(N,n,k,j) ;
	}
   }
return;
}

void somme_matrice_matrice(double M[], double N[], double R[], int n){
	/* Les matrice M et N sont carr�es,  (n x n): 
	   On calcule la matrice R �gale � la somme M+N
	*/
   int i,j,k;

   for (i=1; i<=n ; i++)
	for (j=1; j<=n ; j++)	mat(R,n,i,j) = mat(M,n,i,j) + mat(N,n,i,j) ;
return;
}

void produit_matrice_scalaire(double M[], double a, double N[], int n ){
	/* Les matrice M et N sont carr�es,  (n x n): 
	   On calcule N = a M , "a" est un r�el double
	*/
   int i,j,k;

   for (i=1; i<=n ; i++)
	for (j=1; j<=n ; j++)	mat(N,n,i,j) = a*mat(M,n,i,j) ;
return;
}


void gauss_LU(double A[], double L[], double U[], int n){
	/* Les matrice A, L et U sont carr�es,  (n x n): 
	   Cette fonction transforme la matrice A en un produit L.U,
	   U est triangulaire sup�rieure,
	   L est triangulaire inf�rieure.
	*/
   int i,j,k;
   
   /* On initialise d'abord la matrice L � la matrice Identit� */
   for (i=1; i<=n ; i++){
   	for (j=1; j<=n ; j++){
   		if (i==j) 	mat(L,n,i,j) = 1. ;
   		else		mat(L,n,i,j) = 0. ;
   	}
   }

   /* On copy A dans U, et on modifie U au fur et � mesure. A n'est pas modifi�!!!! */
   copy_matrice(A,U,n);


   /* passage de la matrice U^{k} � U^{k+1} , pour k=1,...,(n-1) (sans interchange de lignes)*/
   for(k=1; k<n ; k++) {
	for (i=k+1; i<=n ; i++){
		mat(L,n,i,k) = mat(U,n,i,k) / mat(U,n,k,k) ;
		mat(U,n,i,k) = 0. ;
		for (j=k+1; j<=n ; j++)		mat(U,n,i,j) -= mat(L,n,i,k) * mat(U,n,k,j) ;
	}
   }
return;
}



double determinant(double A[], int n) {
	double z , L[n*n], U[n*n]; /* On d�compose A = L.U */
	int i;

	gauss_LU(A,L,U,n);
	z = 1. ; /* determinant de A initialis� � 1.   Remarque: det(L)=1. */ 
	for (i=1; i<=n ; i++) 	z *= mat(U,n,i,i); /* on multiplie les �l�ments diagonaux */
	return z;	
}

void resoudre_avec_gauss(double A[], double b[], double x[], int n){
	/*
	La matrice A est carr�e,  (n x n): b est donn�. on r�soud A . x = b avec la m�thode
	A=L.U, puis on r�soud:
	ICI, nous n'interchangeons aucune ligne .
	*/
	double L[n*n], U[n*n], y[n]; /* On d�compose A = L.U */
	int i,j;
	gauss_LU(A,L,U,n);

	/* On r�soud L.y = b */
	for (i=1; i<=n ; i++) {
		vec(y,i) = vec(b,i) ;
		for (j=1; j<i ; j++) 		vec(y,i) -= mat(L,n,i,j) * vec(y,j) ;
		/* Ici on ne divise pas par L_ii, car L_ii = 1 pour la d�composition Doolittle */
	}		

	/* ensuite U.x = y */
	for (i=n; i>=1 ; i--) {
		vec(x,i) = vec(y,i) ;
		for (j=i+1; j<=n ; j++) 	vec(x,i) -= mat(U,n,i,j) * vec(x,j) ;
		vec(x,i) /= mat(U,n,i,i) ;
	}		
}


void gauss_multi(double A[], int n, double B[], int m){
	/*
	La matrice A est carr�e,  (n x n). Nous supposons que son d�terminant est diff�rent de z�ro.
	La matrice B a n lignes et m colonnes, et contient exactement "m" vecteurs
	b_k, et l'on veut r�soudre A . x = b_k , pour k=1,..., k, en utilisant
	la m�thode de Gauss (avec peut-�tre interchange de lignes si n�cessaire).
	Les r�sultats sont mis dans la matrice B, et A n'est pas modifi�!!!!
	*/
	double U[n*n], y[n]; /* On copiera A dans U pour ne pas modifier A! */
	int i,j,k, o[n];

	/* On copy A dans U, et on modifie U au fur et � mesure. A n'est pas modifi�!!!! */
	copy_matrice(A,U,n);

	/* On initialise le vecteur d'ordre o[] */
	/* vec(o,i) repr�sentera la ligne "i", m�me avec interchange de lignes */
	for (i=1; i<=n ; i++) vec(o,i) = i ;  


	/* 1�re �tape: TRIANGULARISATION de "U". 	"U" deviendra triangulaire sup�rieure */
	for(k=1; k<n ; k++) {/* passage de la matrice U^{k} � U^{k+1} , pour k=1,...,(n-1)  */
		double L_ik ;
		for (i=k+1; i<=n ; i++){
		    if ( mat(U,n,vec(o,k),k) == 0. ){ /* Le pivot est nul !!!! */
			/* Cherchons un pivot diff�rent de z�ro en dessous de la ligne vec(o,k) */
			/* Si cela n'est pas possible, cela veut dire que det(U)=det(A)=0. , et on sort alors! */
			j = k+1;
			while ( mat(U,n,vec(o,j),k) == 0. && j <= n ) j++;
			if ( j > n ) {/* cela veut dire que l'on n'a pas trouv� de pivot != 0. ON SORT! */
				printf("\nCette matrice a un d�terminant nul. La m�thode gauss ne peut\n"
					"s'appliquer � cette matrice.\nVoici la matrice obtenue � ce stade:\n");
				afficher_matrice(U,n);
				return;
			}
			else {/* la ligne vec(o,j) a un pivot diff�rent de z�ro. On interchange les lignes
				vec(o,k) et vec(o,j) */
				double tmp=vec(o,k);
				vec(o,k) = vec(o,j);
				vec(o,j) = tmp;
			}
		    }
		    /* On continue avec la m�thode de Gauss pour triangulariser U */
		    L_ik = mat(U,n,vec(o,i),k) / mat(U,n,vec(o,k),k) ;
		    mat(U,n,vec(o,i),k) = 0. ; /* cet �l�ment doit s'annuler. On ne le calcule donc pas. */
		    for (j=k+1; j<=n ; j++)	mat(U,n,vec(o,i),j) -= L_ik * mat(U,n,vec(o,k),j) ;
		    for (j=1;   j<=m ; j++)	mat(B,m,vec(o,i),j) -= L_ik * mat(B,m,vec(o,k),j) ;
		}
	}

	/* 2�me �tape: On r�soud U . x = b_k, pour chaque colonne k=1,...,m */
	/* La solution pour U.x = b est : 
		pour i=n,...,1  (i est le num�ro de la ligne)
			x_i = 1/U_ii ( b_i - sum_{j=i+1}^n   U_ij   x_j  )
	Les r�sultats seront mis dans la matrice B, avec la correspondance suivante:
		b_i  initialement �gal �  B( vec(o,i) , k )
		le r�sultat x_i  est ensuite mis dans B( vec(o,i) , k )
	Ce qui revient � faire le remplacement 
	B( vec(o,i) , k )  -->   1/U_ii * ( B(vec(o,i) , k ) - sum_{j=i+1}^n   U_ij * B( vec(o,j) , k ) )
	*/
	for(k=1; k<=m ; k++) {
	    for (i=n; i>=1 ; i--) {
		for (j=i+1; j<=n ; j++) 
			mat(B,m,vec(o,i),k) -= mat(U,n,vec(o,i),j) * mat(B,m,vec(o,j),k) ;
		mat(B,m,vec(o,i),k) /= mat(U,n,vec(o,i),i) ;
	    }	
	}
return;
}

int gauss_PLU(double A[], int n, int o[], double L[], double U[]){
	/*
	La matrice A est carr�e,  (n x n). On la d�compose comme P.A = A' = L.U
	La matrice de permutation est repr�sent�e par le vecteur d'ordre o[].
	La fonction gauss_PLU  retourne  
		 0 si le det(A)=0, sinon 
		 1 si det(A) = det(U), sinon 
		-1 si det(A) = - det(U)
	Donc, le syst�me AX=B s'�crit A'.X = B' (=P.B) o� A' et B' sont obtenues
	de A et B par interchanges de lignes. Et A' s'�crit L.U == A'
	*/
	int i,j,k,det=1;

	/* On initialise le vecteur d'ordre o[] */
	/* vec(o,i) repr�sentera la ligne "i", m�me avec interchange de lignes */
	for (i=1; i<=n ; i++) vec(o,i) = i ;  

	/* On initialise L � l'identit�. */
	for (i=1; i<=n ; i++)
		for (j=1; j<=n ; j++) if (i==j) mat(L,n,i,j)=1. ; else mat(L,n,i,j)=0. ;

	/* On copie A dans U */
	for (i=0; i<n*n ; i++) U[i] = A[i];

	for(k=1; k<n ; k++) {/* passage de la matrice U^{k} � U^{k+1} , pour k=1,...,(n-1)  */
		double L_ik ;
		for (i=k+1; i<=n ; i++){
		    if ( mat(U,n,vec(o,k),k) == 0. ){ /* Le pivot est nul !!!! */
			/* Cherchons un pivot diff�rent de z�ro en dessous de la ligne vec(o,k) */
			/* Si cela n'est pas possible, cela veut dire que det(U)=det(A)=0. */
			j = k+1;
			while ( mat(U,n,vec(o,j),k) == 0. && j <= n ) j++;
			if ( j > n ) {/* cela veut dire que l'on n'a pas trouv� de pivot != 0. donc det==0 */
				mat(L,n,i,k) = 0. ;
				det = 0. ;
				break; /* on passe au k suivant */
			}
			else {/* la ligne vec(o,j) a un pivot diff�rent de z�ro. On interchange les lignes
				vec(o,k) et vec(o,j) */
				double tmp=vec(o,k);
				vec(o,k) = vec(o,j);
				vec(o,j) = tmp;
				det = - det ;
			}
		    }
		    /* On continue avec la m�thode de Gauss pour triangulariser U */
		    L_ik = mat(U,n,vec(o,i),k) / mat(U,n,vec(o,k),k) ;
		    mat(L,n,i,k) = L_ik ;

		    mat(U,n,vec(o,i),k) = 0. ; /* cet �l�ment doit s'annuler. On ne le calcule donc pas. */
		    for (j=k+1; j<=n ; j++)	mat(U,n,vec(o,i),j) -= L_ik * mat(U,n,vec(o,k),j) ;
		}
	}
	/*
	A ce niveau, c'est la matrice V_ij === U_{ vec(o,i) , j } qui
	   est triangulaire sup�rieure.
	Pour r�soudre A.X=B, (X et B matrices) nous faisons ceci:
	1. Nous transformons d'abord B en B' par la matrice L_ik, puis on r�soud U.X=B' .
	*/
	/* Les �quations seront: B en B'. B' est remis dans B */
/*
	for(k=1; k<n ; k++) {
		for (i=k+1; i<=n ; i++){
		    L_ik = mat(L,n,i,k) ;
		    for (j=1;   j<=m ; j++)	mat(B,m,vec(o,i),j) -= mat(L,n,i,k) * mat(B,m,vec(o,k),j) ;
		}
	}
*/
/* Les �quations seront: U.X=B' . Le r�sultat pour X dans B' */
/*
	for(k=1; k<=m ; k++) {
	    for (i=n; i>=1 ; i--) {
		for (j=i+1; j<=n ; j++) 
			mat(B,m,vec(o,i),k) -= mat(U,n,vec(o,i),j) * mat(B,m,vec(o,j),k) ;
		mat(B,m,vec(o,i),k) /= mat(U,n,vec(o,i),i) ;
	    }	
	}
*/
return det;
}

void produit_aAB_plus_bC(int m, int n, int k, double alpha, double A[], double B[], double beta,double C[]){
	/* On calcule ici $\alpha A B + \beta C$ et on met le r�sultat dans C.
	   C: mxn
	   A: mxk  ,  B: kxn
	   alpha et beta r�els
	Les �quations sont:
		C_ij <-- beta C_ij + alpha sum_{l=1}^k A_il B_lj
	*/
   int i,j,l;
   double sum;
   for (i=1; i<=m ; i++){
	for (j=1; j<=n ; j++){
		if (beta == 0.) 	mat(C,n,i,j)=0. ;
		else if (beta == 1. )	;
		else if (beta == -1. )	mat(C,n,i,j) = - mat(C,n,i,j) ;
		else 			mat(C,n,i,j) *= beta ;
		if(alpha != 0.){
			sum=0.;
			for(l=1;l<=k;l++) 	sum += mat(A,k,i,l) * mat(B,n,l,j) ;
			if (alpha == -1. )	mat(C,n,i,j) -= sum ;
			else if (alpha == 1. )	mat(C,n,i,j) += sum ;
			else 			mat(C,n,i,j) += alpha*sum ;
		}
	}
   }
return;
}

void resoudre_avec_gauss_PLU_keep(double A[], int n, double B[], int m, double X[], int ITERATIONS){
	/* A est nxn ; B et X sont nxm ; On veut r�soudre A.X = B avec la m�thode PLU,
	puis en am�liorant la pr�cision de mani�re it�rative. Voir plus bas. */
	double det, L[n*n], U[n*n], DB[n*m], L_ik;
	int i,j,k,o[n],iii;

	det = gauss_PLU(A,n,o,L,U); /* On d�compose d'abord A */
	if (det==0) {
		printf("\n\n det(A) = 0. Exit \n\n");
		return; /* Ce n'est pas la peine de r�soudre car   det(A) = 0 */
	}

	/*
	Pour r�soudre A.X=B, (X et B matrices) nous faisons ceci:
	1. Nous transformons d'abord B en B' par la matrice L_ik, puis on r�soud U.X=B' .
	*/

	/* Nous copions d'abord B dans X */
	for (i=0; i<n*m ; i++)	X[i] = B[i] ;

	/* Les �quations seront: B en B'. B' est remis dans B (c'est dire dans X=B) */
	for(k=1; k<n ; k++) {
		for (i=k+1; i<=n ; i++){
		    L_ik = mat(L,n,i,k) ;
		    for (j=1;   j<=m ; j++)	mat(X,m,vec(o,i),j) -= L_ik * mat(X,m,vec(o,k),j) ;
		}
	}
	/* Les �quations seront: U.X=B' . Le r�sultat pour X dans B' (donc dans X) */
	for(k=1; k<=m ; k++) {
	    for (i=n; i>=1 ; i--) {
		for (j=i+1; j<=n ; j++) 
			mat(X,m,vec(o,i),k) -= mat(U,n,vec(o,i),j) * mat(X,m,vec(o,j),k) ;
		mat(X,m,vec(o,i),k) /= mat(U,n,vec(o,i),i) ;
	    }	
	}
	
	/* Jusqu'ici, nous avons une solution approch�e X que nous appellerons Xo ici
	DB = -A.X + B
	 */
/* ************************************************************************************/
/* Voici comment on peut it�rer pour rendre plus pr�cis A.X = B :
	On pose X = Xo + DX (toutes des matrices)
	A ( Xo + DX ) = B  ==> A . DX = B - A.Xo = DB � r�soudre.
	La solution sera: X <--  Xo + DX. Cette m�thode pourra �tre it�r�e comme ceci:
	A et B sont donn�s.
		On copy B dans X et on r�soud A.X=B; la solution est dans X maintenant.
	     1.	
		On calcule  DB = B - A.X,
		puis on r�soud A.DX = DB et DX esr mis dans DB ;
		Ensuite X <- X + DX, et on revient � l'�tape 1.
*/
/* ************************************************************************************/

    for (iii=1 ; iii<=ITERATIONS ; iii++){
	/* Nous copions B dans DB */
	for (i=0; i<n*m ; i++)	DB[i] = B[i] ;

	produit_aAB_plus_bC(n, m, n, -1., A, X, 1., DB) ; /* DB <--   -A.X + DB */
	/* On r�soud maintenant A.DX = DB . Le resultat DX est ensuite mis dans DB */
	for(k=1; k<n ; k++) {
		for (i=k+1; i<=n ; i++){
		    L_ik = mat(L,n,i,k) ;
		    for (j=1;   j<=m ; j++)	mat(DB,m,vec(o,i),j) -= L_ik * mat(DB,m,vec(o,k),j) ;
		}
	}
	for(k=1; k<=m ; k++) {
	    for (i=n; i>=1 ; i--) {
		for (j=i+1; j<=n ; j++) 
			mat(DB,m,vec(o,i),k) -= mat(U,n,vec(o,i),j) * mat(DB,m,vec(o,j),k) ;
		mat(DB,m,vec(o,i),k) /= mat(U,n,vec(o,i),i) ;
	    }	
	}
	/* Etape finale: on fait X <-- X + DB */
	for (i=0; i<n*m ; i++)	X[i] += DB[i] ;
    }
return;
}

void resoudre_avec_gauss_PLU(double A[], int n, double B[], int m, double X[], int ITERATIONS){
	/* A est nxn ; B et X sont nxm ; On veut r�soudre A.X = B avec la m�thode PLU,
	puis en am�liorant la pr�cision de mani�re it�rative. Voir plus bas. */
/*
	double det, L[n*n], U[n*n], DB[n*m], L_ik;
	int i,j,k,o[n],iii;
*/
	double *L, *U, *DB, det, L_ik;
	int *o,i,j,k,iii;

	L = 		calloc( sizeof(double) , n*n );
	U = 		calloc( sizeof(double) , n*n );
	DB = 		calloc( sizeof(double) , n*m );
	o = 		calloc( sizeof(int) , n );


	det = gauss_PLU(A,n,o,L,U); /* On d�compose d'abord A */
	if (det==0) {
		printf("\n\n det(A) = 0. Exit \n\n");
		return; /* Ce n'est pas la peine de r�soudre car   det(A) = 0 */
	}

	/*
	Pour r�soudre A.X=B, (X et B matrices) nous faisons ceci:
	1. Nous transformons d'abord B en B' par la matrice L_ik, puis on r�soud U.X=B' .
	*/

	/* Nous copions d'abord B dans X */
	for (i=0; i<n*m ; i++)	X[i] = B[i] ;

	/* Les �quations seront: B en B'. B' est remis dans B (c'est dire dans X=B) */
	for(k=1; k<n ; k++) {
		for (i=k+1; i<=n ; i++){
		    L_ik = mat(L,n,i,k) ;
		    for (j=1;   j<=m ; j++)	mat(X,m,vec(o,i),j) -= L_ik * mat(X,m,vec(o,k),j) ;
		}
	}
	/* Les �quations seront: U.X=B' . Le r�sultat pour X dans B' (donc dans X) */
	for(k=1; k<=m ; k++) {
	    for (i=n; i>=1 ; i--) {
		for (j=i+1; j<=n ; j++) 
			mat(X,m,vec(o,i),k) -= mat(U,n,vec(o,i),j) * mat(X,m,vec(o,j),k) ;
		mat(X,m,vec(o,i),k) /= mat(U,n,vec(o,i),i) ;
	    }	
	}
	
	/* Jusqu'ici, nous avons une solution approch�e X que nous appellerons Xo ici
	DB = -A.X + B
	 */
/* ************************************************************************************/
/* Voici comment on peut it�rer pour rendre plus pr�cis A.X = B :
	On pose X = Xo + DX (toutes des matrices)
	A ( Xo + DX ) = B  ==> A . DX = B - A.Xo = DB � r�soudre.
	La solution sera: X <--  Xo + DX. Cette m�thode pourra �tre it�r�e comme ceci:
	A et B sont donn�s.
		On copy B dans X et on r�soud A.X=B; la solution est dans X maintenant.
	     1.	
		On calcule  DB = B - A.X,
		puis on r�soud A.DX = DB et DX esr mis dans DB ;
		Ensuite X <- X + DX, et on revient � l'�tape 1.
*/
/* ************************************************************************************/

    for (iii=1 ; iii<=ITERATIONS ; iii++){
	/* Nous copions B dans DB */
	for (i=0; i<n*m ; i++)	DB[i] = B[i] ;

	produit_aAB_plus_bC(n, m, n, -1., A, X, 1., DB) ; /* DB <--   -A.X + DB */
	/* On r�soud maintenant A.DX = DB . Le resultat DX est ensuite mis dans DB */
	for(k=1; k<n ; k++) {
		for (i=k+1; i<=n ; i++){
		    L_ik = mat(L,n,i,k) ;
		    for (j=1;   j<=m ; j++)	mat(DB,m,vec(o,i),j) -= L_ik * mat(DB,m,vec(o,k),j) ;
		}
	}
	for(k=1; k<=m ; k++) {
	    for (i=n; i>=1 ; i--) {
		for (j=i+1; j<=n ; j++) 
			mat(DB,m,vec(o,i),k) -= mat(U,n,vec(o,i),j) * mat(DB,m,vec(o,j),k) ;
		mat(DB,m,vec(o,i),k) /= mat(U,n,vec(o,i),i) ;
	    }	
	}
	/* Etape finale: on fait X <-- X + DB */
	for (i=0; i<n*m ; i++)	X[i] += DB[i] ;
    }
	free(L);
	free(U);
	free(DB);
	free(o);
return;
}



void inverse_matrice_simple( double A[], double A_inverse[], int n){
	/* calcule l'inverse de la matrice A, et met le resultat dans A_inverse */
	matrice_identite(A_inverse,n);
	gauss_multi(A,n,A_inverse,n);
}

void inverse_matrice( double A[], double A_inverse[], int n, int ITERATIONS){
	/* calcule l'inverse de la matrice A, et met le resultat dans A_inverse */
/*	double B[n*n]; */
	double *B;
	B = 		calloc( sizeof(double) , n*n );

	matrice_identite(B,n);
	resoudre_avec_gauss_PLU(A,n,B,n,A_inverse,ITERATIONS);
	free(B);
	return;
}
