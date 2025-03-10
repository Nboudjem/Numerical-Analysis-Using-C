/*
	Ce fichier contient des fonctions qui manipulent des matrices
	et des vecteurs (qui sont assimilés a des matrices a 1 seule colonne).
	Une matrice a_ij (n x m) sera représentée par un tableau A[n*m],
	où 	a[1][1] = A[0],		a[1][2] = A[1],	  ...,	 a[1][m] = A[m-1]
	 	a[2][1] = A[m],		a[2][2] = A[m+1],  ...,	 a[2][m] = A[2m-1]
		......................
		a[i][j] = A[ (i-1)*m + (j-1) ]
		......................

	Donc, un vecteur b_i qui a "n" composantes, sera assimilé a une
	matrice (n x 1)  et sera représenté par la matrice B[n*1],
	avec 
		b[i] = B[ (i-1)*1 + (1-1) ] == B[ i-1 ]

   Nous définirons nos matrices comme ceci (ci dessous, une matrice A, 3x3):
	double A[] = 
		{
			2,  3, -1 ,
			1, -2,  5 ,
			6, -3, 15
		} ;

   Pour les vecteurs, nous definirons de même (ci-dessous, un vecteur à 
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
#include <stdlib.h>

#define mat(A,m,i,j) 	A[ (i)*(m) - (m) + (j)-1 ] /* element A_ij , avec "m" colonnes */
#define vec(b,i) 	b[ (i)-1 ]		/* composante b_i */
/* Avec les définitions ci-dessus, on a évidemment: vec(B,i) == mat(B,1,i,1) */



void afficher_matrice( double A[], int n){
	/*
	Cette fonction affiche la matrice A (CARRÉE) sous une forme qui peut
	être ré-utilisée dans un programme C.
	les éléments sont A_ij avec i,j=1, ..., n 
	*/
	int i,j;
	printf("\t{"); /* cette chaine de caractères sera affichée telle quelle */
	for (i=1; i<n ; i++){/* on affiche les (n-1) premières lignes */
		printf("\n");
		for (j=1; j<=n ; j++)	printf("   %10.4g,", mat(A,n,i,j) );
	}

	/* pour i=n (dernière ligne), on affiche différemment */
	printf("\n");
	for (j=1; j<n ; j++)	printf("   %10.4g,", mat(A,n,i,j) );
	printf("   %10.4g\n\t};\n", mat(A,n,i,j) ); /* pour j=n */
}

void afficher_m( double A[], int n, int digits){
	/*
	Cette fonction affiche la matrice A (CARRÉE) sous une forme qui peut
	être ré-utilisée dans un programme C.
	les éléments sont A_ij avec i,j=1, ..., n 
	*/
	int i,j;
	char sss[50];
	sprintf(sss,"%%%d.%dg," , digits+8, digits);

	printf("\t{"); /* cette chaine de caractères sera affichée telle quelle */
	for (i=1; i<n ; i++){/* on affiche les (n-1) premières lignes */
		printf("\n");
		for (j=1; j<=n ; j++)	printf(sss, mat(A,n,i,j) );
	}

	/* pour i=n (dernière ligne), on affiche différemment */
	printf("\n");
	for (j=1; j<n ; j++)	printf(sss, mat(A,n,i,j) );
	sprintf(sss,"%%%d.%dg\n\t};\n" , digits+8, digits);
	printf(sss, mat(A,n,i,j) ); /* pour j=n */
}


void afficher_vecteur( double x[], int n){
	/* Cette fonction affiche le vecteur x qui a n éléments.	*/
	int i,j;
	printf("\t{"); /* cette chaine de caractères sera affichée telle quelle */
	for (i=1; i<n ; i++){ /* on affiche les (n-1) premières composantes */
		printf("\n   %10.4g,", vec(x,i) );
	}

	/* pour i=n (dernière composante), on affiche différemment */
	printf("\n   %10.4g\n\t};\n", vec(x,i) );
}
void afficher_v( double x[], int n, int digits){
	/* Cette fonction affiche le vecteur x qui a n éléments, avec digits chiffres significatifs */
	int i,j;
	char sss[50];
	sprintf(sss,"\n%%%d.%dg," , digits+8, digits);
	printf("\t{"); /* cette chaine de caractères sera affichée telle quelle */
	for (i=1; i<n ; i++){ /* on affiche les (n-1) premières composantes */
		printf(sss, vec(x,i) );
	}

	/* pour i=n (dernière composante), on affiche différemment */
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
	/* La matrice M est carrée,  (n x n), les vecteurs x et y ont n éléments.
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
	/* Les matrice M et N sont carrées,  (n x n): 
	   On calcule la matrice R égale au produit M.N
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
	/* Les matrice M et N sont carrées,  (n x n): 
	   On calcule la matrice R égale à la somme M+N
	*/
   int i,j,k;

   for (i=1; i<=n ; i++)
	for (j=1; j<=n ; j++)	mat(R,n,i,j) = mat(M,n,i,j) + mat(N,n,i,j) ;
return;
}

void produit_matrice_scalaire(double M[], double a, double N[], int n ){
	/* Les matrice M et N sont carrées,  (n x n): 
	   On calcule N = a M , "a" est un réel double
	*/
   int i,j,k;

   for (i=1; i<=n ; i++)
	for (j=1; j<=n ; j++)	mat(N,n,i,j) = a*mat(M,n,i,j) ;
return;
}


void gauss_LU(double A[], double L[], double U[], int n){
	/* Les matrice A, L et U sont carrées,  (n x n): 
	   Cette fonction transforme la matrice A en un produit L.U,
	   U est triangulaire supérieure,
	   L est triangulaire inférieure.
	*/
   int i,j,k;
   
   /* On initialise d'abord la matrice L à la matrice Identité */
   for (i=1; i<=n ; i++){
   	for (j=1; j<=n ; j++){
   		if (i==j) 	mat(L,n,i,j) = 1. ;
   		else		mat(L,n,i,j) = 0. ;
   	}
   }

   /* On copy A dans U, et on modifie U au fur et à mesure. A n'est pas modifié!!!! */
   copy_matrice(A,U,n);


   /* passage de la matrice U^{k} à U^{k+1} , pour k=1,...,(n-1) (sans interchange de lignes)*/
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
	double z , L[n*n], U[n*n]; /* On décompose A = L.U */
	int i;

	gauss_LU(A,L,U,n);
	z = 1. ; /* determinant de A initialisé à 1.   Remarque: det(L)=1. */ 
	for (i=1; i<=n ; i++) 	z *= mat(U,n,i,i); /* on multiplie les éléments diagonaux */
	return z;	
}

void resoudre_avec_gauss(double A[], double b[], double x[], int n){
	/*
	La matrice A est carrée,  (n x n): b est donné. on résoud A . x = b avec la méthode
	A=L.U, puis on résoud:
	ICI, nous n'interchangeons aucune ligne .
	*/
	double L[n*n], U[n*n], y[n]; /* On décompose A = L.U */
	int i,j;
	gauss_LU(A,L,U,n);

	/* On résoud L.y = b */
	for (i=1; i<=n ; i++) {
		vec(y,i) = vec(b,i) ;
		for (j=1; j<i ; j++) 		vec(y,i) -= mat(L,n,i,j) * vec(y,j) ;
		/* Ici on ne divise pas par L_ii, car L_ii = 1 pour la décomposition Doolittle */
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
	La matrice A est carrée,  (n x n). Nous supposons que son déterminant est différent de zéro.
	La matrice B a n lignes et m colonnes, et contient exactement "m" vecteurs
	b_k, et l'on veut résoudre A . x = b_k , pour k=1,..., k, en utilisant
	la méthode de Gauss (avec peut-être interchange de lignes si nécessaire).
	Les résultats sont mis dans la matrice B, et A n'est pas modifié!!!!
	*/
	double U[n*n], y[n]; /* On copiera A dans U pour ne pas modifier A! */
	int i,j,k, o[n];

	/* On copy A dans U, et on modifie U au fur et à mesure. A n'est pas modifié!!!! */
	copy_matrice(A,U,n);

	/* On initialise le vecteur d'ordre o[] */
	/* vec(o,i) représentera la ligne "i", même avec interchange de lignes */
	for (i=1; i<=n ; i++) vec(o,i) = i ;  


	/* 1ère étape: TRIANGULARISATION de "U". 	"U" deviendra triangulaire supérieure */
	for(k=1; k<n ; k++) {/* passage de la matrice U^{k} à U^{k+1} , pour k=1,...,(n-1)  */
		double L_ik ;
		for (i=k+1; i<=n ; i++){
		    if ( mat(U,n,vec(o,k),k) == 0. ){ /* Le pivot est nul !!!! */
			/* Cherchons un pivot différent de zéro en dessous de la ligne vec(o,k) */
			/* Si cela n'est pas possible, cela veut dire que det(U)=det(A)=0. , et on sort alors! */
			j = k+1;
			while ( mat(U,n,vec(o,j),k) == 0. && j <= n ) j++;
			if ( j > n ) {/* cela veut dire que l'on n'a pas trouvé de pivot != 0. ON SORT! */
				printf("\nCette matrice a un déterminant nul. La méthode gauss ne peut\n"
					"s'appliquer à cette matrice.\nVoici la matrice obtenue à ce stade:\n");
				afficher_matrice(U,n);
				return;
			}
			else {/* la ligne vec(o,j) a un pivot différent de zéro. On interchange les lignes
				vec(o,k) et vec(o,j) */
				double tmp=vec(o,k);
				vec(o,k) = vec(o,j);
				vec(o,j) = tmp;
			}
		    }
		    /* On continue avec la méthode de Gauss pour triangulariser U */
		    L_ik = mat(U,n,vec(o,i),k) / mat(U,n,vec(o,k),k) ;
		    mat(U,n,vec(o,i),k) = 0. ; /* cet élément doit s'annuler. On ne le calcule donc pas. */
		    for (j=k+1; j<=n ; j++)	mat(U,n,vec(o,i),j) -= L_ik * mat(U,n,vec(o,k),j) ;
		    for (j=1;   j<=m ; j++)	mat(B,m,vec(o,i),j) -= L_ik * mat(B,m,vec(o,k),j) ;
		}
	}

	/* 2ème étape: On résoud U . x = b_k, pour chaque colonne k=1,...,m */
	/* La solution pour U.x = b est : 
		pour i=n,...,1  (i est le numéro de la ligne)
			x_i = 1/U_ii ( b_i - sum_{j=i+1}^n   U_ij   x_j  )
	Les résultats seront mis dans la matrice B, avec la correspondance suivante:
		b_i  initialement égal à  B( vec(o,i) , k )
		le résultat x_i  est ensuite mis dans B( vec(o,i) , k )
	Ce qui revient à faire le remplacement 
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
	La matrice A est carrée,  (n x n). On la décompose comme P.A = A' = L.U
	La matrice de permutation est représentée par le vecteur d'ordre o[].
	La fonction gauss_PLU  retourne  
		 0 si le det(A)=0, sinon 
		 1 si det(A) = det(U), sinon 
		-1 si det(A) = - det(U)
	Donc, le système AX=B s'écrit A'.X = B' (=P.B) où A' et B' sont obtenues
	de A et B par interchanges de lignes. Et A' s'écrit L.U == A'
	*/
	int i,j,k,det=1;

	/* On initialise le vecteur d'ordre o[] */
	/* vec(o,i) représentera la ligne "i", même avec interchange de lignes */
	for (i=1; i<=n ; i++) vec(o,i) = i ;  

	/* On initialise L à l'identité. */
	for (i=1; i<=n ; i++)
		for (j=1; j<=n ; j++) if (i==j) mat(L,n,i,j)=1. ; else mat(L,n,i,j)=0. ;

	/* On copie A dans U */
	for (i=0; i<n*n ; i++) U[i] = A[i];

	for(k=1; k<n ; k++) {/* passage de la matrice U^{k} à U^{k+1} , pour k=1,...,(n-1)  */
		double L_ik ;
		for (i=k+1; i<=n ; i++){
		    if ( mat(U,n,vec(o,k),k) == 0. ){ /* Le pivot est nul !!!! */
			/* Cherchons un pivot différent de zéro en dessous de la ligne vec(o,k) */
			/* Si cela n'est pas possible, cela veut dire que det(U)=det(A)=0. */
			j = k+1;
			while ( mat(U,n,vec(o,j),k) == 0. && j <= n ) j++;
			if ( j > n ) {/* cela veut dire que l'on n'a pas trouvé de pivot != 0. donc det==0 */
				mat(L,n,i,k) = 0. ;
				det = 0. ;
				break; /* on passe au k suivant */
			}
			else {/* la ligne vec(o,j) a un pivot différent de zéro. On interchange les lignes
				vec(o,k) et vec(o,j) */
				double tmp=vec(o,k);
				vec(o,k) = vec(o,j);
				vec(o,j) = tmp;
				det = - det ;
			}
		    }
		    /* On continue avec la méthode de Gauss pour triangulariser U */
		    L_ik = mat(U,n,vec(o,i),k) / mat(U,n,vec(o,k),k) ;
		    mat(L,n,i,k) = L_ik ;

		    mat(U,n,vec(o,i),k) = 0. ; /* cet élément doit s'annuler. On ne le calcule donc pas. */
		    for (j=k+1; j<=n ; j++)	mat(U,n,vec(o,i),j) -= L_ik * mat(U,n,vec(o,k),j) ;
		}
	}
	/*
	A ce niveau, c'est la matrice V_ij === U_{ vec(o,i) , j } qui
	   est triangulaire supérieure.
	Pour résoudre A.X=B, (X et B matrices) nous faisons ceci:
	1. Nous transformons d'abord B en B' par la matrice L_ik, puis on résoud U.X=B' .
	*/
	/* Les équations seront: B en B'. B' est remis dans B */
/*
	for(k=1; k<n ; k++) {
		for (i=k+1; i<=n ; i++){
		    L_ik = mat(L,n,i,k) ;
		    for (j=1;   j<=m ; j++)	mat(B,m,vec(o,i),j) -= mat(L,n,i,k) * mat(B,m,vec(o,k),j) ;
		}
	}
*/
	/* Les équations seront: U.X=B' . Le résultat pour X dans B' */
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
	/* On calcule ici $\alpha A B + \beta C$ et on met le résultat dans C.
	   C: mxn
	   A: mxk  ,  B: kxn
	   alpha et beta réels
	Les équations sont:
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
	/* A est nxn ; B et X sont nxm ; On veut résoudre A.X = B avec la méthode PLU,
	puis en améliorant la précision de maniére itérative. Voir plus bas. */
	double det, L[n*n], U[n*n], DB[n*m], L_ik;
	int i,j,k,o[n],iii;

	det = gauss_PLU(A,n,o,L,U); /* On décompose d'abord A */
	if (det==0) {
		printf("\n\n det(A) = 0. Exit \n\n");
		return; /* Ce n'est pas la peine de résoudre car   det(A) = 0 */
	}

	/*
	Pour résoudre A.X=B, (X et B matrices) nous faisons ceci:
	1. Nous transformons d'abord B en B' par la matrice L_ik, puis on résoud U.X=B' .
	*/

	/* Nous copions d'abord B dans X */
	for (i=0; i<n*m ; i++)	X[i] = B[i] ;

	/* Les équations seront: B en B'. B' est remis dans B (c'est dire dans X=B) */
	for(k=1; k<n ; k++) {
		for (i=k+1; i<=n ; i++){
		    L_ik = mat(L,n,i,k) ;
		    for (j=1;   j<=m ; j++)	mat(X,m,vec(o,i),j) -= L_ik * mat(X,m,vec(o,k),j) ;
		}
	}
	/* Les équations seront: U.X=B' . Le résultat pour X dans B' (donc dans X) */
	for(k=1; k<=m ; k++) {
	    for (i=n; i>=1 ; i--) {
		for (j=i+1; j<=n ; j++) 
			mat(X,m,vec(o,i),k) -= mat(U,n,vec(o,i),j) * mat(X,m,vec(o,j),k) ;
		mat(X,m,vec(o,i),k) /= mat(U,n,vec(o,i),i) ;
	    }	
	}
	
	/* Jusqu'ici, nous avons une solution approchée X que nous appellerons Xo ici
	DB = -A.X + B
	 */
/* ************************************************************************************/
/* Voici comment on peut itérer pour rendre plus précis A.X = B :
	On pose X = Xo + DX (toutes des matrices)
	A ( Xo + DX ) = B  ==> A . DX = B - A.Xo = DB à résoudre.
	La solution sera: X <--  Xo + DX. Cette méthode pourra être itérée comme ceci:
	A et B sont donnés.
		On copy B dans X et on résoud A.X=B; la solution est dans X maintenant.
	     1.	
		On calcule  DB = B - A.X,
		puis on résoud A.DX = DB et DX esr mis dans DB ;
		Ensuite X <- X + DX, et on revient à l'étape 1.
*/
/* ************************************************************************************/

    for (iii=1 ; iii<=ITERATIONS ; iii++){
	/* Nous copions B dans DB */
	for (i=0; i<n*m ; i++)	DB[i] = B[i] ;

	produit_aAB_plus_bC(n, m, n, -1., A, X, 1., DB) ; /* DB <--   -A.X + DB */
	/* On résoud maintenant A.DX = DB . Le resultat DX est ensuite mis dans DB */
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
	/* A est nxn ; B et X sont nxm ; On veut résoudre A.X = B avec la méthode PLU,
	puis en améliorant la précision de maniére itérative. Voir plus bas. */
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


	det = gauss_PLU(A,n,o,L,U); /* On décompose d'abord A */
	if (det==0) {
		printf("\n\n det(A) = 0. Exit \n\n");
		return; /* Ce n'est pas la peine de résoudre car   det(A) = 0 */
	}

	/*
	Pour résoudre A.X=B, (X et B matrices) nous faisons ceci:
	1. Nous transformons d'abord B en B' par la matrice L_ik, puis on résoud U.X=B' .
	*/

	/* Nous copions d'abord B dans X */
	for (i=0; i<n*m ; i++)	X[i] = B[i] ;

	/* Les équations seront: B en B'. B' est remis dans B (c'est dire dans X=B) */
	for(k=1; k<n ; k++) {
		for (i=k+1; i<=n ; i++){
		    L_ik = mat(L,n,i,k) ;
		    for (j=1;   j<=m ; j++)	mat(X,m,vec(o,i),j) -= L_ik * mat(X,m,vec(o,k),j) ;
		}
	}
	/* Les équations seront: U.X=B' . Le résultat pour X dans B' (donc dans X) */
	for(k=1; k<=m ; k++) {
	    for (i=n; i>=1 ; i--) {
		for (j=i+1; j<=n ; j++) 
			mat(X,m,vec(o,i),k) -= mat(U,n,vec(o,i),j) * mat(X,m,vec(o,j),k) ;
		mat(X,m,vec(o,i),k) /= mat(U,n,vec(o,i),i) ;
	    }	
	}
	
	/* Jusqu'ici, nous avons une solution approchée X que nous appellerons Xo ici
	DB = -A.X + B
	 */
/* ************************************************************************************/
/* Voici comment on peut itérer pour rendre plus précis A.X = B :
	On pose X = Xo + DX (toutes des matrices)
	A ( Xo + DX ) = B  ==> A . DX = B - A.Xo = DB à résoudre.
	La solution sera: X <--  Xo + DX. Cette méthode pourra être itérée comme ceci:
	A et B sont donnés.
		On copy B dans X et on résoud A.X=B; la solution est dans X maintenant.
	     1.	
		On calcule  DB = B - A.X,
		puis on résoud A.DX = DB et DX esr mis dans DB ;
		Ensuite X <- X + DX, et on revient à l'étape 1.
*/
/* ************************************************************************************/

    for (iii=1 ; iii<=ITERATIONS ; iii++){
	/* Nous copions B dans DB */
	for (i=0; i<n*m ; i++)	DB[i] = B[i] ;

	produit_aAB_plus_bC(n, m, n, -1., A, X, 1., DB) ; /* DB <--   -A.X + DB */
	/* On résoud maintenant A.DX = DB . Le resultat DX est ensuite mis dans DB */
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
