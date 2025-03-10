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
	/* o[i] représentera la ligne "i", même avec interchange de lignes */
	for (i=1; i<=n ; i++) o[i] = i ;  


	/* 1ère étape: TRIANGULARISATION de "U". 	"U" devient triangulaire supérieure */
	for(k=1; k<n ; k++) {/* passage de la matrice U^{k} à U^{k+1} , pour k=1,...,(n-1)  */
		double L_ik ;
		for (i=k+1; i<=n ; i++){
		    if ( mat(U,n,o[k],k) == 0. ){ /* Le pivot est nul !!!! */
			/* Cherchons un pivot différent de zéro en dessous de la ligne o[k] */
			/* Si cela n'est pas possible, cela veut dire que det(U)=det(A)=0. , et on sort alors! */
			j = k+1;
			while ( mat(U,n,o[j],k) == 0. && j <= n ) j++;
			if ( j > n ) {/* cela veut dire que l'on n'a pas trouvé de pivot != 0. ON SORT! */
				printf("\nCette matrice a un déterminant nul. La méthode gauss ne peut\n"
					"s'appliquer à cette matrice.\nVoici la matrice obtenue à ce stade:\n");
				afficher_matrice(U,n);
				return;
			}
			else {/* la ligne o[j] a un pivot différent de zéro. On interchange les lignes
				o[k] et o[j] */
				double tmp=o[k];
				o[k] = o[j] ;
				o[j] = tmp ;
			}
		    }
		    /* On continue avec la méthode de Gauss pour triangulariser U */
		    L_ik = mat(U,n,o[i],k) / mat(U,n,o[k],k) ;
		    mat(U,n,o[i],k) = 0. ; /* cet élément doit s'annuler. On ne le calcule donc pas. */
		    for (j=k+1; j<=n ; j++)		mat(U,n,o[i],j) -= L_ik * mat(U,n,o[k],j) ;
		    for (j=1;   j<=m ; j++)		mat(B,m,o[i],j) -= L_ik * mat(B,m,o[k],j) ;
		}
	}

	/* 2ème étape: On résoud U . x = b_k, pour chaque colonne k=1,...,m */
	/* La solution pour U.x = b est : 
		pour i=n,...,1  (i est le numéro de la ligne)
			x_i = 1/U_ii ( b_i - sum_{j=i+1}^n   U_ij   x_j  )
	Les résultats seront mis dans la matrice B, avec la correspondance suivante:
		b_i  initialement égal à  B( o[i] , k )
		le résultat x_i  est ensuite mis dans B( o[i] , k )
	Ce qui revient à faire le remplacement 
	B( o[i] , k )  -->   1/U_ii * ( B(o[i] , k ) - sum_{j=i+1}^n   U_ij * B( o[j] , k ) )
	*/

	for(k=1; k<=m ; k++) {
	    for (i=n; i>=1 ; i--) {
		for (j=i+1; j<=n ; j++) 	mat(B,m,o[i],k) -= mat(U,n,o[i],j) * mat(B,m,o[j],k) ;
		mat(B,m,o[i],k) /= mat(U,n,o[i],i) ;
	    }	
	}
return;
}

void inverse_matrice( double A[], double A_inverse[], int n){
	/* calcule l'inverse de la matrice A, et met le resultat dans A_inverse */
	matrice_identite(A_inverse,n);
	gauss_multi(A,n,A_inverse,n);
}
