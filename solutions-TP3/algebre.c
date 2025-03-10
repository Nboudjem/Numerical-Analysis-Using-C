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
	/* o[i] repr�sentera la ligne "i", m�me avec interchange de lignes */
	for (i=1; i<=n ; i++) o[i] = i ;  


	/* 1�re �tape: TRIANGULARISATION de "U". 	"U" devient triangulaire sup�rieure */
	for(k=1; k<n ; k++) {/* passage de la matrice U^{k} � U^{k+1} , pour k=1,...,(n-1)  */
		double L_ik ;
		for (i=k+1; i<=n ; i++){
		    if ( mat(U,n,o[k],k) == 0. ){ /* Le pivot est nul !!!! */
			/* Cherchons un pivot diff�rent de z�ro en dessous de la ligne o[k] */
			/* Si cela n'est pas possible, cela veut dire que det(U)=det(A)=0. , et on sort alors! */
			j = k+1;
			while ( mat(U,n,o[j],k) == 0. && j <= n ) j++;
			if ( j > n ) {/* cela veut dire que l'on n'a pas trouv� de pivot != 0. ON SORT! */
				printf("\nCette matrice a un d�terminant nul. La m�thode gauss ne peut\n"
					"s'appliquer � cette matrice.\nVoici la matrice obtenue � ce stade:\n");
				afficher_matrice(U,n);
				return;
			}
			else {/* la ligne o[j] a un pivot diff�rent de z�ro. On interchange les lignes
				o[k] et o[j] */
				double tmp=o[k];
				o[k] = o[j] ;
				o[j] = tmp ;
			}
		    }
		    /* On continue avec la m�thode de Gauss pour triangulariser U */
		    L_ik = mat(U,n,o[i],k) / mat(U,n,o[k],k) ;
		    mat(U,n,o[i],k) = 0. ; /* cet �l�ment doit s'annuler. On ne le calcule donc pas. */
		    for (j=k+1; j<=n ; j++)		mat(U,n,o[i],j) -= L_ik * mat(U,n,o[k],j) ;
		    for (j=1;   j<=m ; j++)		mat(B,m,o[i],j) -= L_ik * mat(B,m,o[k],j) ;
		}
	}

	/* 2�me �tape: On r�soud U . x = b_k, pour chaque colonne k=1,...,m */
	/* La solution pour U.x = b est : 
		pour i=n,...,1  (i est le num�ro de la ligne)
			x_i = 1/U_ii ( b_i - sum_{j=i+1}^n   U_ij   x_j  )
	Les r�sultats seront mis dans la matrice B, avec la correspondance suivante:
		b_i  initialement �gal �  B( o[i] , k )
		le r�sultat x_i  est ensuite mis dans B( o[i] , k )
	Ce qui revient � faire le remplacement 
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
