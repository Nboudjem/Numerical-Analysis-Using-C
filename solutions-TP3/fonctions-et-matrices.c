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

#define mat(A,m,i,j) 	A[ ((i)-1)*(m) + (j)-1 ]
#define vec(B,i) 	B[ (i)-1 ]
/* Avec les d�finitions ci-dessus, on a �videmment: vec(B,i) == mat(B,1,i,1) */

void afficher_matrice( double A[], int n, int m){
/*
	Cette fonction affiche la matrice A sous une forme qui peut
	�tre r�-utilis�e dans un programme C.
	les �l�ments sont A_ij avec i=1, ..., n , et j=1, ..., m
*/
   int i,j;

   printf("\n...[] = {"); /* cette chaine de caract�res sera affich�e telle quelle */
	for (i=1; i<n ; i++){/* on affiche les (n-1) premi�res lignes */
		printf("\n");
		for (j=1; j<=m ; j++)	printf("   %10g,", mat(A,m,i,j) );
	}

	/* pour i=n (derni�re ligne), on affiche diff�remment */
	printf("\n");
	for (j=1; j<m ; j++)	printf("   %10g,", mat(A,m,i,j) );
	printf("   %10g\n};\n", mat(A,m,i,j) ); /* pour j=m */
}

void afficher_vecteur( double x[], int n){
	/* Cette fonction affiche le vecteur x qui a n �l�ments.
	*/
	afficher_matrice(x,n,1);
}

void afficher_vecteur_version_2( double x[], int n){
	/* Cette fonction affiche le vecteur x qui a n �l�ments.
	*/
	int i;

	for (i=1; i<=n ; i++)	printf("\n  %10g", vec(x,i) );
	printf("\n\n");
}


void produit_matrices(double A[], double B[], double C[], int n){
	/* Les matrice A et B sont carr�es,  (n x n): 
	   On calcule le produit A.B 
	   Les �l�ments de A.B sont stock�s dans C.
	*/
   int i,j,k;

   for (i=1; i<=n ; i++){
	for (j=1; j<=n ; j++){
		mat(C,n,i,j)  = 0. ;
		for(k=1;k<=n;k++) 	mat(C,n,i,j) += mat(A,n,i,k) * mat(B,n,i,j) ;
	}
   }
return;
}

void produit_matrice_vecteur(double A[], int n, double x[], double y[]){
	/* La matrice A est carr�e,  (n x n), les vecteurs x et y ont n �l�ments.
	   On calcule le produit y = A.x
	*/
   int i,k;

   for (i=1; i<=n ; i++){
	y[i]  = 0. ;
	for(k=1;k<=n;k++) y[i] += mat(A,n,i,k) * vec(x,k) ;
   }
return;
}

void concat_matrices(double A[], int n, int m, double B[], int n_prime, int m_prime, double C[]){
	/* Les matrices A (n x m) et B (n x m_prime) sont copi�es
	   dans une matrice plus grande C qui contient (m+m_prime) colonnes.
	   Donc on peut ecrire C =  A | B, voulant dire que
		C_ij = A_ij si j <= m
		C_ij = B_ip  pour j=m+1, ..., (m+m_prime), avec p = j - m
	Cela revient juste � agrandir par la droite une matrice (souvent carr�e)
	avec une autre matrice, et donc avoir un programme simple pour gauss
	et gauss-jordan, qui  manipulera les �l�ments d'une seule matrice
	au lieu de 2 (dans le cas par exemple du calcul de l'inverse d'une matrice),
	ou dans le cas de la r�solution de A.x = b , avec A=matrice et b=vecteur:
	dans ce dernier cas, on �crit une seule matrice C = A | b, qui a une colonne
	de plus que A.
	Remarque:	la matrice C doit avoir �t� d�clar�e convenablement avant
			d'appeler cette fonction
	*/
   int i,j,M=m+m_prime;
   for (i=1; i<=n ; i++){
	for (j=1; j<=m ; j++) 		mat(C,M,i,j) = mat(A,   m   , i ,  j )  ;
	for (j=m+1; j<=M ; j++) 	mat(C,M,i,j) = mat(B,m_prime, i , j-m)  ;
   }
return;
}

void split_matrices(double A[], int n, int m, double B[], int n_prime, int m_prime, double C[]){
	/* Cette fonction est l'inverse de concat_matrices d�finie plus haut.
	   Les matrices A (n x m) et B (n x m_prime) sont copi�es
	   a partir d'une matrice plus grande C qui contient (m+m_prime) colonnes.
	   Donc on peut ecrire C =  A | B, voulant dire que
		C_ij = A_ij si j <= m
		C_ij = B_ip  pour j=m+1, ..., (m+m_prime), avec p = j - m
	*/
   int i,j,M=m+m_prime;
   for (i=1; i<=n ; i++){
	for (j=1; j<=m ; j++) 		mat(A,   m   , i ,  j ) = mat(C,M,i,j) ;
	for (j=m+1; j<=M ; j++) 	mat(B,m_prime, i , j-m) = mat(C,M,i,j)  ;
   }
return;
}

void sous_matrices(double A[], int n, int m, double B[], int l1, int l2, int c1, int c2){
/*
	La matrice A (n x m) est donn�e. Nous copions un rectangle que nous mettons dans 
	la matrice B. B aura (lll = l2-l1+1) lignes et (ccc =c2-c1+1) colonnes, avec:
		mat(B,ccc,p,q)= mat(A,m, i , j )
		o� p=1,..., lll    et  q=1,..., ccc
		et: p = i-l1 + 1	; 	q = j-c1 + 1
	*/
   int i,j, lll=l2-l1+1, ccc=c2-c1+1;
   for (i=l1; i<=l2 ; i++){
	for (j=c1; j<=c2 ; j++) mat(B, ccc , i-l1+1 ,  j-c1+1 ) = mat(A,m,i,j) ;
   }
return;
}



void gauss_LU(double A[], double L[], int n, char affichage){
	/* Les matrice A et L sont carr�es,  (n x n): 
	   Cette fonction transforme la matrice A en un produit L.U,
	   U est triangulaire sup�rieure,
	   L est triangulaire inf�rieure.
	   Les �l�ments de U sont stock�s dans A.
	Si ( affichage == 'y') , alors on affiche la matrice A^{k}, pour k=1,...,n
	*/
   int i,j,k;
   
   /* On initialise d'abord la matrice L � la matrice Identit� */
   for (i=1; i<=n ; i++){
   	for (j=1; j<=n ; j++){
   		if (i==j) 	mat(L,n,i,j) = 1. ;
   		else		mat(L,n,i,j) = 0. ;
   	}
   }

   if (affichage == 'y'){
	printf("A^1 = \n");
	afficher_matrice(A,n,n);/* affichons la matrice originale A^1 */
   }

   for(k=1; k<n ; k++) {/* passage de la matrice A^{k} � A^{k+1} (sans interchange de lignes)*/
	for (i=k+1; i<=n ; i++){
		mat(L,n,i,k) = mat(A,n,i,k) / mat(A,n,k,k) ;
		for (j=k; j<=n ; j++)	mat(A,n,i,j) -= mat(L,n,i,k) * mat(A,n,k,j) ;
		/*
		La ligne pr�c�dente peut etre remplac�e par les 2 lignes suivantes:
		mat(A,n,i,k) = 0. ;
		for (j=k+1; j<=n ; j++)		mat(A,n,i,j) -= mat(L,n,i,k) * mat(A,n,k,j) ;
		*/
	}

	if (affichage == 'y'){/* affichons la matrice A^{k+1} que nous venons de calculer */
		printf("A^%d = \n" , k+1);
		afficher_matrice(A,n,n);
	}
   }
return;
}

void gauss_jordan_simple_v1(double A[], int n, double b[], char affichage){
	/* La matrice A est carr�e,  (n x n): 
	   Cette fonction transforme la matrice A en une matrice diagonale,
	   Les �l�ments (resulat) sont stock�s dans A.
	Ensuite, on r�soud le syst�me, et on met le r�sultat dans le vecteur b
	ICI, nous n'interchangeons aucune ligne .
	Si (affichage == 'y') nous affichons A^k et b^k pour k=1,...,n
	*/
   int i,j,k;
   double L_ik;
   
   if (affichage == 'y'){
	printf("A^1 = \n");
	afficher_matrice(A,n,n);/* affichons la matrice originale A^1 */
	printf("b^1 = \n");
	afficher_vecteur(b,n);/* affichons le vecteur original b^1 */
   }

   for(k=1; k<=n ; k++) {/* passage de la matrice A^{k} � A^{k+1} */
	for (i=1; i<=n ; i++){
		L_ik = mat(A,n,i,k) / mat(A,n,k,k) ;

		if (i==k) continue;
		for (j=1; j<=n ; j++){
			if (j==k)	mat(A,n,i,j)  = 0. ;
			else 		mat(A,n,i,j) -= L_ik * mat(A,n,k,j) ;
		}

		/* On tansforme aussi le vecteur b */
		b[i] -= L_ik * b[k] ;
	}
	if (affichage == 'y'){
		printf("A^%d = \n",k+1);
		afficher_matrice(A,n,n);/* affichons la matrice A^{k+1} */
		printf("b^%d = \n",k+1);
		afficher_vecteur(b,n);/* affichons le vecteur original b^{k+1} */
	}
   }

   /* On r�soud maintenant le syst�me qui est diagonal */
   for (i=1; i<=n ; i++)	b[i] = b[i] / mat(A,n,i,i) ;

return;
}

void main(int argc, char *argv[]){
	int i,j;
	long n=0;
	double a,b,c, fa,fb,fc, Delta, Delta_relative , eps=1.e-15;
	double L[9], U[9] , xxx[5] , yyy[9];

	double A[] = 
		{
			2,  3, -1 ,
			1, -2,  5 ,
			6, -3, 15
		} ;


	double B[] = 
		{
			2,  3, -1 , 4,
			1, -2,  5 , 3,
			6, -3, 15 , -1,
			-2, 3, 0 , 1
		} ;


	double C[] = 
		{
			2,  3, -1 , 4,
			1, -2,  5 , 3,
			6, -3, 15 , -1,
			-2, 3, 0 , 1
		} ;

	double bbb[] ={
		0.,
		1,
		3,
		5,
		-1
	};

	double ccc[] ={
		0.,
		1,
		3,
		5,
		-1
	};


	double D[] = 
		{
			 2,  3, -1 ,   4,  1, -2,  5, 32,
			 1, -2,  5 ,   3, -2,  3,  0, 10,
			 6, -3, 15 ,  -1,  1,  2, -5,-11,
			-2,  3,  0 ,   1, +3,  4, -2,-37,
			 1, -1,  7 ,   0, 13,  1,  2, 27,
			 2,  2, 13 ,   7,-31,  4, -6,  17,
			 3, -3, 11 ,  18, 43,  7, -8, -14,
			 4, -4,-17,  -23, 37, -9,  3,  15

		} ;
	double ddd[] ={
		 1,
		 3,
		 5,
		-1,
		23,
		13,
		-67,
		99
	};
	double E[64] = 
		{
			 2,  3, -1 ,   4,  1, -2,  5, 32,
			 1, -2,  5 ,   3, -2,  3,  0, 10,
			 6, -3, 15 ,  -1,  1,  2, -5,-11,
			-2,  3,  0 ,   1, +3,  4, -2,-37,
			 1, -1,  7 ,   0, 13,  1,  2, 27,
			 2,  2, 13 ,   7,-31,  4, -6,  17,
			 3, -3, 11 ,  18, 43,  7, -8, -14,
			 4, -4,-17,  -23, 37, -9,  3,  15

		} ;
	double eee[] ={
		 1,
		 3,
		 5,
		-1,
		23,
		13,
		-67,
		99
	};
	/*
		VOICI comment on peut initialiser une matrice 4x4 avec des zeros:
		for (i=1; i<=4 ; i++){
			for (j=1; j<=4 ; j++)	mat(A,4,i,j) = 0. ;
		}

		VOICI comment on peut initialiser une matrice (n x n) � la matrice identit�:
		for (i=1; i<=n ; i++){
			for (j=1; j<=n ; j++){
				if (i==j) 	mat(A,n,i,j) = 1. ;
				else 		mat(A,n,i,j) = 0. ;
			}
		}
	*/

	double DDD[12]; /* 3 lignes, 4 colonnes */
	afficher_matrice(D,8,8);
	sous_matrices(D,8,8,ddd,1,8,8,8);
	afficher_matrice(ddd,8,1);

	double F[72];

	afficher_matrice(D,8,8);
	afficher_vecteur(ddd,8);
	concat_matrices(D,8,8,ddd,8,1,F);
	afficher_matrice(F,8,9);

	split_matrices(D,8,8,ddd,8,1,F);
	afficher_matrice(D,8,8);
	afficher_vecteur(ddd,8);



	gauss_jordan_simple_v1(D,8,ddd,'y');
	/*	afficher_matrice(D,8,8); */
	printf("\nLa solution est le vecteur x =\n");
	afficher_vecteur(ddd,8);

	/* V�rifions les r�sultats pr�c�dents */
	printf("\nV�rifions notre solution D.x =\n");
	produit_matrice_vecteur(E,8,ddd,yyy);

	afficher_vecteur(yyy,8);



	afficher_matrice(B,4,4);
	afficher_vecteur(bbb,4);

	gauss_jordan_simple_v1(B,4,bbb,'y');
	afficher_matrice(B,4,4);
	afficher_vecteur(bbb,4);

	/* V�rifions les r�sultats pr�c�dents */
	produit_matrice_vecteur(C,4,bbb,xxx);

	afficher_vecteur(xxx,4);

	return;
}

