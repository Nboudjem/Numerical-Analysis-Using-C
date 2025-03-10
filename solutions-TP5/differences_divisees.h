#include <stdio.h>
#include <math.h>
#include <stdarg.h>

#define mat(B,n,i,j) 	B[ (i)*(n) + (j) ] /* element B_ij , B matrice (n+1)x(n+1) , i,j = 0,...,n */

void differences_divisees(int n, double X[], double Y[], double B[]){
/* ***************************************************************************
	Calcul de la table des diff�rences divis�es B(i,j) lorsque les X(i) et les Y(i) 
	sont suppos�s donn�s: i,j=0,...,n. Les r�sultats sont dans la matrice B, (n+1)x(n+1),
	que l'on aura auparavant d�clar�e correctement dans la fonction  principale "main".

	RAPPEL: la table des diff. divis�es est calcul�e comme ceci :
		B(i,j) == [Y_i, Y_{i+1}, ... , Y_{j-1}, Y_j] 
		Relation de r�currence pour les B(i,j):
			B(i,j) = ( B(i+1,j) - B(i,j-1) )   /   ( X(j) - X(i) )  avec i < j.
			avec, par d�finition, B(i,i) = Y(i)

	1. On calcule d'abord B(i,i+1) pour i=0,...,n-1
	2. On calcule ensuite B(i,i+2) pour i=0,...,n-2
	3. ........
	4. On calcule finalement B(i,n). pour i=0,...,0 (une seule valeur)
   ***************************************************************************/

	int i,j,l;

	for(i=0 ; i<= n ; i++)	mat(B,n,i,i) = Y[i] ; /* initialisation de B(i,i) = Y_i */

	for(l=1 ; l<= n ; l++){
		for(i=0 ; i<= n-l ; i++){
			j = i+l ;
			mat(B,n,i,j) = ( mat(B,n,i+1,j) - mat(B,n,i,j-1) )/ ( X[j] - X[i] ) ;
		}
	}
	return;
}


void afficher_differences_divisees(int n, double B[], char *sss){
	/* Affichons la table B_ij (uniquement les �lements calcul�s: i<= j ) 
		le dernier argument "sss", par exemple de la forme "%16.8g",
		indiquera comment formater les B_ij lors de l'affichage sur �cran.
	*/
	int i,j;
	printf("\nB = ");
	for(i=0 ; i<= n ; i++){
	printf("\n\t");
		for(j=0 ; j<= n ; j++) {
			if (i>j)	printf("%16.0f", 0. ); /* Les �l�ments non calcul�s sont pos�s egaux � z�ro */
			else 		printf(sss,mat(B,n,i,j) );
		}
	}
	printf("\n");
	return;
}


double Polynome_Newton1(double x, double X[], int n, double B[] ){
	/*
	Cette fonction calcule P_n(x), 1�re formule d'interpolation polynomiale de Newton.
	La matrice B, (n+1)x(n+1), est la table des diff�rences divis�es, suppos�ee donn�,
	et les X[i] sont aussi suppos�s donn�s.

	P_n(x) = a0     + a1 * (x - X[0])
			+ a2 * (x - X[0]) * (x - X[1])
			+ ...
			+ an * (x - X[0]) * (x - X[1]) * ... * (x-X[n-1])

	Les coefficients ai sont les differences divis�es, �gales �: a_i = B_{0,i} pour la
	1�re formule d'interpolation de Newton.
	*/
	int i;
	double y=1. , sum = mat(B,n,0,0) ; /* sum = a0 : initialisation */
	for (i=1 ; i<= n ; i++){
		y *= x-X[i-1] ;
		sum += mat(B,n,0,i) * y ;
	}	
    return sum;
}

double Polynome_Newton2(double x, double X[], int n, double B[] ){
	/*
	Cette fonction calcule P_n(x), 2�me formule d'interpolation polynomiale de Newton.
	La matrice B, (n+1)x(n+1), est la table des diff�rences divis�es, suppos�ee donn�,
	et les X[i] sont aussi suppos�s donn�s.

	P_n(x) = b0     + b1 * (x - X[n])
			+ b2 * (x - X[n]) * (x - X[n-1])
			+ ...
			+ bn * (x - X[n]) * (x - X[n-1]) * ... * (x-X[1])

	Les coefficients bi sont les differences divis�es, �gales �: b_i = B_{i,n} pour la
	2�me formule d'interpolation de Newton.
	*/
	int i;
	double y=1. , sum = mat(B,n,n,n) ; /* sum = b0 : initialisation */
	for (i=n ; i>= 1 ; i--){
		y *= x-X[i] ;
		sum += mat(B,n,i-1,n) * y ;
	}	
    return sum;
}

double Polynome_Legendre(double x, double X[], double Y[], int n){
	/*
	Cette fonction calcule le polyn�me de Legendre au point "x".
	Nous supposons donn�s les X[i] and Y[i], i=0,...,n.

	P_n(x) = \sum_{i=0}^{n} Y[i] * \Pi_{j=0,...,n} (x-X[j])/(X[i]-X[j]) ;
		dans le produit \Pi, "j" doit �tre diff�rent de "i". 

	*/
	int i,j;
	double z , sum = 0. ; /* sum = 0 : initialisation ; "z" contiendra le produit \Pi */
	for (i=0 ; i<= n ; i++){ /* On �crit 2 sommes pour �viter j=i */
		z = 1. ;
		for (j=0 ; j<i ; j++)		z *= (x-X[j])/(X[i]-X[j]) ;
		for (j=i+1 ; j<=n ; j++)	z *= (x-X[j])/(X[i]-X[j]) ;
		/* les 2 boucles "for" pr�c�dentes auraient pu �tre �crite
		comme une seule boucle comme ceci:
		for (j=0 ; j<=n ; j++)	if (j!=i)	z *= (x-X[j])/(X[i]-X[j]) ;
		*/
		sum += Y[i] * z ;
	}	
    return sum;
}


int Read_X_Y(char *sss, double X[], double Y[] , int N){
	/*
	Read at most (N+1) data points from file "sss" and put results into
	arrays X[i] and Y[i], i=0,1,...,n.  X and Y are supposed to have dimension >= (N+1).
	Each line is supposed to contain at least 2 "double" numbers, else it will simply be ignored.
	This function returns "n" if (n+1) data points are found, else returns -1 if no data or if file cannot
	be opened for reading.
	*/
	FILE *in;
	int i, n, m, c;
	char line[1000]; /* chaque ligne lue sera mise dans le tableau line. On scan ensuite chaque ligne
			pour voir si cette ligne contient deux float, sinon, cette ligne est ignor�e */

	in = fopen(sss, "r") ;
	if ( in == NULL) return -1;

	n=0;
	while( n <= N ){
		/* we read each line and put it to line[1000], 
		we then scan the numbers and put them into arrays X[n] and Y[n] */
		i=0;
		while ( (c = fgetc(in) ) != '\n' && c != EOF && i <1000) line[i++] = c;
		line[i++] = 0 ;

		/* we now scan line[] */
		m = sscanf(line, "%lf %lf" , &X[n] , &Y[n] );
		if ( m==2 ) n++; /* we found 2 numbers */
		if ( c==EOF) break ;
	}

	fclose(in);
	
    return n-1;
}

int Read_X_Y_DY(char *sss, double X[], double Y[], double DY[] , int N){
	/*
	Read at most (N+1) data points from file "sss" and put results into
	arrays X[i], Y[i] and DY[i] (les DY[i] sont les erreurs sur Y[i], i=0,1,...,n.
	X, Y and DY are supposed to have dimension >= (N+1).
	Each line is supposed to contain at least 3 float numbers, else it will simply be ignored.
	This function returns "n" if (n+1) data points are found, else returns -1 if no data or if file cannot
	be opened for reading.
	*/
	FILE *in;
	int i, n, m, c; 
	char line[1000]; /* chaque ligne lue sera mise dans le tableau line. On scan ensuite chaque ligne
			pour voir si cette ligne contient 3 "doubles", sinon, cette ligne est ignor�e */

	in = fopen(sss, "r") ;
	if ( in == NULL) return -1;

	n=0;
	while( n <= N ){
		/* we read each line and put it to line[1000], 
		we then scan the numbers and put them into arrays X[n] and Y[n] */
		i=0;
		while ( (c = fgetc(in) ) != '\n' && c != EOF && i <1000) line[i++] = c;
		line[i++] = 0 ;

		/* we now scan line[] */
		m = sscanf(line, "%lf %lf %lf" , &X[n] , &Y[n] , &DY[n]);
		if ( m==3 ) n++; /* we found 3 numbers */
		if ( c==EOF) break ;
	}

	fclose(in);
	
    return n-1;

}





void moindres_carres_ax(int n, double x[], double y[], double *a ){
	/* 
	Nous supposons que les tableaux x[i] et y[i], qui contiennent des "doubles"
	sont donn�s. Nous utiliserons uniquement les premiers points i=0,...,n , car
	"n-1" est suppos� <= � la dimension de x[] et � celle de y[].
	Nous cherchons ici:
		f(x) = a*x
	en utilisant les moindres carr�s.
	R�sultat : a = ( sum_0^n { x[i]*y[i] }  ) / ( sum_0^n { x[i]*x[i] }  )
	*/
	int i;
	double xy_moyen , x2_moyen ;
	
	xy_moyen = x[0]*y[0] ;
	x2_moyen = x[0]*x[0] ; /* initialisation des sommes � calculer */
	for(i=1 ; i <= n ; i++) 
		{
		xy_moyen += x[i]*y[i] ;
		x2_moyen += x[i]*x[i] ;
		}
	*a = xy_moyen / x2_moyen ;
return;
}

void moindres_carres_ax_b(int n, double x[], double y[], double *a , double *b){
	/* 
	Nous supposons que les tableaux x[i] et y[i], qui contiennent des "doubles"
	sont donn�s. Nous utiliserons uniquement les "n+1" premiers points i=0,...,n , car
	les dimensions de x[] et y[] sont suppos�es >= "n+1" .
	Nous cherchons ici:
		f(x) = a*x + b
	en utilisant les moindres carr�s.
	*/
	int i;
	double xy_moyen , x2_moyen, x_moyen , y_moyen ; 
		/* D�finitions utilis�es ici, en g�n�ral: X_moyen == \sum_0^n X_i */
	
	x_moyen = x[0] ;
	y_moyen = y[0] ;
	xy_moyen = x[0]*y[0] ;
	x2_moyen = x[0]*x[0] ; /* initialisation des sommes � calculer */
	for(i=1 ; i <= n ; i++) 
		{
		x_moyen += x[i] ;
		y_moyen += y[i] ;
		xy_moyen += x[i]*y[i] ;
		x2_moyen += x[i]*x[i] ;
		}
	*a = ( xy_moyen - x_moyen * y_moyen/(double)(n+1) ) 
		/ ( x2_moyen - x_moyen * x_moyen/(double)(n+1) ) ;
	*b = ( y_moyen - (*a) * x_moyen ) / (double)(n+1) ; 
return;
}

void moindres_carres_ax2_bx_c(int n, double x[], double y[], double *a , double *b, double *c){
	/* 
	Nous supposons que les tableaux x[i] et y[i], qui contiennent des "doubles"
	sont donn�s. Nous utiliserons uniquement les premiers points i=0,...,n , car
	"n-1" est suppos� <= � la dimension de x[] et � celle de y[].
	Nous cherchons ici:
		f(x) = a*x^2 + b*x + c
	en utilisant les moindres carr�s.
	Nous devons r�soudre le syst�me de 3 �quations � 3 inconnues (a,b,c):
		U11 a + U12 b + U13 c = V1		(U est sym�trique et sp�ciale)
		U21 a + U22 b + U23 c = V2
		U31 a + U32 b + U33 c = V3
	avec les d�finitions suivantes :
		--------------------------------------------------------------------
		U11 = \sum x_i^4 ;	U12 = \sum x_i^3 ;	U13 = \sum x_i^2	
		U21 = U12 ;		U22 = U13 ;		U23 = \sum x_i
		U31 = U13 ;		U32 = U23 ;		U33 = \sum 1 = (n+1)
		--------------------------------------------------------------------
		V1 = \sum 	x_i^2	y_i
		V2 = \sum 	x_i	y_i
		V3 = \sum		y_i
		--------------------------------------------------------------------
	*/
	int i;
	double x_1, x_2, x_3, x_4, y_1,Delta;
	double  U11,U12,U13,V1 ; 
	double  U21,U22,U23,V2 ; 
	double  U31,U32,U33,V3 ; 
	
	U11 = 0. ;	U12 = 0. ;	U13 = 0. ;	U23 = 0. ;
	V1 = 0. ;
	V2 = 0. ;
	V3 = 0. ; /* initialisation des sommes � calculer */
	for(i=0 ; i <= n ; i++) 
		{
		x_1 = x[i];		y_1 = y[i] ;
		x_2 = x_1*x_1;
		x_3 = x_2*x_1;
		x_4 = x_2*x_2;
		/* -------------------------------------------*/
		U11 += x_4 ;	U12 += x_3 ;	U13 += x_2 ;	U23 += x_1 ;
		V1 += x_2 * y_1 ;
		V2 += x_1 * y_1 ;
		V3 +=       y_1 ;
		}
	/* ------ On calcule maintenant les autres U_ij  ---------*/
	U21 = U12 ;		U22 = U13 ;		
	U31 = U13 ;		U32 = U23 ;		U33 = (double)(n+1);
	
	/* On r�soud maintenant pour a, b, c , avec les d�terminants */
	/* determinant_principal = Delta =
		| U11  U12  U13 |
		| U21  U22  U23 |
		| U31  U32  U33 |
		= 	U11 * (U22*U33 - U23*U32) 
		      + U12 * (U23*U31 - U21*U33) 
		      + U13 * (U21*U32 - U22*U31) 
	DONC: a = 1/Delta *
		| V1  U12  U13 |
		| V2  U22  U23 |
		| V3  U32  U33 |
	DONC: b = 1/Delta *
		| U11  V1  U13 |
		| U21  V2  U23 |
		| U31  V3  U33 |
	De m�me pour c.
	*/

	Delta	= 	U11 * (U22*U33 - U23*U32) 
		      + U12 * (U23*U31 - U21*U33) 
		      + U13 * (U21*U32 - U22*U31) ;

	/* -------------------------------------*/
	*a	= 	V1 * (U22*U33 - U23*U32) 
		      + V2 * (U23*U31 - U21*U33) 
		      + V3 * (U21*U32 - U22*U31) ;
	*a /= Delta ;

	/* -------------------------------------*/
	*b	= 	V1 * (U23*U31 - U21*U33) 
		      + V2 * (U11*U33 - U13*U31) 
		      + V3 * (U13*U21 - U11*U23) ;
	*b /= Delta ;

	/* -------------------------------------*/
	/* pour c, il est plus simple d'utiliser la 3�me �quation */
	*c = ( V3 - (*a)*U31 - (*b)*U32 )/(double)(n+1);
return;
}

void ki2_ax(int n, double x[], double y[], double Dy[], double *a ){
	/* 
	Nous supposons que les tableaux x[i] et y[i], ainsi que les erreurs Dy[i]
	qui contiennent des "doubles" sont donn�s. Nous utiliserons uniquement
	les premiers points i=0,...,n , car
	"n-1" est suppos� <= aux dimensions de x[], y[] et Dy
	Nous cherchons ici:
		f(x) = a*x
	en utilisant les moindres carr�s pond�r�s (m�thode du ki^2)
	*/
	int i;
	double xy_moyen , x2_moyen , w_i;
	
	xy_moyen = 0. ;
	x2_moyen = 0. ; /* initialisation des sommes � calculer */
	for(i=0 ; i <= n ; i++) 
	{
		if ( Dy[i] != 0. ) w_i = 1./Dy[i]/Dy[i] ;
		else {
			printf("\n\nERREUR: Division par z�ro dans le calcul des poids.\nExiting ...\n\n");
			exit(1); 
		}
		xy_moyen += x[i]*y[i]*w_i ;
		x2_moyen += x[i]*x[i]*w_i ;
	}
	*a = xy_moyen / x2_moyen ;
return;
}


void ki2_ax_b(int n, double x[], double y[], double Dy[], double *a , double *b){
	/* 
	Nous supposons que les tableaux x[i] et y[i] (et les erreurs Dy[]), 
	qui contiennent des "doubles" sont donn�s. Nous utiliserons uniquement
	les premiers points i=0,...,n , car "n-1" est suppos� <= 
	aux dimensions de x[], y[] et Dy[].
	Nous cherchons ici:
		f(x) = a*x + b
	en utilisant les moindres carr�s.
	Nous devons r�soudre le syst�me de 3 �quations � 3 inconnues (a,b,c):
		U11 a + U12 b = V1		(U est sym�trique et sp�ciale)
		U21 a + U22 b = V2
	avec les d�finitions suivantes :
		--------------------------------------------------------------------
		U11 = \sum w_i*x_i^2 ;	U12 = \sum w_i*x_i
		U21 = U12 ;		U22 = \sum w_i
		--------------------------------------------------------------------
		V1 = \sum 	w_i	x_i	y_i
		V2 = \sum	w_i	y_i
		--------------------------------------------------------------------
	*/
	int i;
	double x_1, x_2, y_1, w_i,Delta;
	double  U11,U12,V1 ; 
	double  U21,U22,V2 ; 
	
	U11 = 0. ;	U12 = 0. ;	U22 = 0. ;
	V1 = 0. ;
	V2 = 0. ;	/* initialisation des sommes � calculer */
	for(i=0 ; i <= n ; i++) 
		{
		if ( Dy[i] != 0. ) w_i = 1./Dy[i]/Dy[i] ;
		else {
			printf("\n\nERREUR: Division par z�ro dans le calcul des poids.\nExiting ...\n\n");
			exit(1); 
		}
		x_1 = x[i];		y_1 = y[i] ;
		x_2 = x_1*x_1;
		/* -------------------------------------------*/
		U11 += w_i*x_2 ;
		U12 += w_i*x_1 ;
		U22 += w_i ;
		V1 += w_i* x_1 * y_1 ;
		V2 += w_i*       y_1 ;
		}
	/* ------ On calcule maintenant les autres U_ij  ---------*/
	U21 = U12 ;
	
	/* On r�soud maintenant pour a, b , avec les d�terminants */
	/* determinant_principal = Delta =
		| U11  U12 |
		| U21  U22 |
		= U11 * U22  - U12*U21 ;
	DONC: a = 1/Delta *
		| V1  U12 |
		| V2  U22 |
		= V1*U22 - V2*U12 ;
	DONC: b = 1/Delta *
		| U11  V1 |
		| U21  V2 |
		= V2*U11 - V1*U21 ;
	*/

	Delta	=  U11 * U22  - U12*U21 ;

	/* -------------------------------------*/
	*a	=  V1*U22 - V2*U12 ;
	*a /= Delta ;

	/* -------------------------------------*/
	*b	=  V2*U11 - V1*U21 ;

	*b /= Delta ;
return;
}


void ki2_ax2_bx_c(int n, double x[], double y[], double Dy[], double *a , double *b, double *c){
	/* 
	Nous supposons que les tableaux x[i] et y[i] (et les erreurs Dy[]), 
	qui contiennent des "doubles" sont donn�s. Nous utiliserons uniquement
	les premiers points i=0,...,n , car "n-1" est suppos� <= 
	aux dimensions de x[], y[] et Dy[].
	Nous cherchons ici:
		f(x) = a*x^2 + b*x + c
	en utilisant les moindres carr�s.
	Nous devons r�soudre le syst�me de 3 �quations � 3 inconnues (a,b,c):
		U11 a + U12 b + U13 c = V1		(U est sym�trique et sp�ciale)
		U21 a + U22 b + U23 c = V2
		U31 a + U32 b + U33 c = V3
	avec les d�finitions suivantes :
		--------------------------------------------------------------------
		U11 = \sum w_i*x_i^4 ;	U12 = \sum w_i*x_i^3 ;	U13 = \sum w_i*x_i^2	
		U21 = U12 ;		U22 = U13 ;		U23 = \sum w_i*x_i
		U31 = U13 ;		U32 = U23 ;		U33 = \sum w_i
		--------------------------------------------------------------------
		V1 = \sum 	w_i	x_i^2	y_i
		V2 = \sum 	w_i	x_i	y_i
		V3 = \sum	w_i	y_i
		--------------------------------------------------------------------
	*/
	int i;
	double x_1, x_2, x_3, x_4, y_1, w_i, Delta;
	double  U11,U12,U13,V1 ; 
	double  U21,U22,U23,V2 ; 
	double  U31,U32,U33,V3 ; 
	
	U11 = 0. ;	U12 = 0. ;	U13 = 0. ;	U23 = 0. ;	U33 = 0. ;
	V1 = 0. ;
	V2 = 0. ;
	V3 = 0. ; /* initialisation des sommes � calculer */
	for(i=0 ; i <= n ; i++) 
		{
		if ( Dy[i] != 0. ) w_i = 1./Dy[i]/Dy[i] ;
		else {
			printf("\n\nERREUR: Division par z�ro dans le calcul des poids.\nExiting ...\n\n");
			exit(1); 
		}
		x_1 = x[i];		y_1 = y[i] ;
		x_2 = x_1*x_1;
		x_3 = x_2*x_1;
		x_4 = x_2*x_2;
		/* -------------------------------------------*/
		U11 += w_i*x_4 ;
		U12 += w_i*x_3 ;
		U13 += w_i*x_2 ;
		U23 += w_i*x_1 ;
		U33 += w_i ;
		V1 += w_i* x_2 * y_1 ;
		V2 += w_i* x_1 * y_1 ;
		V3 += w_i*       y_1 ;
		}
	/* ------ On calcule maintenant les autres U_ij  ---------*/
	U21 = U12 ;		U22 = U13 ;
	U31 = U13 ;		U32 = U23 ;
	
	/* On r�soud maintenant pour a, b, c , avec les d�terminants */
	/* determinant_principal = Delta =
		| U11  U12  U13 |
		| U21  U22  U23 |
		| U31  U32  U33 |
		= 	U11 * (U22*U33 - U23*U32) 
		      + U12 * (U23*U31 - U21*U33) 
		      + U13 * (U21*U32 - U22*U31) 
	DONC: a = 1/Delta *
		| V1  U12  U13 |
		| V2  U22  U23 |
		| V3  U32  U33 |
	DONC: b = 1/Delta *
		| U11  V1  U13 |
		| U21  V2  U23 |
		| U31  V3  U33 |
	DONC: c = 1/Delta *
		| U11  U12  V1 |
		| U21  U22  V2 |
		| U31  U32  V3 |
	*/

	Delta	= 	U11 * (U22*U33 - U23*U32) 
		      + U12 * (U23*U31 - U21*U33) 
		      + U13 * (U21*U32 - U22*U31) ;

	/* -------------------------------------*/
	*a	= 	V1 * (U22*U33 - U23*U32) 
		      + V2 * (U23*U31 - U21*U33) 
		      + V3 * (U21*U32 - U22*U31) ;
	*a /= Delta ;

	/* -------------------------------------*/
	*b	= 	V1 * (U23*U31 - U21*U33) 
		      + V2 * (U11*U33 - U13*U31) 
		      + V3 * (U13*U21 - U11*U23) ;
	*b /= Delta ;

	/* -------------------------------------*/
	*c	= 	V1 * (U21*U32 - U22*U31) 
		      + V2 * (U12*U31 - U11*U32) 
		      + V3 * (U11*U22 - U12*U21) ;
	*c /= Delta ;
return;
}

void moindres_carres_ah1_bh2_ch3(int n, double x[], double y[], 
	double (*h1)(), double (*h2)(), double (*h3)(), 
	double *a , double *b, double *c){
	/* 
	Nous supposons que les tableaux x[i] et y[i], 
	qui contiennent des "doubles" sont donn�s. Nous utiliserons uniquement
	les premiers points i=0,...,n , car "n-1" est suppos� <= 
	aux dimensions de x[] et y[].
	Nous cherchons ici:
		f(x) = a*h1(x) + b*h2(x) + c*h3(x)
	en utilisant les moindres carr�s. Remarquez que h1, h2 et h3 sont
	des pointeurs � des fonctions; les fonction seront (*h1)(x), ....

	Nous devons r�soudre le syst�me de 3 �quations � 3 inconnues (a,b,c):
		U11 a + U12 b + U13 c = V1		(U est sym�trique)
		U21 a + U22 b + U23 c = V2
		U31 a + U32 b + U33 c = V3
	avec les d�finitions suivantes : ici tous les w_i = 1 ;
		--------------------------------------------------------------------
		U11 = \sum  h1(x_i) * h1(x) ;
		U12 = \sum  h1(x_i) * h2(x) ;		
		U13 = \sum  h1(x_i) * h3(x) ;

		U22 = \sum  h2(x_i) * h2(x) ;
		U23 = \sum  h2(x_i) * h3(x) ;

		U33 = \sum  h3(x_i) * h3(x) ;

		Les autres U_ij sont :
		U21 = U12 ;
		U31 = U13 ;
		U32 = U23 ;
		--------------------------------------------------------------------
		V1 = \sum  h1(x_i) * y_i ;
		V2 = \sum  h2(x_i) * y_i ;
		V3 = \sum  h3(x_i) * y_i ;
		--------------------------------------------------------------------
	*/
	int i;
	double x_1, y_1,  h1_i , h2_i , h3_i, Delta;
	double  U11,U12,U13,V1 ; 
	double  U21,U22,U23,V2 ; 
	double  U31,U32,U33,V3 ; 
	
	U11 = 0. ;	U12 = 0. ;	U13 = 0. ;
	U23 = 0. ;	U23 = 0. ;
	U33 = 0. ;
	V1 = 0. ;	V2 = 0. ;	V3 = 0. ; /* initialisation des sommes � calculer */

	for(i=0 ; i <= n ; i++) {
		x_1 = x[i];		y_1 = y[i] ;
		h1_i = (*h1)(x_1) ;
		h2_i = (*h2)(x_1) ;
		h3_i = (*h3)(x_1) ;
		/* -------------------------------------------*/
		U11 +=  h1_i * h1_i ;
		U12 +=  h1_i * h2_i ;
		U13 +=  h1_i * h3_i ;

		U22 +=  h2_i * h2_i ;
		U23 +=  h2_i * h3_i ;

		U33 +=  h3_i * h3_i ;

		V1  +=  h1_i * y_1 ;
		V2  +=  h2_i * y_1 ;
		V3  +=  h3_i * y_1 ;
		}
	/* ------ On calcule maintenant les autres U_ij : U est sym�trique ---------*/
	U21 = U12 ;
	U31 = U13 ;
	U32 = U23 ;
	
	/* On r�soud maintenant pour a, b, c , avec les d�terminants */
	/* determinant_principal = Delta =
		| U11  U12  U13 |
		| U21  U22  U23 |
		| U31  U32  U33 |
		= 	U11 * (U22*U33 - U23*U32) 
		      + U12 * (U23*U31 - U21*U33) 
		      + U13 * (U21*U32 - U22*U31) 
	DONC: a = 1/Delta *
		| V1  U12  U13 |
		| V2  U22  U23 |
		| V3  U32  U33 |
	DONC: b = 1/Delta *
		| U11  V1  U13 |
		| U21  V2  U23 |
		| U31  V3  U33 |
	DONC: c = 1/Delta *
		| U11  U12  V1 |
		| U21  U22  V2 |
		| U31  U32  V3 |
	*/

	Delta	= 	U11 * (U22*U33 - U23*U32) 
		      + U12 * (U23*U31 - U21*U33) 
		      + U13 * (U21*U32 - U22*U31) ;

	/* -------------------------------------*/
	*a	= 	V1 * (U22*U33 - U23*U32) 
		      + V2 * (U23*U31 - U21*U33) 
		      + V3 * (U21*U32 - U22*U31) ;
	*a /= Delta ;

	/* -------------------------------------*/
	*b	= 	V1 * (U23*U31 - U21*U33) 
		      + V2 * (U11*U33 - U13*U31) 
		      + V3 * (U13*U21 - U11*U23) ;
	*b /= Delta ;

	/* -------------------------------------*/
	*c	= 	V1 * (U21*U32 - U22*U31) 
		      + V2 * (U12*U31 - U11*U32) 
		      + V3 * (U11*U22 - U12*U21) ;
	*c /= Delta ;
return;
}

void ki2_ah1_bh2_ch3(int n, double x[], double y[], double Dy[], 
	double (*h1)(), double (*h2)(), double (*h3)(), 
	double *a , double *b, double *c){
	/* 
	Nous supposons que les tableaux x[i] et y[i] (et les erreurs Dy[]), 
	qui contiennent des "doubles" sont donn�s. Nous utiliserons uniquement
	les premiers points i=0,...,n , car "n-1" est suppos� <= 
	aux dimensions de x[], y[] et Dy[].
	Nous cherchons ici:
		f(x) = a*h1(x) + b*h2(x) + c*h3(x)
	en utilisant les moindres carr�s PONDERES. Remarquez que h1, h2 et h3 sont
	des pointeurs �a des fonctions; les fonction seront (*h1)(x), ....

	Nous devons r�soudre le syst�me de 3 �quations � 3 inconnues (a,b,c):
		U11 a + U12 b + U13 c = V1		(U est sym�trique)
		U21 a + U22 b + U23 c = V2
		U31 a + U32 b + U33 c = V3
	avec les d�finitions suivantes :
		--------------------------------------------------------------------
		U11 = \sum w_i * h1(x_i) * h1(x) ;
		U12 = \sum w_i * h1(x_i) * h2(x) ;		
		U13 = \sum w_i * h1(x_i) * h3(x) ;

		U22 = \sum w_i * h2(x_i) * h2(x) ;
		U23 = \sum w_i * h2(x_i) * h3(x) ;

		U33 = \sum w_i * h3(x_i) * h3(x) ;

		Les autres U_ij sont :
		U21 = U12 ;
		U31 = U13 ;
		U32 = U23 ;
		--------------------------------------------------------------------
		V1 = \sum w_i * h1(x_i) * y_i ;
		V2 = \sum w_i * h2(x_i) * y_i ;
		V3 = \sum w_i * h3(x_i) * y_i ;
		--------------------------------------------------------------------
	*/
	int i;
	double x_1, y_1, w_i, h1_i , h2_i , h3_i, Delta;
	double  U11,U12,U13,V1 ; 
	double  U21,U22,U23,V2 ; 
	double  U31,U32,U33,V3 ; 
	
	U11 = 0. ;	U12 = 0. ;	U13 = 0. ;
	U23 = 0. ;	U23 = 0. ;
	U33 = 0. ;
	V1 = 0. ;	V2 = 0. ;	V3 = 0. ; /* initialisation des sommes � calculer */

	for(i=0 ; i <= n ; i++) 
		{
		if ( Dy[i] != 0. ) w_i = 1./Dy[i]/Dy[i] ;
		else {
			printf("\n\nERREUR: Division par z�ro dans le calcul des poids.\nExiting ...\n\n");
			exit(1); 
		}
		x_1 = x[i];		y_1 = y[i] ;
		h1_i = (*h1)(x_1) ;
		h2_i = (*h2)(x_1) ;
		h3_i = (*h3)(x_1) ;
		/* -------------------------------------------*/
		U11 += w_i * h1_i * h1_i ;
		U12 += w_i * h1_i * h2_i ;
		U13 += w_i * h1_i * h3_i ;

		U22 += w_i * h2_i * h2_i ;
		U23 += w_i * h2_i * h3_i ;

		U33 += w_i * h3_i * h3_i ;

		V1 += w_i * h1_i * y_1 ;
		V2 += w_i * h2_i * y_1 ;
		V3 += w_i * h3_i * y_1 ;
		}
	/* ------ On calcule maintenant les autres U_ij : U est sym�trique ---------*/
	U21 = U12 ;
	U31 = U13 ;
	U32 = U23 ;
	
	/* On r�soud maintenant pour a, b, c , avec les d�terminants */
	/* determinant_principal = Delta =
		| U11  U12  U13 |
		| U21  U22  U23 |
		| U31  U32  U33 |
		= 	U11 * (U22*U33 - U23*U32) 
		      + U12 * (U23*U31 - U21*U33) 
		      + U13 * (U21*U32 - U22*U31) 
	DONC: a = 1/Delta *
		| V1  U12  U13 |
		| V2  U22  U23 |
		| V3  U32  U33 |
	DONC: b = 1/Delta *
		| U11  V1  U13 |
		| U21  V2  U23 |
		| U31  V3  U33 |
	DONC: c = 1/Delta *
		| U11  U12  V1 |
		| U21  U22  V2 |
		| U31  U32  V3 |
	*/

	Delta	= 	U11 * (U22*U33 - U23*U32) 
		      + U12 * (U23*U31 - U21*U33) 
		      + U13 * (U21*U32 - U22*U31) ;

	/* -------------------------------------*/
	*a	= 	V1 * (U22*U33 - U23*U32) 
		      + V2 * (U23*U31 - U21*U33) 
		      + V3 * (U21*U32 - U22*U31) ;
	*a /= Delta ;

	/* -------------------------------------*/
	*b	= 	V1 * (U23*U31 - U21*U33) 
		      + V2 * (U11*U33 - U13*U31) 
		      + V3 * (U13*U21 - U11*U23) ;
	*b /= Delta ;

	/* -------------------------------------*/
	*c	= 	V1 * (U21*U32 - U22*U31) 
		      + V2 * (U12*U31 - U11*U32) 
		      + V3 * (U11*U22 - U12*U21) ;
	*c /= Delta ;
return;
}
