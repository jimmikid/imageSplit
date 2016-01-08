//
//  matrix.h
//  ImageProcessors
//
//  Created by Gianmarco Stinchi on 02/10/15.
//  Copyright © 2015 Gianmarco Stinchi. All rights reserved.
//

#ifndef matrix_h
#define matrix_h

#include <stdbool.h>

//Definizione della struttura della matrice
typedef struct Matrix
{
    int w,h;			//dimensione
    double **buffer;	//contiene i numeri della matrice
}
Matrix;

//autovalori matrice 2x2
typedef struct Eigenvalues2x2
{
    bool   success;
    double lambda1;
    double lambda2;
}
Eigenvalues2x2;

//autovettori matrice 2x2
typedef struct Eigenvectors2x2
{
    bool    success;
    Matrix* v1;
    Matrix* v2;
}
Eigenvctors2x2;

//v/d eig, like mathlab
typedef struct Eigen2x2
{
	bool    success;
	Matrix* V;		//Matrice degli Autovettori
	Matrix* D;		//Matrice diagonale contenente gli Autovalori
}
Eigen2x2;

//Allocazione dinamica del tipo matrice
Matrix* matrix_alloc(int w,int h);
//Alloca ed inizializza da un array di double
Matrix* matrix_init(double* values,int w,int h);
//create a rotation matrix
Matrix* matrix_rotate(double alpha);
//Alloca una matrice di identità
Matrix* matrix_identity(int w, int h);
//Alloca una matrice sqrt(in)
Matrix* matrix_sqrt(const Matrix* in);
//Alloca una matrice diagonale
Matrix* matrix_diagonal(int w,int h,double value);
//copia la matrice
Matrix* matrix_copy(const Matrix* in);
//Moltiplicazione Matrice Matrice
Matrix* matrix_multiply(const Matrix* a, const Matrix* b);
//Somma Matrice Matrice
Matrix* matrix_sum(const Matrix* a, const Matrix* b);
//Media moda
double matrix_mode_row(const Matrix* a,int row);
//fattorizazione di cholesky
Matrix* matrix_cholesky_factorization(const Matrix* a);
//matrix to vector
Matrix* matrix_to_vector(const Matrix* in);
//row vector to matrix
Matrix* matrix_from_vector(const Matrix* in, int row,int w,int h);
//Moltiplicazione Matrice per scalare
Matrix* matrix_multiply_to_scalar(const Matrix* a, double scalar);
//Somma uno scalare
Matrix* matrix_sum_scalar(const Matrix* a, double scalar);
//Somma uno scalare
void matrix_sum_scalar_inplace(Matrix* a, double scalar);
//Moltiplicazione Matrice per scalare
void matrix_multiply_to_scalar_inplace(Matrix* a, double scalar);
//Imposta tutti i valori della matrice
void matrix_set(Matrix* in,int x,int y,double val);
//restituisce tutti i valori della matrice
double matrix_get(const Matrix* in,int x,int y);
//restituisce tutti i valori della matrice (sicuro)
double matrix_get_safe(const Matrix* in,int x,int y);
//get normalize value
double matrix_get_safe_norm(const Matrix* in,int x,int y);
//Libera memoria allocata
void matrix_free(Matrix* in);
//legge da file gimp
Matrix* matrix_from_gimp_file (char *filename);
//salva matrice nel formato desiderato
void matrix_save_only_buffer(const char* path, const Matrix* in);
//print matrix
void matrix_print(const Matrix* in);

//calcolo eigenvalues (autovalori) 2x2
Eigenvalues2x2 matrix_eigenvalues2x2(Matrix* in);

//calcolo degli autovettori su matrice 2x2
Eigenvctors2x2 matrix_eigenvectors2x2(Matrix* in);

//calcolo degli autovettori su matrice 2x2 e li normalizza
Eigenvctors2x2 matrix_eigenvectors2x2_normalized(Matrix* in);

//calcolo gli autovettori (matrice 2x2) e gli auto valori (matrice 2*2, diagonale)
Eigen2x2 matrix_eigen2x2(Matrix* in);

//Libera memoria allocata
void eigenvectors2x2_free(Eigenvctors2x2 values);
//Libera memoria allocata
void eigen2x2_free(Eigen2x2 values);
//Normalizzazione vettore
Matrix* matrix_column_normalized(const Matrix* in, int icolumn);
//normalizza il vettore dato
void matrix_column_normalized_inplace(Matrix* inout, int icolumn);
//Trasposta matrice
Matrix* matrix_transpose(const Matrix* in);
//Inversa matrice 2x2 
Matrix* matrix_inv2x2(const Matrix* in);
#endif /* matrix_h */
