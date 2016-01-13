//
//  matrix.h
//  clean-merge-app
//
//  Created by Gianmarco Stinchi on 28/12/15.
//  Copyright © 2015 Gianmarco Stinchi. All rights reserved.
//


#ifndef matrix_h
#define matrix_h

#include <stdlib.h>
#include <stdbool.h>

//Definizione della struttura della matrice
typedef struct Matrix
{
    size_t w,h;			//dimensione
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
Matrix* matrix_alloc(size_t w,size_t h);
//Alloca ed inizializza da un array di double
Matrix* matrix_init(const double* values, size_t w, size_t h);
//create a rotation matrix
Matrix* matrix_rotate(double alpha);
//Alloca una matrice di identità
Matrix* matrix_identity(size_t w, size_t h);
//Alloca una matrice sqrt(in)
Matrix* matrix_sqrt(const Matrix* in);
//Alloca una matrice diagonale
Matrix* matrix_diagonal(size_t w,size_t h,double value);
//copia la matrice
Matrix* matrix_copy(const Matrix* in);
//Moltiplicazione Matrice Matrice
Matrix* matrix_multiply(const Matrix* a, const Matrix* b);
//Somma Matrice Matrice
Matrix* matrix_sum(const Matrix* a, const Matrix* b);
//Copia i valori di una matrice in un altra matrice
void matrix_put_from_sub_matrix(Matrix* dest,const Matrix* source,size_t pos[2]);
//Limita i valori
void matrix_clamp_inplace(const Matrix* a, double min, double max);
//Media moda
double matrix_mode_row(const Matrix* a,size_t row);
//Media moda
double matrix_mode_col(const Matrix* a,size_t col);
//fattorizazione di cholesky
Matrix* matrix_cholesky_factorization(const Matrix* a);
//matrix to vector
Matrix* matrix_to_vector(const Matrix* in);
//matrix to vector
Matrix* matrix_to_vector_col(const Matrix* in);
//sub matrix to vector
Matrix* matrix_sub_to_vector(const Matrix* in,size_t pos[2],size_t size[2]);
//sub matrix to vector
Matrix* matrix_sub_to_vector_col(const Matrix* in,size_t pos[2],size_t size[2]);
//row vector to matrix
Matrix* matrix_from_vector(const Matrix* in, size_t col,size_t w,size_t h);
//Moltiplicazione Matrice per scalare
Matrix* matrix_multiply_to_scalar(const Matrix* a, double scalar);
//Somma uno scalare
Matrix* matrix_sum_scalar(const Matrix* a, double scalar);
//Somma uno scalare
void matrix_sum_scalar_inplace(Matrix* a, double scalar);
//Moltiplicazione Matrice per scalare
void matrix_multiply_to_scalar_inplace(Matrix* a, double scalar);
//Imposta tutti i valori della matrice
void matrix_set(Matrix* in,size_t x,size_t y,double val);
//restituisce tutti i valori della matrice
double matrix_get(const Matrix* in,size_t x,size_t y);
//restituisce tutti i valori della matrice (sicuro)
double matrix_get_safe(const Matrix* in,size_t x,size_t y);
//get normalize value
double matrix_get_safe_norm(const Matrix* in,size_t x,size_t y);
//Libera memoria allocata
void matrix_free(Matrix* in);
//legge da file gimp
Matrix* matrix_from_gimp_file (char *filename);
//salva matrice nel formato desiderato
void matrix_save_only_buffer(const char* path, const Matrix* in);
//print matrix
void matrix_print(const Matrix* in);
void matrix_print_to_file(const Matrix* in, const char* path);

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
Matrix* matrix_column_normalized(const Matrix* in, size_t icolumn);
//normalizza il vettore dato
void matrix_column_normalized_inplace(Matrix* inout, size_t icolumn);
//Trasposta matrice
Matrix* matrix_transpose(const Matrix* in);
//Inversa matrice 2x2 
Matrix* matrix_inv2x2(const Matrix* in);
#endif /* matrix_h */
