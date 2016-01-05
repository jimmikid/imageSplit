//
//  leggi.c
//
//
//  Created by Gianmarco Stinchi on 07/09/15.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <matrix.h>
#include <assert.h>
#include <math.h>
#include <tga_loader.h>
#define MODE_DEBUG

#if 0

int main()
{
    
    Matrix* A = matrix_alloc(5, 1);
    matrix_set(A,0,0,0.7);
    matrix_set(A,1,0,0.6);
    matrix_set(A,2,0,0.1);
    matrix_set(A,3,0,0.2);
    matrix_set(A,4,0,0.3);
    matrix_print(A);
    //get mode 
    printf("moda: %f\n",matrix_mode_row(A, 0));
    return 0;
}

#elif 0
int main()
{
    //values
    double v[] = {
                   18, 22,  54,  42,
                   22, 70,  86,  62,
                   54, 86, 174, 134,
                   42, 62, 134, 106
                 };
    //init matrix
    Matrix* m=matrix_init(v, 4, 4);
    //execute cholesky factorization
    Matrix* l=matrix_cholesky_factorization(m);
    //print all
    printf("\nMatrix:\n");
    matrix_print(m);
    printf("\nCholesky factorization:\n");
    matrix_print(l);
    return 0;
}
#elif 1

typedef struct ImageMarge
{
	Matrix* front;
	Matrix* back;
}
ImageMarge;



ImageMarge marge_images_init(Matrix* s[2])
{
	ImageMarge output;
	output.front = s[0];
	output.back = s[1];
	return output;
}

ImageMarge marge_images(Matrix* s[2])
{
	//alloc output function
	ImageMarge output;

	//merge images (front)
	Matrix* front_0 = matrix_multiply_to_scalar(s[0], 0.7);
	Matrix* front_1 = matrix_multiply_to_scalar(s[1], 0.3);
	output.front    = matrix_sum(front_0, front_1);
	matrix_free(front_0);
	matrix_free(front_1);
	//merge images (back)
	Matrix* back_0 = matrix_multiply_to_scalar(s[1], 0.8);
	Matrix* back_1 = matrix_multiply_to_scalar(s[0], 0.2);
	output.back    = matrix_sum(back_0, back_1);
	matrix_free(back_0);
	matrix_free(back_1);

	return output;
}

typedef struct EstimateThetaReturn
{ 
	double theta;
	double estimate;
}
EstimateThetaReturn;

EstimateThetaReturn estimate_theta_funz(Matrix* L, Matrix* x[2], size_t nm2, double theta)
{
	EstimateThetaReturn output;
	return output;
}

Matrix* estimate(Matrix* x[2], size_t nm2)
{
	//default angle
	const double theta = M_PI_2;
	//compute C
	Matrix* t_x[] =
	{
		matrix_transpose(x[0]),
		matrix_transpose(x[1])
	};
	//C block version
	Matrix* C_v[] =
	{
		matrix_multiply(x[0], t_x[0]), matrix_multiply(x[0], t_x[1]),
		matrix_multiply(x[1], t_x[0]), matrix_multiply(x[1], t_x[1]) 
	};
	//create c init
	double C_init[]=
	{
		matrix_get(C_v[0],0,0),matrix_get(C_v[1],0,0),
		matrix_get(C_v[2],0,0),matrix_get(C_v[3],0,0)
	};
	Matrix* C = matrix_init(C_init,2, 2);
	//get eigs
	Eigen2x2 VE = matrix_eigen2x2(C);
	//compute L = ( V * sqrt(E) ) * V'
	Matrix* V_t     = matrix_transpose(VE.V);
	Matrix* E_s     = matrix_sqrt(VE.D);
	Matrix* V_E_s   = matrix_multiply(VE.V, E_s);
	Matrix* L       = matrix_multiply(V_E_s, V_t);
	//compute theta
	EstimateThetaReturn estimate_t = estimate_theta_funz(L, x, nm2, theta);
	//compute U
	Matrix* U     = matrix_rotate(estimate_t.theta);
	Matrix* new_L = matrix_multiply(L, U);
	//compute sum1/2
	double sum1 = matrix_get(L, 1, 1) + matrix_get(L, 2, 1);
	double sum2 = matrix_get(L, 1, 2) + matrix_get(L, 2, 2);
	//compute estimate
	double A_estimate_init[] =
	{
		matrix_get(L, 1, 1)/sum1, matrix_get(L, 2, 1)/sum1,
		matrix_get(L, 1, 2)/sum2, matrix_get(L, 2, 2)/sum2
	};
	Matrix* A_estimate = matrix_init(A_estimate_init,2, 2);
	//dealoc all
	int i = 0;
	for (i = 0; i != 2; ++i) matrix_free(t_x[i]);
	for (i = 0; i != 4; ++i) matrix_free(C_v[i]);
	matrix_free(C);
	eigen2x2_free(VE);
	matrix_free(V_t);
	matrix_free(E_s);
	matrix_free(V_E_s);
	matrix_free(L);
	matrix_free(U);
	matrix_free(new_L);
	//return matrix
	return A_estimate;
}

void split_images(ImageMarge images,  size_t nm[2], size_t nm2)
{
	//front/back to vector
	Matrix* front_row = matrix_to_vector(images.front);
	Matrix* back_row = matrix_to_vector(images.back);
    //inv front and back
    Matrix* front_inv=r_matrix_inverse(images.front);
    Matrix* back_inv =r_matrix_inverse(images.back);
    //front/back to vector
    Matrix* front_inv_row = matrix_to_vector(front_inv);
    Matrix* back_inv_row  = matrix_to_vector(back_inv);
    //mode
    double mode_front = matrix_mode_row(front_inv_row, 0);
    double mode_back  = matrix_mode_row(back_inv_row, 0);
	//mode is the zero of all values
	Matrix* x[]=
	{
		matrix_alloc(front_inv_row->w,1),
		matrix_alloc(back_inv_row->w,1)
	};
	//..
	for (int i = 0; i != front_inv_row->w; ++i)
	{
		if (matrix_get(front_inv_row, x, 0) < mode_front)
			matrix_set(x[0], x, 0, 0);
		else
			matrix_set(x[0], x, 0, matrix_get(x[0], x, 0) - mode_front);

		if (matrix_get(back_inv_row, x, 0) < mode_back)
			matrix_set(x[1], x, 0, 0);
		else
			matrix_set(x[1], x, 0, matrix_get(x[1], x, 0) - mode_back);
	}
	//Processing

	//matrix
	Matrix* A = matrix_alloc(2, 2);
	matrix_set(A, 0, 0, 0.7);    matrix_set(A, 1, 0, 0.3);
	matrix_set(A, 0, 1, 0.2);    matrix_set(A, 1, 1, 0.8);

    //to zero
    Matrix* front_inv_row_mode = matrix_sum_scalar(front_inv_row, -mode_front);
    Matrix* back_inv_row_mode = matrix_sum_scalar(back_inv_row, -mode_back);
    //rebuild images
    Matrix* front_inv_mode = matrix_from_vector(front_inv_row_mode, 0, (int)nm[0], (int)nm[1]);
    Matrix* back_inv_mode = matrix_from_vector(back_inv_row_mode, 0,   (int)nm[0], (int)nm[1]);
    //debug
#ifdef MODE_DEBUG
    RGB_Matrix front_image =
    RGB_Matrix_init(matrix_copy(front_inv_mode),
                    matrix_copy(front_inv_mode),
                    matrix_copy(front_inv_mode));
    rgb_matrix_to_tga_file("__front_zero.tga", front_image);
    
    RGB_Matrix back_image =
    RGB_Matrix_init(matrix_copy(back_inv_mode),
                    matrix_copy(back_inv_mode),
                    matrix_copy(back_inv_mode));
    rgb_matrix_to_tga_file("__back_zero.tga", back_image);
#endif
}

int main(int argc,const char* argv[])
{
    //paths
    const char* path_image[] =
    {
        "assets/img11.tga",
        "assets/img12.tga"
    };
    //load images
    RGB_Matrix image[] =
    {
        rgb_matrix_from_tga_file(path_image[0]),
        rgb_matrix_from_tga_file(path_image[1])
    };
    //assets
    assert(image[0].matrix_r->w == image[0].matrix_r->h &&
           image[0].matrix_r->w == image[1].matrix_r->w &&
           image[0].matrix_r->h == image[1].matrix_r->h);
    //sizes
    size_t size[]=
    {
        image[0].matrix_r->w,
        image[0].matrix_r->h
    };
    //elements count
    size_t elements = size[0] * size[1];
    //make arguments
    Matrix* colors[] = { image[0].matrix_r, image[1].matrix_r };
	ImageMarge  images = false ? marge_images_init(colors) : marge_images(colors);
	//call split images (red channel)
    split_images(images, size, elements);
    return 0;
}

#elif 0
double deleteLT128(double value, int x, int y,Matrix* in)
{
    return value < 128 ?  0 : value;
}
double deleteGT128(double value, int x, int y,Matrix* in)
{
    return value > 128 ?  0 : value;
}

int main()
{
    printf("Inserisci nome file: ");
    char filename[64];
    scanf("%s", filename);
    Matrix* mat = matrix_from_gimp_file(filename);
    if(mat)
    {
        matrix_save_only_buffer("gerace_prova.pgm",mat);
        matrix_free(mat);
    }
    else
        printf("errore caricamento file. ciao!\n");

    Matrix *A = matrix_alloc(2, 2);
    matrix_set(A, 0, 0, -5); matrix_set(A, 1, 0, -2);
    matrix_set(A, 0, 1, -2); matrix_set(A, 1, 1, 3);
    
    Eigenvalues2x2 eigen = matrix_eigenvalues2x2(A);
    Matrix *E = matrix_alloc(2, 2);
    matrix_set(E, 0, 0, eigen.lambda1); matrix_set(E, 1, 0, 0);
    matrix_set(E, 0, 1, 0); matrix_set(E, 1, 1, eigen.lambda2);
    
    Eigenvctors2x2 vectors = matrix_eigenvectors2x2(A);
    
    Matrix* V = matrix_alloc(2, 2);
    
    if(vectors.success)
    {
        Matrix* v1n = matrix_column_normalized(vectors.v1,0);
        Matrix* v2n = matrix_column_normalized(vectors.v2,0);
        matrix_set(V, 0, 0, matrix_get(v1n, 0, 0)); matrix_set(V, 1, 0, matrix_get(v2n, 0, 0));
        matrix_set(V, 0, 1, matrix_get(v1n, 0, 1)); matrix_set(V, 1, 1, matrix_get(v2n, 0, 1));
        //dealloc
        matrix_free(v1n);
        matrix_free(v2n);
    }
    else printf("autovettori non trovati");
    
    Matrix* Vt = matrix_transpose(V);
    Matrix* VE = matrix_multiply(V, E);
    Matrix* VEVt = matrix_multiply(VE, Vt);
    
    printf("A:\n");
    matrix_print(A);
    
    printf("VEVt:\n");
    matrix_print(VEVt);
    //dealloc all
    eigenvctors2x2_free(vectors);
    matrix_free(A);
    matrix_free(E);
    matrix_free(V);
    matrix_free(Vt);
    matrix_free(VEVt);
    
    return 0;
}
#endif