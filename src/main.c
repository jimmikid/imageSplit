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
#include <float.h>
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
	Matrix* m = matrix_init(v, 4, 4);
	//execute cholesky factorization
	Matrix* l = matrix_cholesky_factorization(m);
	//print all
	printf("\nMatrix:\n");
	matrix_print(m);
	printf("\nCholesky factorization:\n");
	matrix_print(l);
	return 0;
}
#elif 0
int main()
{
	//values
	double v[] = {
		2, 3, 
		1, 5
	};
	//init matrix
	Matrix* m = matrix_init(v, 2, 2);
	//execute cholesky factorization
	Matrix* l = matrix_inv2x2(m);
	//print all
	printf("\nMatrix:\n");
	matrix_print(m);
	printf("\nInverse:\n");
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

Matrix* put_x_into_matrix_row(Matrix* x[2], size_t N)
{
	//Metodo che trasformi vettore di 2 matrici 
	//a singola colonna in una matrice con 2 colonne...
	Matrix* x2 = matrix_alloc(N, 2);
	//copy all elements into new Matrix
	for (size_t y = 0; y != N; ++y)
	{
		matrix_set(x2, y, 0, matrix_get(x[0], y, 0));
		matrix_set(x2, y, 1, matrix_get(x[1], y, 0));
	}
	return x2;
}

Matrix* put_x_into_matrix_col(Matrix* x[2], size_t N)
{
	//Metodo che trasformi vettore di 2 matrici 
	//a singola colonna in una matrice con 2 colonne...
	Matrix* x2 = matrix_alloc(2, N);
	//copy all elements into new Matrix
	for (size_t y = 0; y != N; ++y)
	{
		matrix_set(x2, 0, y, matrix_get(x[0], y, 0));
		matrix_set(x2, 1, y, matrix_get(x[1], y, 0));
	}
	return x2;
}

double compute_funz_ob_theta(Matrix* L, Matrix* x[2], size_t N, double theta)
{
	Matrix* U = matrix_rotate(theta);
	Matrix* A = matrix_multiply(L, U);
	//debug
#if defined( MODE_DEBUG ) && 0
	printf("U:\n");
	matrix_print(U);
	printf("A:\n");
	matrix_print(A);
#endif
	double sum[] =
	{
		matrix_get(A,0,0) + matrix_get(A,1,0),
		matrix_get(A,0,1) + matrix_get(A,1,1)
	};

	matrix_set(A, 0, 0, matrix_get(A, 0, 0) / sum[0]);
	matrix_set(A, 1, 0, matrix_get(A, 1, 0) / sum[0]);

	matrix_set(A, 0, 1, matrix_get(A, 0, 1) / sum[1]);
	matrix_set(A, 1, 1, matrix_get(A, 1, 1) / sum[1]);

	Matrix* A_t   = /* matrix_transpose */ matrix_copy(A);
	Matrix* A_t_i = matrix_inv2x2(A);

	//put 2 col. matrix into a matrix Nx2
	Matrix* x2 = put_x_into_matrix_col(x, N);

	//debug
#if defined( MODE_DEBUG ) && 0
	printf("A normalized:\n");
	matrix_print(A);
	printf("A':\n");
	matrix_print(A_t);
	printf("inv(A'):\n");
	matrix_print(A_t_i);
#endif

	//compute sk
	Matrix* sk = matrix_multiply(x2, A_t_i);
	const double s_min = 0;
	const double s_max = 255;

	for (size_t y = 0; y != N; ++y)
	{
		if (matrix_get(sk, 0, y) < s_min) matrix_set(sk, 0, y, s_min);
		if (matrix_get(sk, 1, y) < s_min) matrix_set(sk, 1, y, s_min);
		/*
		if (matrix_get(sk, 0, y) > s_max) matrix_set(sk, 0, y, s_max);
		if (matrix_get(sk, 1, y) > s_max) matrix_set(sk, 1, y, s_max);
		*/
	}
#if defined( MODE_DEBUG ) && 0
	printf("SK:\n");
	matrix_print(sk);
#endif
	//compute output
	double n = 0;
	//sum all components
	for (size_t y = 0; y != N; ++y)
	{
		n += abs(matrix_get(sk, 0, y) * matrix_get(sk, 1, y));
	}

	//dealloc all
	matrix_free(U);
	matrix_free(A);
	matrix_free(x2);
	matrix_free(sk);
	matrix_free(A_t);
	matrix_free(A_t_i);
	//return computed value
	return n;
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

	double p = 0.1;
	double nn = DBL_MAX;
	double funz_ob[] = { 0.0, 0.0, 0.0, 0.0 };
	double theta2 = theta;
	double theta_v[] = { theta2 - p, theta2 , theta2 + p };
	//...
	for (size_t i = 0; i != 3; ++i)
	{
		funz_ob[i] = compute_funz_ob_theta(L, x, nm2, theta_v[i]);
	}
	//loop count
	int loop_count = 1;
	//..
	while (true)
	{
		nn = funz_ob[1];

		if (funz_ob[1] < funz_ob[0] && funz_ob[1] < funz_ob[2])
		{
			break;
		}
		else
		{
			if (funz_ob[0] < funz_ob[2])
			{
				//move to left
				theta_v[0] = theta_v[0] - p;
				theta_v[1] = theta_v[0] ;
				theta_v[2] = theta_v[1] ;
				//recompute
				funz_ob[0] = compute_funz_ob_theta(L, x, nm2, theta_v[0]);
				funz_ob[1] = funz_ob[0];
				funz_ob[2] = funz_ob[1];
			}
			else if (funz_ob[2] < funz_ob[0])
			{
				//move to right
				theta_v[0] = theta_v[1];
				theta_v[1] = theta_v[2];
				theta_v[2] = theta_v[2] + p;
				//recompute
				funz_ob[0] = funz_ob[1];
				funz_ob[1] = funz_ob[2];
				funz_ob[2] = compute_funz_ob_theta(L, x, nm2, theta_v[2]);
			}
		}

		++loop_count;
	};

	const double r  = (sqrt(5.) - 1.) / 2.;
	const double r1 = 1 - r;

#if 0
	//x(1) == x[0](0,0);
	//x(2) == x[0](1,0);
	//x(3) == x[0](2,0);
	//x(4) == x[0](3,0);
	matrix_set(x[0], 0, 0, theta_v[0]);
	matrix_set(x[0], 3, 0, theta_v[2]);

	if (fabs(theta_v[2] - theta_v[1]) > fabs(theta_v[1] - theta_v[0]))
	{
		matrix_set(x[0], 1, 0, theta_v[1]);
		matrix_set(x[0], 2, 0, theta_v[1] + r1*(theta_v[2]-theta_v[1]));
	}
	else
	{
		matrix_set(x[0], 2, 0, theta_v[1]);
		matrix_set(x[0], 1, 0, theta_v[1] - r1*(theta_v[1] - theta_v[0]));
	}

	funz_ob[1] = compute_funz_ob_theta(L, x, nm2, matrix_get(x[0], 1, 0));
	funz_ob[3] = compute_funz_ob_theta(L, x, nm2, matrix_get(x[0], 3, 0));
	funz_ob[2] = compute_funz_ob_theta(L, x, nm2, matrix_get(x[0], 2, 0));

	while (fabs(matrix_get(x[0], 2, 0)- matrix_get(x[0], 1, 0)) > 1.0e-9)
	{
		if (funz_ob[2] < funz_ob[1])
		{
			matrix_set(x[0], 0, 0, matrix_get(x[0], 1, 0));
			matrix_set(x[0], 1, 0, matrix_get(x[0], 2, 0));
			matrix_set(x[0], 2, 0, r*matrix_get(x[0], 1, 0)+r1*matrix_get(x[0], 3, 0));
			funz_ob[0] = funz_ob[1];
			funz_ob[1] = funz_ob[2];
			funz_ob[2] = compute_funz_ob_theta(L, x, nm2, matrix_get(x[0], 2, 0));
		}
		else
		{
			matrix_set(x[0], 3, 0, matrix_get(x[0], 2, 0));
			matrix_set(x[0], 2, 0, matrix_get(x[0], 1, 0));
			matrix_set(x[0], 1, 0, r*matrix_get(x[0], 2, 0) + r1*matrix_get(x[0], 0, 0));
			funz_ob[3] = funz_ob[2];
			funz_ob[2] = funz_ob[1];
			funz_ob[1] = compute_funz_ob_theta(L, x, nm2, matrix_get(x[0], 1, 0));
		}
	}

	if ( funz_ob[1] < funz_ob[2] )
	{
		output.theta = matrix_get(x[0], 1, 0);
		output.estimate = funz_ob[1];
	}
	else
	{
		output.theta    = matrix_get(x[0], 2, 0);
		output.estimate = funz_ob[2];
	}
#else
	//ipotizzando che sia stata programmata da 
	//una scimmia che non sa programmare...
	double tmp[4] = { 0.0, 0.0, 0.0, 0.0 };

	tmp[0] = theta_v[0];
	tmp[3] = theta_v[2];

	if (fabs(theta_v[2] - theta_v[1]) > fabs(theta_v[1] - theta_v[0]))
	{
		tmp[1] = theta_v[1];
		tmp[2] = theta_v[1] + +r1*(theta_v[2] - theta_v[1]);
	}
	else
	{
		tmp[2] = theta_v[1];
		tmp[1] = theta_v[1] - r1*(theta_v[1] - theta_v[0]);
	}

	funz_ob[1] = compute_funz_ob_theta(L, x, nm2, tmp[1]);
	funz_ob[3] = compute_funz_ob_theta(L, x, nm2, tmp[3]);
	funz_ob[2] = compute_funz_ob_theta(L, x, nm2, tmp[2]);

	while (fabs(tmp[2] - tmp[1]) > 1.0e-10)
	{

		if (funz_ob[2] < funz_ob[1])
		{
			tmp[0] = tmp[1];
			tmp[1] = tmp[2];
			tmp[2] = r*tmp[1] + r1*tmp[3];
			funz_ob[0] = funz_ob[1];
			funz_ob[1] = funz_ob[2];
			funz_ob[2] = compute_funz_ob_theta(L, x, nm2, tmp[2]);
		}
		else
		{
			tmp[3] = tmp[2];
			tmp[2] = tmp[1];
			tmp[1] = r*tmp[2] + r1*tmp[0];
			funz_ob[3] = funz_ob[2];
			funz_ob[2] = funz_ob[1];
			funz_ob[1] = compute_funz_ob_theta(L, x, nm2, tmp[1]);
		}
	}

	if (funz_ob[1] < funz_ob[2])
	{
		output.theta    = tmp[1];
		output.estimate = funz_ob[1];
	}
	else
	{
		output.theta    = tmp[2];
		output.estimate = funz_ob[2];
	}
#endif

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
#if defined( MODE_DEBUG ) && 0
	printf("new_L:\n");
	matrix_print(new_L);
#endif
	//compute sum1/2
	double sum1 = matrix_get(L, 0, 0) + matrix_get(L, 1, 0);
	double sum2 = matrix_get(L, 0, 1) + matrix_get(L, 1, 1);
	//compute estimate
	double A_estimate_init[] =
	{
		matrix_get(L, 0, 0)/sum1, matrix_get(L, 1, 0)/sum1,
		matrix_get(L, 0, 1)/sum2, matrix_get(L, 1, 1)/sum2
	};
	Matrix* A_estimate = matrix_init(A_estimate_init,2, 2);
#if defined( MODE_DEBUG ) && 0
	printf("A_estimate:\n");
	matrix_print(A_estimate);
#endif
	//dealloc all
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

ImageMarge split_images(ImageMarge images,  size_t nm[2], size_t nm2)
{
	//front/back to vector
	Matrix* front_row = matrix_to_vector(images.front);
	Matrix* back_row = matrix_to_vector(images.back);
	//debug
#if defined( MODE_DEBUG ) && 0
	printf("Front vector:\n");
	matrix_print(front_row);
	printf("Back vector:\n");
	matrix_print(back_row);
#endif
    //inv front and back
    Matrix* front_inv_row =r_matrix_inverse(front_row);
    Matrix* back_inv_row  =r_matrix_inverse(back_row);
	//xorig link to front_inv and back_inv
	Matrix* x_orig[] =
	{
		front_inv_row,
		back_inv_row
	};
#if defined( MODE_DEBUG ) && 0
	printf("Front inverse vector:\n");
	matrix_print(front_inv_row);
	printf("Back inverse vector:\n");
	matrix_print(back_inv_row);
#endif
    //mode
    double mode_front = matrix_mode_row(front_inv_row, 0);
    double mode_back  = matrix_mode_row(back_inv_row, 0);
	//mode is the zero of all values
	Matrix* x[] =
	{
		matrix_alloc(front_inv_row->w,1),
		matrix_alloc(back_inv_row->w,1)
	};
	//..
	for (size_t i = 0; i != front_inv_row->w; ++i)
	{
		if (matrix_get(front_inv_row, i, 0) < mode_front)
			matrix_set(x[0], i, 0, 0);
		else
			matrix_set(x[0], i, 0, matrix_get(front_inv_row, i, 0) - mode_front);

		if (matrix_get(back_inv_row, i, 0) < mode_back)
			matrix_set(x[1], i, 0, 0);
		else
			matrix_set(x[1], i, 0, matrix_get(back_inv_row, i, 0) - mode_back);
	}

#if defined( MODE_DEBUG ) && 0
	printf("Front vector - mode:\n");
	matrix_print(x[0]);
	printf("Back vector - mode:\n");
	matrix_print(x[1]);
#endif
	//Processing
	Matrix* A_stimate     = estimate(x, nm2);
	Matrix* A_stimate_i   = matrix_inv2x2(A_stimate);
	Matrix* A_stimate_i_t = /* matrix_transpose */ matrix_copy(A_stimate_i);
	Matrix* x_orig2		  = put_x_into_matrix_col(x_orig, x_orig[0]->w);
	Matrix* sk			  = matrix_multiply(x_orig2, A_stimate_i_t);
	//get mode of images
	double	new_mode_front = matrix_mode_row(sk, 0);
	double	new_mode_back  = matrix_mode_row(sk, 1);
#if defined( MODE_DEBUG ) && 0
	printf("new_mode_front: %g, new_mode_back:%g\n", new_mode_front, new_mode_back);
	//printf("x_orig2\n");
	//matrix_print(x_orig2);
	printf("A_stimate\n");
	matrix_print(A_stimate);
	printf("A_stimate_i\n");
	matrix_print(A_stimate_i);
	printf("A_stimate_i_t\n");
	matrix_print(A_stimate_i_t);
#endif
	//create matrix to apply to the images matrix
	for (int y = 0; y != x_orig2->h; ++y)
	{
		matrix_set(sk, 0, y, matrix_get(sk, 0, y) + mode_front - new_mode_front);
		matrix_set(sk, 1, y, matrix_get(sk, 1, y) + mode_back - new_mode_back);
	}
	//inverso colore
	Matrix* ss1_inv = matrix_from_vector(sk, 0, nm[0], nm[1]);
	Matrix* ss2_inv = matrix_from_vector(sk, 1, nm[0], nm[1]);
	Matrix* ss1 = r_matrix_inverse(ss1_inv);
	Matrix* ss2 = r_matrix_inverse(ss2_inv);

	ImageMarge output;
	output.front = ss1;
	output.back  = ss2;
	//free all
	matrix_free(ss1_inv);
	matrix_free(ss2_inv);
	matrix_free(sk);
	matrix_free(A_stimate_i_t);
	matrix_free(A_stimate_i);
	matrix_free(A_stimate);
	matrix_free(x[0]);
	matrix_free(x[1]);
	matrix_free(front_inv_row);
	matrix_free(back_inv_row);
	matrix_free(front_row);
	matrix_free(back_row);
	//return
	return output;
#if 0
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
	rgb_matrix_init(matrix_copy(front_inv_mode),
                    matrix_copy(front_inv_mode),
                    matrix_copy(front_inv_mode));
    rgb_matrix_to_tga_file("__front_zero.tga", front_image);
    
    RGB_Matrix back_image =
	rgb_matrix_init(matrix_copy(back_inv_mode),
                    matrix_copy(back_inv_mode),
                    matrix_copy(back_inv_mode));
    rgb_matrix_to_tga_file("__back_zero.tga", back_image);
	#endif
#endif
}

const int k = 100; //step for the grid
const int d1 = 5;  //minimun distance from mode
const int d2 = 5;  //minimun distance from original image
const int c = 1;   //correction factor


void correct_rgb(RGB_Matrix ss,
				 RGB_Matrix orig,
				 size_t pos[2],
				 size_t sub_nm[2],
				 size_t nm[2],
				 size_t nm2)
{
	/*compute mode of input image for the three channels*/
	//mode(red)
	//mode(green)
	//mode(blue)
	Matrix* rgb_rows[3] =
	{
		matrix_sub_to_vector(ss.matrix_r, pos, sub_nm),
		matrix_sub_to_vector(ss.matrix_g, pos, sub_nm),
		matrix_sub_to_vector(ss.matrix_b, pos, sub_nm)
	};

	//matrix_print_to_file(rgb_rows[0],"pisello.txt");

	double mode[3] =
	{
		matrix_mode_row(rgb_rows[0],0),
		matrix_mode_row(rgb_rows[1],0),
		matrix_mode_row(rgb_rows[2],0)
	};

	for (size_t i = 0; i != 3; ++i)
	{
		matrix_free(rgb_rows[i]);
	}

	double mode_md = (mode[0] + mode[1] + mode[2]) / 3.0;

	for (size_t sub_i = 0; sub_i != sub_nm[0]; ++sub_i)
	for (size_t sub_j = 0; sub_j != sub_nm[1]; ++sub_j)
	{
		size_t i = sub_i + pos[0];
		size_t j = sub_j + pos[1];

		if (fabs((matrix_get(ss.matrix_r, i, j) + 
				  matrix_get(ss.matrix_g, i, j) +
				  matrix_get(ss.matrix_b, i, j)) / 3.0 - mode_md) > d1 &&
			fabs(matrix_get(ss.matrix_r, i, j) - matrix_get(orig.matrix_r, i, j)) > d2 &&
			fabs(matrix_get(ss.matrix_g, i, j) - matrix_get(orig.matrix_g, i, j)) > d2 &&
			fabs(matrix_get(ss.matrix_b, i, j) - matrix_get(orig.matrix_b, i, j)) > d2)
		{
#if 1
			matrix_set(ss.matrix_r, i, j, mode[0] + (c*(matrix_get(ss.matrix_r, i, j) - mode[0]) / matrix_get(ss.matrix_r, i, j)));
			matrix_set(ss.matrix_g, i, j, mode[1] + (c*(matrix_get(ss.matrix_g, i, j) - mode[1]) / matrix_get(ss.matrix_g, i, j)));
			matrix_set(ss.matrix_b, i, j, mode[2] + (c*(matrix_get(ss.matrix_b, i, j) - mode[2]) / matrix_get(ss.matrix_b, i, j)));
#else
			matrix_set(ss.matrix_r, i, j, matrix_get(ss.matrix_r, i, j));
			matrix_set(ss.matrix_g, i, j, matrix_get(ss.matrix_g, i, j));
			matrix_set(ss.matrix_b, i, j, matrix_get(ss.matrix_b, i, j));
#endif
		}
	}

}
void correct_color(RGB_Matrix ss, 
				   RGB_Matrix orig,
	               size_t nm[2], 
	               size_t nm2)
{


	if (!k)
	{
		size_t pos[] = { 0,0 };
		correct_rgb(ss, orig, pos, nm, nm, nm2);
	}
	
	size_t i = 0; 
	size_t j = 0;

	for ( i = 0; i < nm[0] - k; i += k)
	{
		for ( j = 0; j < nm[1] - k; j += k)
		{
			size_t pos[2]  = { i, j };
			size_t sub_nm[2] = { k, k };
			correct_rgb(ss,orig,pos,sub_nm,nm,nm2);
		}
		size_t pos[2] = { i, j };
		size_t sub_nm[2] = { k, nm[1] - j };
		correct_rgb(ss, orig, pos, sub_nm, nm, nm2);
	}
	for ( j = 0; j < nm[1] - k; j += k )
	{
		size_t pos[2] = { i, j };
		size_t sub_nm[2] = { nm[0] - i, k };
		correct_rgb(ss, orig, pos, sub_nm, nm, nm2);
	}
	size_t pos[2] = { i, j };
	size_t sub_nm[2] = { nm[0]-i, nm[1]-j };
	correct_rgb(ss, orig, pos, sub_nm, nm, nm2);

}

int main(int argc,const char* argv[])
{
#if 0
    //paths
    const char* path_image[] =
    {
        "assets/mistura_recto.tga",
        "assets/mistura_verso.tga"
    };
    //load images
	RGB_Matrix image[2];
	image[0] = rgb_matrix_from_tga_file(path_image[0]);
	image[1] = rgb_matrix_from_tga_file(path_image[1]);
    //assets
    assert(image[0].matrix_r->w == image[1].matrix_r->w &&
           image[0].matrix_r->h == image[1].matrix_r->h );
    //sizes
    size_t size[]=
    {
        image[0].matrix_r->w,
        image[0].matrix_r->h
    };
	////////////////////////////////////////////////////////////////////////////////////
    //elements count
    size_t elements = size[0] * size[1];
	//marge images?
	bool b_not_merge = true;
	//call split images (red channel)
	Matrix*     colors_r[] = { image[0].matrix_r, image[1].matrix_r };
	ImageMarge  images_r = b_not_merge ? marge_images_init(colors_r) : marge_images(colors_r);
				images_r = split_images(images_r, size, elements);
	//call split images (green channel)
    Matrix*     colors_g[] = { image[0].matrix_g, image[1].matrix_g };
	ImageMarge  images_g = b_not_merge ? marge_images_init(colors_g) : marge_images(colors_g);
				images_g = split_images(images_g, size, elements);
	//call split images (blue channel)
	Matrix*     colors_b[] = { image[0].matrix_b, image[1].matrix_b };
	ImageMarge  images_b = b_not_merge ? marge_images_init(colors_b) : marge_images(colors_b);
				images_b = split_images(images_b, size, elements);

	////////////////////////////////////////////////////////////////////////////////////
	RGB_Matrix front_image = rgb_matrix_init(images_r.front, images_g.front, images_b.front);
	//correct_color(front_image, image[0], size,  elements);
	RGB_Matrix back_image = rgb_matrix_init(images_r.back, images_g.back, images_b.back);
	//correct_color(back_image, image[1], size, elements);
	////////////////////////////////////////////////////////////////////////////////////
	rgb_matrix_to_tga_file("out_front.tga", front_image);
	rgb_matrix_to_tga_file("out_back.tga", back_image);
#else
	//paths
	const char* path_image[] =
	{
		"assets/mistura_recto.tga",
		"assets/mistura_verso.tga",

		"assets/verso_stimato.tga",
		"assets/recto_stimato.tga"
	};
	//load images
	RGB_Matrix image[2];
	image[0] = rgb_matrix_from_tga_file(path_image[0]);
	image[1] = rgb_matrix_from_tga_file(path_image[1]);
	
	size_t size[] =
	{
		image[0].matrix_r->w,
		image[0].matrix_r->h
	};
	//elements count
	size_t elements = size[0] * size[1];

	RGB_Matrix front_image = rgb_matrix_from_tga_file(path_image[2]);
	correct_color(front_image, image[0], size, elements);

	RGB_Matrix back_image = rgb_matrix_from_tga_file(path_image[3]);
	correct_color(back_image, image[1], size, elements);

	rgb_matrix_to_tga_file("_out_front.tga", front_image);
	rgb_matrix_to_tga_file("_out_back.tga", back_image);
#endif
    return 0;
}

#elif 0
double deleteLT128(double value, size_t x, size_t y,Matrix* in)
{
    return value < 128 ?  0 : value;
}
double deleteGT128(double value, size_t x, size_t y,Matrix* in)
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