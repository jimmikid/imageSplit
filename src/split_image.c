//
//  split_image.c
//  ImageProcessors
//
//  Created by Gianmarco Stinchi on 02/10/15.
//  Copyright © 2015 Gianmarco Stinchi. All rights reserved.
//
#include <stdlib.h>
#include <split_image.h>
#include <tga_loader.h>
#include <matrix.h>
#include <math.h>
#include <float.h>
#include <string.h>

#define MODE_DEBUG

static long s_loop_limit = 10000;

void set_loop_limit(long loop_limit)
{
    s_loop_limit = loop_limit;
}

long get_loop_limit()
{
    return s_loop_limit;
}

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
    //copy all
    memcpy(x2->buffer[0],x[0]->buffer[0],sizeof(double)*N);
    memcpy(x2->buffer[1],x[1]->buffer[0],sizeof(double)*N);
    //Return
    return x2;
}

double compute_funz_ob_theta(Matrix* L, Matrix* x[2], size_t N, double theta)
{
    Matrix* U = matrix_rotate(theta);
    Matrix* A = matrix_multiply(L, U);
    //debug
    double sum[] =
    {
        matrix_get(A,0,0) + matrix_get(A,1,0),
        matrix_get(A,0,1) + matrix_get(A,1,1)
    };
    
    matrix_set(A, 0, 0, matrix_get(A, 0, 0) / sum[0]);
    matrix_set(A, 1, 0, matrix_get(A, 1, 0) / sum[0]);
    
    matrix_set(A, 0, 1, matrix_get(A, 0, 1) / sum[1]);
    matrix_set(A, 1, 1, matrix_get(A, 1, 1) / sum[1]);
    
    Matrix* A_t   = matrix_transpose(A);
    Matrix* A_t_i = matrix_inv2x2(A_t);
    //put 2 col. matrix into a matrix Nx2
    Matrix* x2 = put_x_into_matrix_col(x, N);
    //compute sk
    Matrix* sk = matrix_multiply(x2, A_t_i);
    
    const double s_min = 0;
    const double s_max = 255;
    
    for (size_t y = 0; y != N; ++y)
    {
        if (matrix_get(sk, 0, y) < s_min) matrix_set(sk, 0, y, s_min);
        if (matrix_get(sk, 1, y) < s_min) matrix_set(sk, 1, y, s_min);
#if 0
         if (matrix_get(sk, 0, y) > s_max) matrix_set(sk, 0, y, s_max);
         if (matrix_get(sk, 1, y) > s_max) matrix_set(sk, 1, y, s_max);
#endif
    }
    //compute output
    double n = 0;
    //sum all components
    for (size_t y = 0; y != N; ++y)
    {
        n += fabs(matrix_get(sk, 0, y) * matrix_get(sk, 1, y));
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
    long loop_count = 1;
    //..
    while (loop_count <= s_loop_limit)
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
    
    while (fabs(tmp[2] - tmp[1]) > 0.00000000001)
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
        matrix_multiply(t_x[0], x[0]), matrix_multiply(t_x[0], x[1]),
        matrix_multiply(t_x[1], x[0]), matrix_multiply(t_x[1], x[1])
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
    Matrix* U  = matrix_rotate(estimate_t.theta);
    Matrix* LU = matrix_multiply(L, U);
    //compute sum1/2
    double sum1 = matrix_get(LU, 0, 0) + matrix_get(LU, 1, 0);
    double sum2 = matrix_get(LU, 0, 1) + matrix_get(LU, 1, 1);
    //compute estimate
    double A_estimate_init[] =
    {
        matrix_get(LU, 0, 0)/sum1, matrix_get(LU, 1, 0)/sum1,
        matrix_get(LU, 0, 1)/sum2, matrix_get(LU, 1, 1)/sum2
    };
    Matrix* A_estimate = matrix_init(A_estimate_init,2, 2);
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
    matrix_free(LU);
    //return matrix
    return A_estimate;
}

ImageMarge split_images(ImageMarge images,  size_t nm[2], size_t nm2)
{
    //front/back to vector
    Matrix* front_col = matrix_to_vector_col(images.front);
    Matrix* back_col = matrix_to_vector_col(images.back);
    //inv front and back
    Matrix* front_inv_col =r_matrix_inverse(front_col);
    Matrix* back_inv_col  =r_matrix_inverse(back_col);
    //xorig link to front_inv and back_inv
    Matrix* x_orig[] =
    {
        front_inv_col,
        back_inv_col
    };
    //mode
    double mode_front = matrix_mode_col(front_inv_col, 0);
    double mode_back  = matrix_mode_col(back_inv_col, 0);
    //mode is the zero of all values
    Matrix* x[] =
    {
        matrix_alloc(1,front_inv_col->h),
        matrix_alloc(1,back_inv_col->h)
    };
    //..
    for (size_t i = 0; i != front_inv_col->h; ++i)
    {
        if (matrix_get(front_inv_col, 0, i) < mode_front)
            matrix_set(x[0], 0, i, 0);
        else
            matrix_set(x[0], 0, i, matrix_get(front_inv_col, 0, i) - mode_front);
        
        if (matrix_get(back_inv_col, 0, i) < mode_back)
            matrix_set(x[1], 0, i, 0);
        else
            matrix_set(x[1], 0, i, matrix_get(back_inv_col, 0, i) - mode_back);
    }
    //Processing
    Matrix* A_stimate     = estimate(x, nm2);
    Matrix* A_stimate_i   = matrix_inv2x2(A_stimate);
    Matrix* A_stimate_i_t = matrix_transpose(A_stimate_i);
    Matrix* x_orig2		  = put_x_into_matrix_col(x_orig, x_orig[0]->h);
    Matrix* sk			  = matrix_multiply(x_orig2, A_stimate_i_t);
    //get mode of images
    double	new_mode_front = matrix_mode_col(sk, 0);
    double	new_mode_back  = matrix_mode_col(sk, 1);
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
    matrix_free(front_inv_col);
    matrix_free(back_inv_col);
    matrix_free(front_col);
    matrix_free(back_col);
    //return
    return output;
}



