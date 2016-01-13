//
//  old.split_image.c
//  ImageProcessors
//
//  Created by Gianmarco Stinchi on 02/10/15.
//  Copyright © 2015 Gianmarco Stinchi. All rights reserved.
//
#include <matrix.h>
#include <tga_loader.h>

/*
example
const int k = 100; //step for the grid
const double d1 = 5;  //minimun distance from mode
const double d2 = 5;  //minimun distance from original image
const double c = 1;   //correction factor
front_image = correct_color(front_image, image[0], k, d1, d2, c);
back_image = correct_color(back_image, image[1], k, d1, d2, c);
*/
void correct_rgb(RGB_Matrix ss,
                 RGB_Matrix orig,
                 size_t pos[2],
                 double d1,
                 double d2,
                 double c)
{
    
    size_t nm[]=
    {
        ss.matrix_r->w,
        ss.matrix_r->h
    };
    
    size_t nm2 = nm[0] * nm[1];
    
    Matrix* r_row = matrix_to_vector(ss.matrix_r);
    Matrix* g_row = matrix_to_vector(ss.matrix_g);
    Matrix* b_row = matrix_to_vector(ss.matrix_b);
    
    double mode_r = matrix_mode_row(r_row, 0);
    double mode_g = matrix_mode_row(g_row, 0);
    double mode_b = matrix_mode_row(b_row, 0);
    
    matrix_free(r_row);
    matrix_free(g_row);
    matrix_free(b_row);
    
    double mode_media = (mode_r + mode_g + mode_b)/3.0;
    size_t k = pos[0];
    size_t l = pos[1];
    
    for(size_t i = 0; i < nm[0]; ++i)
        for(size_t j = 0; j < nm[1]; ++j)
        {
            /*
             if(    abs(  (ss(i,j,1)+ss(i,j,2)+ss(i,j,3) )/3 - moda_s_media  ) > d1 &&
             abs(  ss(i,j,1)-orig(k+i-1,l+j-1,1)  )>d2 &&
             abs(  ss(i,j,2)-orig(k+i-1,l+j-1,2)  )>d2 &&
             abs(  ss(i,j,3)-orig(k+i-1,l+j-1,3)  )>d2 )
             */
            if(fabs((matrix_get(ss.matrix_r, i, j)+matrix_get(ss.matrix_g, i, j)+matrix_get(ss.matrix_b, i, j))/3.0 - mode_media) > d1 &&
               fabs(matrix_get(ss.matrix_r, i, j) - matrix_get(orig.matrix_r, k+i, l+j)) > d2 &&
               fabs(matrix_get(ss.matrix_g, i, j) - matrix_get(orig.matrix_g, k+i, l+j)) > d2 &&
               fabs(matrix_get(ss.matrix_b, i, j) - matrix_get(orig.matrix_b, k+i, l+j)) > d2
               )
            {
                /*
                 
                 ss(i,j,1) = moda_s_r + (c*(ss(i,j,1)-moda_s_r) / (ss(i,j,1)));
                 ss(i,j,2) = moda_s_g + (c*(ss(i,j,2)-moda_s_g) / (ss(i,j,2)));
                 ss(i,j,3) = moda_s_b + (c*(ss(i,j,3)-moda_s_b) / (ss(i,j,3)));
                 
                 */
                
                matrix_set(ss.matrix_r, i, j, mode_r + (c*(matrix_get(ss.matrix_r, i, j)-mode_r))/matrix_get(ss.matrix_r, i, j) );
                matrix_set(ss.matrix_g, i, j, mode_g + (c*(matrix_get(ss.matrix_g, i, j)-mode_g))/matrix_get(ss.matrix_g, i, j) );
                matrix_set(ss.matrix_b, i, j, mode_b + (c*(matrix_get(ss.matrix_b, i, j)-mode_b))/matrix_get(ss.matrix_b, i, j) );
                
            }
        }
    
}
RGB_Matrix correct_color(RGB_Matrix orig_s,
                         RGB_Matrix orig,
                         size_t kk,
                         double d1,
                         double d2,
                         double c)
{
    RGB_Matrix ss = rgb_matrix_copy(orig_s);
    
    size_t nm[]=
    {
        ss.matrix_r->w,
        ss.matrix_r->h
    };
    
    size_t nm2 = nm[0] * nm[1];
    
    if (kk==0)
    {
        size_t pos[] = { 0,0 };
        RGB_Matrix ss_sub     = rgb_matrix_sub(ss, pos, nm);
        correct_rgb(ss_sub, orig, pos, d1, d2, c);
        rgb_matrix_put_from_sub_matrix(ss,ss_sub,pos);
        matrix_rgb(ss_sub);
    }
    
    size_t i = 0;
    size_t j = 0;
    
    for ( i = 0; i < nm[0] - kk + 1; i += kk)
    {
        for ( j = 0; j < nm[1] - kk + 1; j += kk)
        {
            size_t pos[2]  = { i, j };
            size_t sub_nm[2] = { kk - 1, kk - 1 };
            RGB_Matrix ss_ij = rgb_matrix_sub(ss, pos, sub_nm);
            correct_rgb(ss_ij, orig, pos, d1, d2, c);
            rgb_matrix_put_from_sub_matrix(ss,ss_ij,pos);
            matrix_rgb(ss_ij);
            
        }
        size_t pos[2] = { i, j };
        size_t sub_nm[2] = { kk - 1, nm[1] - j };
        RGB_Matrix ss_ij = rgb_matrix_sub(ss, pos, sub_nm);
        correct_rgb(ss_ij, orig, pos, d1, d2, c);
        rgb_matrix_put_from_sub_matrix(ss,ss_ij,pos);
        matrix_rgb(ss_ij);
    }
    for ( j = 0; j < nm[1] - kk + 1; j += kk )
    {
        size_t pos[2] = { i, j };
        size_t sub_nm[2] = { nm[0] - i, kk-1 };
        RGB_Matrix ss_ij = rgb_matrix_sub(ss, pos, sub_nm);
        correct_rgb(ss_ij, orig, pos, d1, d2, c);
        rgb_matrix_put_from_sub_matrix(ss,ss_ij,pos);
        matrix_rgb(ss_ij);
    }
    size_t pos[2] = { i, j };
    size_t sub_nm[2] = { nm[0]-i, nm[1]-j };
    RGB_Matrix ss_ij = rgb_matrix_sub(ss, pos, sub_nm);
    correct_rgb(ss_ij, orig, pos, d1, d2, c);
    rgb_matrix_put_from_sub_matrix(ss,ss_ij,pos);
    matrix_rgb(ss_ij);
    
    return ss;
    
}
