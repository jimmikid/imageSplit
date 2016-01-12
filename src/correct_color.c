//
//  correct_color.c
//  ImageProcessors
//
//  Created by Gianmarco Stinchi on 02/10/15.
//  Copyright © 2015 Gianmarco Stinchi. All rights reserved.
//
#include <matrix.h>
#include <tga_loader.h>
#include <correct_color.h>
#if 1
int    k = 100;  //step for the grid
double d1 = 5;  //minimun distance from mode
double d2 = 5; //minimun distance from original image
double c = 1;   //correction factor
#else
int    k = 90;  //step for the grid
double d1 = 5;  //minimun distance from mode
double d2 = 12; //minimun distance from original image
double c = 1;   //correction factor
#endif
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
    Matrix* rgb_cols[3] =
    {
        matrix_sub_to_vector_col(ss.matrix_r, pos, sub_nm),
        matrix_sub_to_vector_col(ss.matrix_g, pos, sub_nm),
        matrix_sub_to_vector_col(ss.matrix_b, pos, sub_nm)
    };
    
    //matrix_print_to_file(rgb_rows[0],"pisello.txt");
    
    double mode[3] =
    {
        matrix_mode_col(rgb_cols[0],0),
        matrix_mode_col(rgb_cols[1],0),
        matrix_mode_col(rgb_cols[2],0)
    };
    
    //printf("%g, %g, %g \n",mode[0],mode[1],mode[2]);
    
    for (size_t i = 0; i != 3; ++i)
    {
        matrix_free(rgb_cols[i]);
    }
    
    double mode_md = (mode[0] + mode[1] + mode[2]) / 3.0;
    //if(mode_md > 0.0)
    {
        for (size_t sub_i = 0; sub_i != sub_nm[0]; ++sub_i)
        for (size_t sub_j = 0; sub_j != sub_nm[1]; ++sub_j)
            {
                size_t i = sub_i + pos[0];
                size_t j = sub_j + pos[1];
                
                if( (fabs((matrix_get(ss.matrix_r, i, j) +
                           matrix_get(ss.matrix_g, i, j) +
                           matrix_get(ss.matrix_b, i, j)) / 3.0 - mode_md) > d1 &&
                     fabs(matrix_get(ss.matrix_r, i, j) - matrix_get(orig.matrix_r, i, j)) > d2 &&
                     fabs(matrix_get(ss.matrix_g, i, j) - matrix_get(orig.matrix_g, i, j)) > d2 &&
                     fabs(matrix_get(ss.matrix_b, i, j) - matrix_get(orig.matrix_b, i, j)) > d2) )
                    
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
