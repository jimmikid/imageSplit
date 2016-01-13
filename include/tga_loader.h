//
//  tga_loader.h
//  clean-merge-app
//
//  Created by Gianmarco Stinchi on 28/12/15.
//  Copyright © 2015 Gianmarco Stinchi. All rights reserved.
//

#ifndef tga_loader_h
#define tga_loader_h
#include <stdlib.h>
#include <stdbool.h>

struct Matrix;
typedef struct Matrix Matrix;

typedef struct RGB_Matrix
{
    bool    success;
    Matrix* matrix_r;
    Matrix* matrix_g;
    Matrix* matrix_b;
}
RGB_Matrix;
//rgb matrixs ops
RGB_Matrix rgb_matrix_init(Matrix* matrix_r,Matrix* matrix_g,Matrix* matrix_b);
RGB_Matrix rgb_matrix_from_tga_file(const char* filepath);
RGB_Matrix rgb_matrix_inverse(RGB_Matrix rgb_matrix);
RGB_Matrix rgb_matrix_normalize(RGB_Matrix rgb_matrix);
RGB_Matrix rgb_matrix_denormalize(RGB_Matrix rgb_matrix);
RGB_Matrix rgb_matrix_copy(RGB_Matrix rgb_matrix);
RGB_Matrix rgb_matrix_sub(RGB_Matrix rgb_matrix, size_t pos[2], size_t size[2]);
void rgb_matrix_put_from_sub_matrix(RGB_Matrix dest,const RGB_Matrix source,size_t pos[2]);
void rgb_matrix_clamp_inplace(RGB_Matrix rgb_matrix, double min, double max);
bool rgb_matrix_to_tga_file(const char* filepath,RGB_Matrix rgb_matrix);
void rgb_matrix_free(RGB_Matrix rgb_matrix);
//single channel ops
Matrix*    r_matrix_inverse(Matrix* r_matrix);
Matrix*    r_matrix_sub(Matrix* r_matrix, size_t pos[2], size_t size[2]);

#endif /* tga_loader_h */
