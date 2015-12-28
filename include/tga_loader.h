//
//  tga_loader.h
//  TesiProject
//
//  Created by Gianmarco Stinchi on 28/12/15.
//  Copyright © 2015 Gianmarco Stinchi. All rights reserved.
//

#ifndef tga_loader_h
#define tga_loader_h

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

RGB_Matrix RGB_Matrix_init(Matrix* matrix_r,Matrix* matrix_g,Matrix* matrix_b);
RGB_Matrix matrix_from_tga_file(const char* filepath);
bool matrix_to_tga_file(const char* filepath,RGB_Matrix rgb_matrix);
void matrix_rgb(RGB_Matrix rgb_matrix);

#endif /* tga_loader_h */
