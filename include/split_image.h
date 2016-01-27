//
//  split_image.h
//  ImageProcessors
//
//  Created by Gianmarco Stinchi on 02/10/15.
//  Copyright © 2015 Gianmarco Stinchi. All rights reserved.
//
#ifndef split_image_h
#define split_image_h
#include <stdlib.h>

struct Matrix;
typedef struct Matrix Matrix;

typedef struct ImageMarge
{
    Matrix* front;
    Matrix* back;
}
ImageMarge;


void set_loop_limit(long loop_limit);
long get_loop_limit();
ImageMarge split_images(ImageMarge images,  size_t nm[2], size_t nm2);
ImageMarge marge_images_init(Matrix* s[2]);
ImageMarge marge_images(Matrix* s[2]);

#endif