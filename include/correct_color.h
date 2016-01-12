//
//  correct_color.h
//  TesiProject
//
//  Created by Gianmarco Stinchi on 28/12/15.
//  Copyright © 2015 Gianmarco Stinchi. All rights reserved.
//

#ifndef correct_color_h
#define correct_color_h

#include <stdlib.h>

struct  RGB_Matrix;
typedef struct RGB_Matrix RGB_Matrix;

extern int    k;   //step for the grid
extern double d1;  //minimun distance from mode
extern double d2 ; //minimun distance from original image
extern double c;   //correction factor

void correct_color(RGB_Matrix ss, RGB_Matrix orig, size_t nm[2], size_t nm2);

#endif /* correct_color_h */
