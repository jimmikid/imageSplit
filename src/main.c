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
#include <split_image.h>
#include <correct_color.h>


int main(int argc,const char* argv[])
{
#if 1
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
    Matrix*     colors[2]={ NULL, NULL };
    //elements count
    size_t elements = size[0] * size[1];
	//marge images?
	bool b_not_merge = true;
    //call split images (red channel)
    colors[0] = image[0].matrix_r;
    colors[1] = image[1].matrix_r;
	ImageMarge  images_r = b_not_merge ? marge_images_init(colors) : marge_images(colors);
				images_r = split_images(images_r, size, elements);
    //call split images (green channel)
    colors[0] = image[0].matrix_g;
    colors[1] = image[1].matrix_g;
	ImageMarge  images_g = b_not_merge ? marge_images_init(colors) : marge_images(colors);
				images_g = split_images(images_g, size, elements);
    //call split images (blue channel)
    colors[0] = image[0].matrix_b;
    colors[1] = image[1].matrix_b;
	ImageMarge  images_b = b_not_merge ? marge_images_init(colors) : marge_images(colors);
				images_b = split_images(images_b, size, elements);

	////////////////////////////////////////////////////////////////////////////////////
    RGB_Matrix front_image = rgb_matrix_init(images_r.front, images_g.front, images_b.front);
    RGB_Matrix back_image = rgb_matrix_init(images_r.back, images_g.back, images_b.back);
    ////////////////////////////////////////////////////////////////////////////////////
    rgb_matrix_to_tga_file("out_front.tga", front_image);
    rgb_matrix_to_tga_file("out_back.tga", back_image);
	////////////////////////////////////////////////////////////////////////////////////
    correct_color(front_image, image[0], size, elements);
    correct_color(back_image, image[1], size, elements);
    ////////////////////////////////////////////////////////////////////////////////////
	rgb_matrix_to_tga_file("cr_out_front.tga", front_image);
	rgb_matrix_to_tga_file("cr_out_back.tga", back_image);
#else
	//paths
	const char* path_image[] =
	{
		"assets/mistura_recto.tga",
		"assets/mistura_verso.tga",

		"assets/verso_stimato.tga",
		"assets/recto_stimato.tga"
	};
    const int k = 100; //step for the grid
    const double d1 = 5;  //minimun distance from mode
    const double d2 = 5;  //minimun distance from original image
    const double c = 1;   //correction factor
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
	RGB_Matrix back_image = rgb_matrix_from_tga_file(path_image[3]);
    
    correct_color(front_image, image[0], size, elements);
    correct_color(back_image, image[1], size, elements);

	rgb_matrix_to_tga_file("_out_front.tga", front_image);
	rgb_matrix_to_tga_file("_out_back.tga", back_image);
#endif
    return 0;
}