//
//  main.c
//  clean-merge-app
//
//  Created by Gianmarco Stinchi on 28/12/15.
//  Copyright © 2015 Gianmarco Stinchi. All rights reserved.
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

bool compare_str(const char* arg, const char* command)
{
	return strncmp(arg, command, strlen(command))==0;
}

int main(int argc, const char* argv[])
{
	bool b_valid_arguments = true;
	bool b_show_help = false;
	bool b_merge = false;
	bool b_clean = false;
	const char* path_images_merge[] = { NULL, NULL };
	const char* path_images_merge_output[] = { NULL, NULL };
	const char* path_images_clean[] = { NULL, NULL };
	const char* path_images_clean_output[] = { NULL, NULL };

	for (int i = 1; i != argc; ++i)
	{
		if (compare_str(argv[i], "-h") || compare_str(argv[i], "--help"))
		{
			b_show_help = true;
		}
		else if (compare_str(argv[i], "-m") || compare_str(argv[i], "--merge"))
		{
			if (argc <= i + 4)
			{
				b_valid_arguments = false;
				break;
			}
			b_merge = true;
			path_images_merge[0] = argv[i + 1];
			path_images_merge[1] = argv[i + 2];
			path_images_merge_output[0] = argv[i + 3];
			path_images_merge_output[1] = argv[i + 4];
			i += 4;
		}
		else if (compare_str(argv[i], "-c") || compare_str(argv[i], "--clean"))
		{
			if (argc <= i + 4)
			{
				b_valid_arguments = false;
				break;
			}
			b_clean = true;
			path_images_clean[0] = argv[i + 1];
			path_images_clean[1] = argv[i + 2];
			path_images_clean_output[0] = argv[i + 3];
			path_images_clean_output[1] = argv[i + 4];
			i += 4;
		}
		else
		{
			b_valid_arguments = false;
			break;
		}
	}
	//print errors
	if (!b_valid_arguments)
	{
		printf("Not valid arguments\n");
		return -1;
	}

	//execute commands
	if (b_show_help)
	{
		printf("%s [options]\n", argv[0]);
		printf("Options:\n");
		printf("\t--help/-h help\n");
		printf("\t--clean/-c <path image 1> <path image 2> <path output image 1> <path ouput image 2> clean a merged images\n");
		printf("\t--merge/-m <path image 1> <path image 2> <path output image 1> <path ouput image 2> merge images\n");
	}

	if (b_merge)
	{
		//load images
		RGB_Matrix image[2];
		image[0] = rgb_matrix_from_tga_file(path_images_merge[0]);
		image[1] = rgb_matrix_from_tga_file(path_images_merge[1]);
		//marge images
		Matrix* colors[2] = { NULL, NULL };
		//call split images (red channel)
		colors[0] = image[0].matrix_r;
		colors[1] = image[1].matrix_r;
		ImageMarge  images_r = marge_images(colors);
		//call split images (green channel)
		colors[0] = image[0].matrix_g;
		colors[1] = image[1].matrix_g;
		ImageMarge  images_g = marge_images(colors);
		//call split images (blue channel)
		colors[0] = image[0].matrix_b;
		colors[1] = image[1].matrix_b;
		ImageMarge  images_b = marge_images(colors);
		////////////////////////////////////////////////////////////////////////////////////
		RGB_Matrix front_image = rgb_matrix_init(images_r.front, images_g.front, images_b.front);
		RGB_Matrix back_image = rgb_matrix_init(images_r.back, images_g.back, images_b.back);
		////////////////////////////////////////////////////////////////////////////////////
		rgb_matrix_to_tga_file(path_images_merge_output[0], front_image);
		rgb_matrix_to_tga_file(path_images_merge_output[1], back_image);
		//dealloc
		rgb_matrix_free(image[0]);
		rgb_matrix_free(image[1]);
		rgb_matrix_free(front_image);
		rgb_matrix_free(back_image);
	}

	if (b_clean)
	{
		//load images
		RGB_Matrix image[2];
		image[0] = rgb_matrix_from_tga_file(path_images_clean[0]);
		image[1] = rgb_matrix_from_tga_file(path_images_clean[1]);
		//assets
		assert(image[0].matrix_r->w == image[1].matrix_r->w &&
			image[0].matrix_r->h == image[1].matrix_r->h);
		//sizes
		size_t size[] =
		{
			image[0].matrix_r->w,
			image[0].matrix_r->h
		};
		////////////////////////////////////////////////////////////////////////////////////
		Matrix*     colors[2] = { NULL, NULL };
		//elements count
		size_t elements = size[0] * size[1];
		//call split images (red channel)
		colors[0] = image[0].matrix_r;
		colors[1] = image[1].matrix_r;
		ImageMarge  images_r = marge_images_init(colors);
		images_r = split_images(images_r, size, elements);
		//call split images (green channel)
		colors[0] = image[0].matrix_g;
		colors[1] = image[1].matrix_g;
		ImageMarge  images_g = marge_images_init(colors);
		images_g = split_images(images_g, size, elements);
		//call split images (blue channel)
		colors[0] = image[0].matrix_b;
		colors[1] = image[1].matrix_b;
		ImageMarge  images_b = marge_images_init(colors);
		images_b = split_images(images_b, size, elements);
		////////////////////////////////////////////////////////////////////////////////////
		RGB_Matrix front_image = rgb_matrix_init(images_r.front, images_g.front, images_b.front);
		RGB_Matrix back_image = rgb_matrix_init(images_r.back, images_g.back, images_b.back);
		////////////////////////////////////////////////////////////////////////////////////
		correct_color(front_image, image[0], size, elements);
		correct_color(back_image, image[1], size, elements);
		////////////////////////////////////////////////////////////////////////////////////
		rgb_matrix_to_tga_file(path_images_clean_output[0], front_image);
		rgb_matrix_to_tga_file(path_images_clean_output[1], back_image);
		////////////////////////////////////////////////////////////////////////////////////
		//dealloc
		rgb_matrix_free(image[0]);
		rgb_matrix_free(image[1]);
		rgb_matrix_free(front_image);
		rgb_matrix_free(back_image);
	}

	return 0;
}