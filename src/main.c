//
//  main.c
//  clean-merge-app
//
//  Created by Gianmarco Stinchi on 28/12/15.
//  Copyright � 2015 Gianmarco Stinchi. All rights reserved.
//


#include <stdio.h>
#include <stdlib.h>
#include <threads.h>
#include <string.h>
#include <matrix.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <tga_loader.h>
#include <split_image.h>
#include <correct_color.h>

typedef struct context
{
    ImageMarge images;
    size_t size[2];
    size_t elements;
}
context;

void* images_split_task(void* void_ptr_ctx)
{
    context* ctx = (context*)void_ptr_ctx;
    ctx->images = split_images(ctx->images, ctx->size, ctx->elements);
    return ctx;
}


thread* images_split_thread(ImageMarge images,
                            size_t size[2],
                            size_t elements)
{
    context* ctx = (context*)malloc(sizeof(context));
    ctx->elements = elements;
    ctx->images   = images;
    ctx->size[0]  = size[0];
    ctx->size[1]  = size[1];
    
    return execute_task(images_split_task, ctx);
}

ImageMarge images_split_join(thread* th)
{
    void* void_ptr_ctx=joint(th);
    assert(void_ptr_ctx);
    
    context* ctx = (context*)void_ptr_ctx;
    ImageMarge images = ctx->images;
    free(ctx);
    
    return images;
}

bool compare_str(const char* arg, const char* command)
{
	return strncmp(arg, command, strlen(command))==0;
}

int main(int argc, const char* argv[])
{
	bool b_valid_arguments = true;
	bool b_show_help       = false;
	bool b_merge           = false;
	bool b_clean           = false;
    bool b_parallel        = false;
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
        else if (compare_str(argv[i], "-p") || compare_str(argv[i], "--parallel"))
        {
            b_parallel = true;
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
        printf("\t--parallel/-p help\n");
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
        Matrix*     colors_r[2] = { NULL, NULL };
        Matrix*     colors_g[2] = { NULL, NULL };
        Matrix*     colors_b[2] = { NULL, NULL };
		//elements count
		size_t elements = size[0] * size[1];
        //channels
        ImageMarge  images_r;
        ImageMarge  images_g;
        ImageMarge  images_b;
        //call split images (red channel)
        colors_r[0] = image[0].matrix_r;
        colors_r[1] = image[1].matrix_r;
        images_r = marge_images_init(colors_r);
        //call split images (green channel)
        colors_g[0] = image[0].matrix_g;
        colors_g[1] = image[1].matrix_g;
        images_g = marge_images_init(colors_g);
        //call split images (blue channel)
        colors_b[0] = image[0].matrix_b;
        colors_b[1] = image[1].matrix_b;
        images_b = marge_images_init(colors_b);
        //execute
        if(b_parallel)
        {
            thread* th_r = images_split_thread(images_r, size, elements);
            thread* th_g = images_split_thread(images_g, size, elements);
            thread* th_b = images_split_thread(images_b, size, elements);
            images_r = images_split_join(th_r);
            images_g = images_split_join(th_g);
            images_b = images_split_join(th_b);
        }
        else
        {
            images_r = split_images(images_r, size, elements);
            images_g = split_images(images_g, size, elements);
            images_b = split_images(images_b, size, elements);
        }
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