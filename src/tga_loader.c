//
//  tga_loader.c
//  TesiProject
//
//  Created by Gabriele Di Bari on 28/12/15.
//  Copyright © 2015 Gianmarco Stinchi. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <matrix.h>
#include <tga_loader.h>
#include <memory.h>
#include <assert.h>

#ifdef _MSC_VER
#define ASPACKED( __Declaration__ ) __pragma( pack(push,1) ) __Declaration__   __pragma( pack(pop) )
#else
#define ASPACKED( __Declaration__ ) __Declaration__ __attribute__((packed))
#endif

typedef unsigned char Byte;

ASPACKED(struct TgaHeader
{
    Byte  idlength;                // size of ID field that follows 18 byte header (0 usually)
    Byte  colourmaptype;           // type of colour map 0=none, 1=has palette
    Byte  imagetype;               // type of image 0=none,1=indexed,2=rgb,3=grey,+8=rle packed
    short int colourmaporigin;     // first colour map entry in palette
    short int colourmaplength;     // number of colours in palette
    Byte  colourmapbits;           // number of bits per palette entry 15,16,24,32
    short int xorigin;             // image x origin
    short int yorigin;             // image x origin
    short width;                   // image width in pixels
    short height;                  // image height in pixels
    Byte  bitsperpixel;            // image bits per pixel 8,16,24,32
    Byte  descriptor;              // image descriptor bits (vh flip bits)
});

RGB_Matrix rgb_matrix_from_tga_file(const char* filepath)
{
    RGB_Matrix output;
    output.success  = false;
    output.matrix_r = NULL;
    output.matrix_g = NULL;
    output.matrix_b = NULL;
    //image struct
    struct TgaHeader header;
    memset(&header,0,sizeof(struct TgaHeader));
    //open file
    FILE* pfile=fopen(filepath, "rb");
    //errors?
    if(!pfile) return output;
    //read header
    if(!fread(&header, sizeof(struct TgaHeader), 1, pfile))
    {
        //fail to read header
        fclose(pfile);
        return output;
    }
    //is rgb?
    if(header.imagetype!=2){ fclose(pfile); return output; }
    //24 or 32 bit image?
    if(header.bitsperpixel/8 < 3){ fclose(pfile); return output; }
    //alloc
    output.matrix_r = matrix_alloc(header.width, header.height);
    output.matrix_g = matrix_alloc(header.width, header.height);
    output.matrix_b = matrix_alloc(header.width, header.height);
    //channels
    size_t channels = header.bitsperpixel/8;
    //size of image
    size_t size = header.width * header.height * channels;
    //alloc
    Byte* data = (Byte*) calloc(size,1);
    //read data
    if(!fread(data, size, 1, pfile))
    {
        //fail to read data
        fclose(pfile);
        //safe delete
        free(data);
        //fail
        return output;
    }
    //close file
    fclose(pfile);
    //success to read
    output.success = true;
    //from data to matrixes
	size_t x = 0;
	size_t y = 0;
	//tga has flip image
    for(y=0;y!=header.height; ++y)
    {
        for(x=0;x!=header.width; ++x)
        {
            //get pixel
            Byte* pixel = &data[(x+y*header.width)*channels];
            //pixel[0] B
            //pixel[1] G
            //pixel[2] R
            //pixel[3] A
#if 0
			matrix_set(output.matrix_r, x, header.height - y - 1, ((double)pixel[2])/*/255.0*/);
			matrix_set(output.matrix_g, x, header.height - y - 1, ((double)pixel[1])/*/255.0*/);
			matrix_set(output.matrix_b, x, header.height - y - 1, ((double)pixel[0])/*/255.0*/); 
#else
			matrix_set(output.matrix_r, x, y , ((double)pixel[2])/*/255.0*/);
			matrix_set(output.matrix_g, x, y , ((double)pixel[1])/*/255.0*/);
			matrix_set(output.matrix_b, x, y , ((double)pixel[0])/*/255.0*/);
#endif
        };
    }
    //dealloc
    free(data);
    //close file
    fclose(pfile);
    
    return  output;
}

static double d_max(double a, double b)
{
	return a > b ? a : b;
}
static double d_min(double a, double b)
{
	return a < b ? a : b;
}
static double clamp(double value, double v_min, double v_max)
{
	return d_min(d_max(value, v_min), v_max);
}

void rgb_matrix_clamp_inplace(RGB_Matrix rgb_matrix, double min, double max)
{
	for (size_t y = 0; y != rgb_matrix.matrix_r->h; ++y)
	for (size_t x = 0; x != rgb_matrix.matrix_r->w; ++x)
	{
		matrix_set(rgb_matrix.matrix_r, x, y, clamp(matrix_get(rgb_matrix.matrix_r, x, y), min, max));
		matrix_set(rgb_matrix.matrix_g, x, y, clamp(matrix_get(rgb_matrix.matrix_g, x, y), min, max));
		matrix_set(rgb_matrix.matrix_b, x, y, clamp(matrix_get(rgb_matrix.matrix_b, x, y), min, max));
	}
}
bool rgb_matrix_to_tga_file(const char* filepath,RGB_Matrix rgb_matrix)
{
    //image struct
    struct TgaHeader header;
    memset(&header,0,sizeof(struct TgaHeader));
    //safe size
    if(rgb_matrix.matrix_r->w != rgb_matrix.matrix_g->w ||
       rgb_matrix.matrix_r->w != rgb_matrix.matrix_b->w ||
       rgb_matrix.matrix_r->h != rgb_matrix.matrix_b->h ||
       rgb_matrix.matrix_r->h != rgb_matrix.matrix_b->h) return false;
    //open file
    FILE* pfile=fopen(filepath, "wb");
    //errors?
    if(!pfile) return false;
    //channels
    size_t channels = 3;
    //size
    header.width        = rgb_matrix.matrix_r->w;
    header.height       = rgb_matrix.matrix_r->h;
    header.bitsperpixel = channels*8;//8 8 8
    header.imagetype    = 2;         //r g b
    //size of image
    size_t size = header.width * header.height * channels;
    //alloc
    Byte* data = (Byte*) calloc(size,1);
    //from matrixes to data
    size_t x = 0;
	size_t y = 0;
    for(y=0;y!=header.height; ++y)
    for(x=0;x!=header.width; ++x)
    {
#if 0
        data[(x+y*header.width)*channels+2] = (Byte)clamp(matrix_get(rgb_matrix.matrix_r, x, header.height - y - 1),0.0,255.0);
        data[(x+y*header.width)*channels+1] = (Byte)clamp(matrix_get(rgb_matrix.matrix_g, x, header.height - y - 1),0.0,255.0);
        data[(x+y*header.width)*channels+0] = (Byte)clamp(matrix_get(rgb_matrix.matrix_b, x, header.height - y - 1),0.0,255.0);
#else
		data[(x + y*header.width)*channels + 2] = (Byte)clamp(matrix_get(rgb_matrix.matrix_r, x, header.height - y - 1), 0.0, 255.0);
		data[(x + y*header.width)*channels + 1] = (Byte)clamp(matrix_get(rgb_matrix.matrix_g, x, header.height - y - 1), 0.0, 255.0);
		data[(x + y*header.width)*channels + 0] = (Byte)clamp(matrix_get(rgb_matrix.matrix_b, x, header.height - y - 1), 0.0, 255.0);
#endif
    }
    //write header
    fwrite(&header, sizeof(struct TgaHeader), 1, pfile);
    //write raw image
    fwrite(data, size, 1, pfile);
    //dealloc
    free(data);
    //close file
    fclose(pfile);
    return true;
}

RGB_Matrix rgb_matrix_init(Matrix* matrix_r,Matrix* matrix_g,Matrix* matrix_b)
{
    RGB_Matrix output;
    output.success = true;
    output.matrix_r = matrix_r;
    output.matrix_g = matrix_g;
    output.matrix_b = matrix_b;
    return output;
}

RGB_Matrix rgb_matrix_copy(RGB_Matrix rgb_matrix)
{
    RGB_Matrix output_rgb_matrix = rgb_matrix;
    output_rgb_matrix.matrix_r = matrix_copy(rgb_matrix.matrix_r);
    output_rgb_matrix.matrix_g = matrix_copy(rgb_matrix.matrix_g);
    output_rgb_matrix.matrix_b = matrix_copy(rgb_matrix.matrix_b);
    return output_rgb_matrix;
}

RGB_Matrix rgb_matrix_inverse(RGB_Matrix rgb_matrix)
{
    RGB_Matrix output = rgb_matrix_copy(rgb_matrix);
    
    for(size_t y=0;y!=output.matrix_r->h;++y)
    for(size_t x=0;x!=output.matrix_r->w;++x)
    {
        matrix_set(output.matrix_r, x, y, 255.0-matrix_get(output.matrix_r, x, y));
        matrix_set(output.matrix_g, x, y, 255.0-matrix_get(output.matrix_g, x, y));
        matrix_set(output.matrix_b, x, y, 255.0-matrix_get(output.matrix_b, x, y));
    }
    
    return output;
}

Matrix* r_matrix_inverse(Matrix* r_matrix)
{
    Matrix* output = matrix_copy(r_matrix);
    
    for(size_t y=0;y!=output->h;++y)
    for(size_t x=0;x!=output->w;++x)
    {
        matrix_set(output, x, y, 255.0-matrix_get(r_matrix, x, y));
    }
    
    return output;
}

Matrix* r_matrix_sub(Matrix* r_matrix, size_t pos[2], size_t size[2])
{
	assert(pos[0] + size[0] <= r_matrix->w);
	assert(pos[1] + size[1] <= r_matrix->h);
	Matrix* output = matrix_alloc(size[0], size[1]);
#if 1
	for (size_t x = 0; x != size[0]; ++x)
	{
		memcpy(output->buffer[x], &r_matrix->buffer[x+pos[0]][pos[1]], sizeof(double)*size[1]);
	}
#else
	for (size_t x = 0; x != size[0]; ++x)
	for (size_t y = 0; y != size[1]; ++y)
	{
		matrix_set(output, x, y, matrix_get(r_matrix, x+pos[0], y+pos[1]));
	}
#endif
	return output;
}

RGB_Matrix rgb_matrix_normalize(RGB_Matrix rgb_matrix)
{
    RGB_Matrix output = rgb_matrix_copy(rgb_matrix);
    
    for(size_t y=0;y!=output.matrix_r->h;++y)
        for(size_t x=0;x!=output.matrix_r->w;++x)
        {
            matrix_set(output.matrix_r, x, y, matrix_get(output.matrix_r, x, y) / 255.0);
            matrix_set(output.matrix_g, x, y, matrix_get(output.matrix_g, x, y) / 255.0);
            matrix_set(output.matrix_b, x, y, matrix_get(output.matrix_b, x, y) / 255.0);
        }
    
    return output;
}

RGB_Matrix rgb_matrix_denormalize(RGB_Matrix rgb_matrix)
{
    RGB_Matrix output = rgb_matrix_copy(rgb_matrix);
    
    for(size_t y=0;y!=output.matrix_r->h;++y)
    for(size_t x=0;x!=output.matrix_r->w;++x)
    {
        matrix_set(output.matrix_r, x, y, matrix_get(output.matrix_r, x, y) * 255.0);
        matrix_set(output.matrix_g, x, y, matrix_get(output.matrix_g, x, y) * 255.0);
        matrix_set(output.matrix_b, x, y, matrix_get(output.matrix_b, x, y) * 255.0);
    }
    
    return output;
}

RGB_Matrix rgb_matrix_sub(RGB_Matrix rgb_matrix,size_t pos[2],size_t size[2])
{
	RGB_Matrix output;
	output.matrix_r = r_matrix_sub(rgb_matrix.matrix_r,pos,size);
	output.matrix_g = r_matrix_sub(rgb_matrix.matrix_g,pos,size);
	output.matrix_b = r_matrix_sub(rgb_matrix.matrix_b,pos,size);
	return output;
}

void matrix_rgb(RGB_Matrix rgb_matrix)
{
    matrix_free(rgb_matrix.matrix_r);
    matrix_free(rgb_matrix.matrix_g);
    matrix_free(rgb_matrix.matrix_b);
}