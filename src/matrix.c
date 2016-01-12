//
//  matrix.c
//  ImageProcessors
//
//  Created by Gianmarco Stinchi on 02/10/15.
//  Copyright © 2015 Gianmarco Stinchi. All rights reserved.
//

#include <matrix.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>
#include <math.h>

//Allocazione dinamica del tipo matrice
Matrix* matrix_alloc(size_t w,size_t h)
{
    //allocazione della struttura matrice
    Matrix* output = (Matrix*) malloc(sizeof(Matrix));
    output->w=w;
    output->h=h;
    
    //allocazione puntatori a colonne
    output->buffer = (double**)malloc(sizeof(double*)*w);
    
    //allocazione delle colonne
    for (size_t i = 0; i < w; i++)
    {
        output->buffer[i] = calloc(sizeof(double),h);
    }
    return output;
}
//Alloca ed inizializza da un array di double
Matrix* matrix_init(double* values,size_t w,size_t h)
{
    Matrix* output = matrix_alloc(w,h);
    for(size_t x=0;x!=w;++x)
    for(size_t y=0;y!=h;++y)
    {
        matrix_set(output, x, y, values[y*w+x]);
    }
    return output;
}
//create a rotation matrix
Matrix* matrix_rotate(double alpha)
{
	const double s = sin(alpha);
	const double c = cos(alpha);
	const double r_init[] =
	{
		s,-c,
		c, s
	};
	return matrix_init(r_init,2,2);
}
//Alloca una matrice di identità
Matrix* matrix_identity(size_t w,size_t h)
{
    //alloc
    Matrix* output = matrix_alloc(w,h);
    //get min size
    size_t mwh = w < h ? w : h;
    //fill all
    for(size_t xy=0;xy!=mwh;++xy)
        matrix_set(output, xy, xy, 1.0);
    //return matrix
    return output;
}
//Alloca una matrice sqrt(in)
Matrix* matrix_sqrt(const Matrix* in)
{
	Matrix* output = matrix_alloc(in->w, in->h);
	for (size_t x = 0; x != in->w; ++x)
	for (size_t y = 0; y != in->h; ++y)
	{
		matrix_set(output, x, y, sqrt(matrix_get(in,x,y)));
	}
	return output;
}
int compare_matrix_values(const void *x, const void *y)
{
    double xx = *(double*)x;
    double yy = *(double*)y;
    if (xx < yy) return -1;
    if (xx > yy) return  1;
    return 0;
}
static void merge_sort(Matrix* matrix, size_t col)
{
    qsort(matrix->buffer[col],
          matrix->h,
          sizeof(double),
          compare_matrix_values);
}
//Media moda
double matrix_mode_row(const Matrix* a,size_t row)
{
    //assert...
    assert(a->w);
    //alloc
    Matrix* tmp_m = matrix_alloc(1, a->w);
    //copy
    for(size_t x=0;x!=a->w;++x)
    {
        matrix_set(tmp_m,0,x,matrix_get(a,x,row));
    }
    //sort
	merge_sort(tmp_m,0);
    //comput mode
    size_t count      = 1;
    double output     = matrix_get(tmp_m,0,0);
    size_t tmp_count  = 0;
    double tmp_output = matrix_get(tmp_m,0,0);
    
    for(size_t y=1;y!=tmp_m->h;++y)
    {
        if( tmp_output == matrix_get(tmp_m,0,y) )
        {
            ++tmp_count;
        }
        else
        {
            if( tmp_count > count )
            {
                count  = tmp_count;
                output = tmp_output;
            }
            //release
            tmp_count  = 1;
            tmp_output = matrix_get(tmp_m,0,y);
        }
    }
    //last case
    if( tmp_count > count )
    {
        count  = tmp_count;
        output = tmp_output;
    }
    //free
    matrix_free(tmp_m);
    //return
    return output;
}

//Media moda
double matrix_mode_col(const Matrix* a,size_t col)
{
    //assert...
    assert(a->h);
    //alloc
    Matrix* tmp_m = matrix_alloc(1, a->h);
    //copy
    memcpy(tmp_m->buffer[0],a->buffer[col],sizeof(double)*a->h);
    //sort
    merge_sort(tmp_m,0);
    //comput mode
    size_t count      = 1;
    double output     = matrix_get(tmp_m,0,0);
    size_t tmp_count  = 0;
    double tmp_output = matrix_get(tmp_m,0,0);
    //compute 
    for(size_t y=1;y!=tmp_m->h;++y)
    {
        if( tmp_output == matrix_get(tmp_m,0,y) )
        {
            ++tmp_count;
        }
        else
        {
            if( tmp_count > count )
            {
                count  = tmp_count;
                output = tmp_output;
            }
            //release
            tmp_count  = 1;
            tmp_output = matrix_get(tmp_m,0,y);
        }
    }
    //last case
    if( tmp_count > count )
    {
        count  = tmp_count;
        output = tmp_output;
    }
    //free
    matrix_free(tmp_m);
    //return
    return output;
}

//fattorizazione di cholesky
Matrix* matrix_cholesky_factorization(const Matrix* a)
{
    if(a->w!=a->h) return NULL;
    Matrix* L = matrix_alloc(a->w,a->w);
    
    for (size_t i = 0; i != a->w; ++i)
    {
        for (size_t j = 0; j < (i+1); ++j)
        {
            double s = 0;
            
            for (size_t k = 0; k < j; ++k)
            {
                s += matrix_get(L, k, i) *  matrix_get(L, k, j);
                /* s += L[i * n + k] * L[j * n + k]; */
            }
            
            matrix_set(L, j, i,  (i == j) ?
                       sqrt(matrix_get(a, i, i)-s) :
                       1.0/ matrix_get(L, j, j) * (matrix_get(a, j, i)-s));
            /*
             L[i * n + j] = (i == j) ?
             sqrt(A[i * n + i] - s) :
             (1.0 / L[j * n + j] * (A[i * n + j] - s));
            */
        }
    }
    return L;
}

//Copia i valori di una matrice in un altra matrice
void matrix_put_from_sub_matrix(Matrix* dest,const Matrix* source,size_t pos[2])
{
    assert(source->w+pos[0] <= dest->w);
    assert(source->h+pos[1] <= dest->h);
    
    for(size_t x=0; x!=source->w ; ++x)
    for(size_t y=0; y!=source->h ; ++y)
    {
        matrix_set(dest, x + pos[0], y + pos[1], matrix_get(source,x,y));
    }
}
//alloca una matrice diagonale
Matrix* matrix_diagonal(size_t w,size_t h,double value)
{
    //alloc
    Matrix* output = matrix_alloc(w,h);
    //get min size
    size_t mwh = w < h ? w : h;
    //fill all
    for(size_t xy=0;xy!=mwh;++xy)
        matrix_set(output, xy, xy, value);
    //return matrix
    return output;
}

//copia la matrice
Matrix* matrix_copy(const Matrix* in)
{
    Matrix* output = matrix_alloc(in->w, in->h);
    for(size_t x=0;x!=in->w;++x)
    {
        memcpy(output->buffer[x], in->buffer[x], sizeof(double)*in->h);
    }
    return output;
}

//Moltiplicazione Matrice Matrice
Matrix* matrix_multiply(const Matrix* a, const Matrix* b)
{
    if(a->w != b->h) return NULL;
    
    Matrix* output = matrix_alloc(b->w, a->h);
    
    for (size_t y = 0; y != output->h; ++y)
    for (size_t x = 0; x != output->w; ++x)
    {
        double sum = 0;
        for (size_t k = 0; k != b->h; ++k)
            sum += matrix_get(a, k, y) * matrix_get(b, x, k);
        
        matrix_set(output, x, y, sum);
    }
    return output;
}
//Somma Matrice Matrice
Matrix* matrix_sum(const Matrix* a, const Matrix* b)
{
    if(a->w != b->w) return NULL;
    if(a->h != b->h) return NULL;
    
    Matrix* output = matrix_alloc(a->w, a->h);
    
    for(size_t y=0;y!=output->h;++y)
    for(size_t x=0;x!=output->w;++x)
    {
        matrix_set(output, x, y, matrix_get(a, x, y) + matrix_get(b, x, y) );
    }
    
    return output;
}
//matrix to vector
Matrix* matrix_to_vector(const Matrix* in)
{
    Matrix* output = matrix_alloc(in->w*in->h, 1);

	for (size_t x = 0; x != in->w; ++x)
	for (size_t y = 0; y != in->h; ++y)
    {
        matrix_set(output, y+x*in->h, 0, matrix_get(in,x,y));
    }
    
    return output;
}
//matrix to vector
Matrix* matrix_to_vector_col(const Matrix* in)
{
    Matrix* output = matrix_alloc(1,in->w*in->h);
    
    for (size_t x = 0; x != in->w; ++x)
    for (size_t y = 0; y != in->h; ++y)
    {
        matrix_set(output,  0 , y+x*in->h, matrix_get(in,x,y));
    }
    
    return output;
}
//sub vector
Matrix* matrix_sub_to_vector(const Matrix* in, size_t pos[2], size_t size[2])
{
    Matrix* output = matrix_alloc(size[0] * size[1], 1);
    
    for (size_t x = 0; x != size[0]; ++x)
        for (size_t y = 0; y != size[1]; ++y)
        {
            matrix_set(output, y + x*size[1], 0, matrix_get(in, x + pos[0], y + pos[1]));
        }
    
    return output;
}
//sub vector
Matrix* matrix_sub_to_vector_col(const Matrix* in, size_t pos[2], size_t size[2])
{
    Matrix* output = matrix_alloc(1, size[0] * size[1]);
    
    for (size_t x = 0; x != size[0]; ++x)
    for (size_t y = 0; y != size[1]; ++y)
    {
        matrix_set(output, 0, y + x*size[1], matrix_get(in, x + pos[0], y + pos[1]));
    }
    
    return output;
}
//row vector to matrix
Matrix* matrix_from_vector(const Matrix* in, size_t col,size_t w,size_t h)
{
    //filter
    if(w*h > in->h) return NULL;
    //alloc matrix
    Matrix* output = matrix_alloc(w,h);

	for (size_t x = 0; x != w; ++x)
	for (size_t y = 0; y != h; ++y)
    {
        matrix_set(output, x, y, matrix_get(in,col,y+x*h));
    }
    
    return output;
}
//Moltiplicazione Matrice pre scalare
Matrix* matrix_multiply_to_scalar(const Matrix* a, double scalar)
{
    Matrix* output = matrix_alloc(a->w, a->h);
    
    for (size_t y = 0; y != output->h; ++y)
    for (size_t x = 0; x != output->w; ++x)
    {
        matrix_set(output, x, y, matrix_get(a, x, y) * scalar);
    }
    
    return output;
}
//Moltiplicazione Matrice per scalare
void matrix_multiply_to_scalar_inplace(Matrix* a, double scalar)
{
    for (size_t y = 0; y != a->h; ++y)
    for (size_t x = 0; x != a->w; ++x)
    {
        matrix_set(a, x, y, matrix_get(a, x, y) * scalar);
    }
}
//Somma uno scalare
Matrix* matrix_sum_scalar(const Matrix* a, double scalar)
{
    Matrix* output = matrix_alloc(a->w, a->h);
    
    for (size_t y = 0; y != output->h; ++y)
    for (size_t x = 0; x != output->w; ++x)
    {
        matrix_set(output, x, y, matrix_get(a, x, y) + scalar);
    }
    
    return output;
}
//Somma uno scalare
void matrix_sum_scalar_inplace(Matrix* a, double scalar)
{
    for (size_t y = 0; y != a->h; ++y)
    for (size_t x = 0; x != a->w; ++x)
    {
        matrix_set(a, x, y, matrix_get(a, x, y) + scalar);
    }
}
//Imposta tutti i valori della matrice
void matrix_set(Matrix* in,size_t x,size_t y,double val)
{
    in->buffer[x][y] = val;
}

//restituisce tutti i valori della matrice
double matrix_get(const Matrix* in,size_t x,size_t y)
{
    return in->buffer[x][y];
}
//restituisce tutti i valori della matrice (sicuro)
double matrix_get_safe(const Matrix* in,size_t x,size_t y)
{
    x = x < in->w ? x : in->w-1;
    x = x < 0     ? 0 : x;
    
    y = y < in->h ? y : in->h-1;
    y = y < 0     ? 0 : y;
    
    return in->buffer[x][y];
}
//get normalize value
double matrix_get_safe_norm(const Matrix* in,size_t x,size_t y)
{
    return matrix_get_safe(in,x,y) / 255.0;
}

//Libera memoria allocata
void matrix_free(Matrix* in)
{
    for (size_t i = 0; i < in->w; i++)
    {
        free(in->buffer[i]);
    }
    free(in->buffer);
    free(in);
}

//legge da file gimp
Matrix* matrix_from_gimp_file (char *filename)
{
    char *res = NULL;
    FILE *fp  = NULL;
    char buf[256];                      //Buffer per la lettura dei dati da file
    size_t *dimensione;                    //allocazione di un array in memoria
    dimensione = malloc(2*sizeof(int)); //allocazione dinamica dell'array di due elementi(quelli che verranno letti nel file)
    fp = fopen(filename, "rb");          //apertura del file in lettura
    if(!fp) return NULL;
    
    while(res != NULL)
    {
        res = fgets(buf,256,fp);
        if(buf[0] != '#' && buf[0] != 'P')
        {
            sscanf(buf, "%d %d", &dimensione[0], &dimensione[1]); //lettura della dimensione della matrice sul file
            res = fgets(buf, 256,fp);                             //legge la riga da saltare dopo le dimensioni della matrice
            break;
        }
    }
    //Inizializzaione e allocazione dinamica della matrice
    Matrix* output = matrix_alloc(dimensione[0],dimensione[1]);
    
    for(size_t x=0;x!=output->w;++x)
        for(size_t y=0;y!=output->h;++y)
        {
            //read
            res = fgets(buf, 256, fp);
            int value = 0;
            sscanf(buf,"%d",&value);
            //set value
            matrix_set(output,x,y,(double)value);
        };
    
    fclose(fp);
    //ouput
    return output;
}

//salva matrice nel formato desiderato
void matrix_save_only_buffer(const char* path, const Matrix* in)
{
    FILE* file = fopen(path,"wb");
    
    if(file)
    {
        //buffer number
        char buffer[256] = { 0 };
        //scale of grey
        const int sgrey = 255;
        //constant header
        const char cheader[]=
        "P2\n"
        "%d %d\n"
        "%d\n";
        //final buffer of header
        char hbuffer[256] = { 0 };
        //final header
        snprintf(hbuffer,254,cheader,in->w,in->h,sgrey);
        //write header
        fwrite(hbuffer,sizeof(char),strlen(hbuffer),file);
        //write image
        for(size_t x=0;x!=in->w;++x)
            for(size_t y=0;y!=in->h;++y)
            {
                //set value
                double value = matrix_get(in,x,y);
                //int to string
                //snprintf (buffer,254,"%g\n",value);
                snprintf (buffer,254,"%d\n",(int)value);
                //write number
                fwrite(buffer,sizeof(char),strlen(buffer),file);
                
            };
    }
    
}


//print matrix
void matrix_print(const Matrix* in)
{
	for (size_t y = 0; y != in->h; ++y)
	{
		for (size_t x = 0; x != in->w; ++x)
		{
			printf("%0.3f\t", matrix_get(in, x, y));
		}
		printf("\n");
	}
}
void matrix_print_to_file(const Matrix* in,const char* path)
{
	FILE* pfile = fopen(path, "w");
	
	if (!pfile) return;

	for (size_t y = 0; y != in->h; ++y)
	{
		for (size_t x = 0; x != in->w; ++x)
		{
			fprintf(pfile,"%0.3f\t", matrix_get(in, x, y));
		}
		fprintf(pfile,"\n");
	}

	fclose(pfile);
}

//calcolo eigenvalues (autovalori) 2x2
Eigenvalues2x2 matrix_eigenvalues2x2(Matrix* in)
{
    Eigenvalues2x2 output;
    
    if(in->w != 2 || in->h != 2)
    {
        output.success = false;
        output.lambda1 = 0.0;
        output.lambda2 = 0.0;
        return output;
    }

    double a = matrix_get(in, 0, 0);
    double b = matrix_get(in, 1, 0);
    double c = matrix_get(in, 0, 1);
    double d = matrix_get(in, 1, 1);
    double value = (a+d);
    double det   = (d*a - c*b);
    double delta = sqrt(value*value-4.0*det);
    
    output.lambda1 = (value+delta)/2.0;
    output.lambda2 = (value-delta)/2.0;
    output.success = true;
    return output;
}

//calcolo degli autovettori su matrice 2x2
Eigenvctors2x2 matrix_eigenvectors2x2(Matrix* in)
{
    Eigenvctors2x2 output;
    
    if(in->w != 2 || in->h != 2)
    {
        output.success = false;
        output.v1 = NULL;
        output.v2 = NULL;
        return output;
    }
    
    Eigenvalues2x2 values = matrix_eigenvalues2x2(in);
    
    if(!values.success)
    {
        output.success = false;
        output.v1 = NULL;
        output.v2 = NULL;
        return output;
    }
    
    output.success = true;
    output.v1      = matrix_alloc(1, 2);
    output.v2      = matrix_alloc(1, 2);
    
    double b = matrix_get(in, 1, 0) ;
    double c = matrix_get(in, 0, 1) ;
    
    if( c != 0.0 )
    {
        matrix_set(output.v1, 0, 0,-matrix_get(in, 1, 1) + values.lambda1);
        matrix_set(output.v1, 0, 1, c);
        
        matrix_set(output.v2, 0, 0,-matrix_get(in, 1, 1) + values.lambda2);
        matrix_set(output.v2, 0, 1, c);
    }
    else if( b != 0.0 )
    {
        matrix_set(output.v1, 0, 0, b);
        matrix_set(output.v1, 0, 1,-matrix_get(in, 0, 0) + values.lambda1);
    
        matrix_set(output.v2, 0, 0, b);
        matrix_set(output.v2, 0, 1,-matrix_get(in, 0, 0) + values.lambda2);
    }
    else
    {
        matrix_set(output.v1, 0, 0, 1);
        matrix_set(output.v1, 0, 1, 0);
        
        matrix_set(output.v2, 0, 1, 0);
        matrix_set(output.v2, 0, 0, 1);
    }
    
    return output;
}
//calcolo degli autovettori su matrice 2x2 e li normalizza
Eigenvctors2x2 matrix_eigenvectors2x2_normalized(Matrix* in)
{
    Eigenvctors2x2 output = matrix_eigenvectors2x2(in);
    matrix_column_normalized_inplace(output.v1,0);
    matrix_column_normalized_inplace(output.v2,0);
    return output;
}
//calcolo gli autovettori (matrice 2x2) e gli auto valori (matrice 2*2, diagonale)
Eigen2x2 matrix_eigen2x2(Matrix* in)
{
	Eigen2x2 output;
	output.success = false;
	output.V = NULL;
	output.D = NULL;

	Eigenvalues2x2 values = matrix_eigenvalues2x2(in);
	if (!values.success) return output;
	Eigenvctors2x2 vectors = matrix_eigenvectors2x2_normalized(in);
	if (!vectors.success) return output;
	//create init matrix
	double init_values_v[] =
	{
		matrix_get(vectors.v2,0,0),matrix_get(vectors.v1,0,0),
		matrix_get(vectors.v2,0,1),matrix_get(vectors.v1,0,1),
	};
	//create init matrix
	double init_values_d[] =
	{
		values.lambda2,           0.0,
				   0.0,values.lambda1
	};
	//dealloc vectors
	eigenvectors2x2_free(vectors);
	//alloc V and D
	output.V = matrix_init(init_values_v, 2, 2);
	output.D = matrix_init(init_values_d, 2, 2);

	return output;
}
//Libera memoria allocata
void eigenvectors2x2_free(Eigenvctors2x2 values)
{
	if (values.success)
	{
		matrix_free(values.v1);
		matrix_free(values.v2);
		values.success = false;
		values.v1 = NULL;
		values.v2 = NULL;
	}
}
//Libera memoria allocata
void eigen2x2_free(Eigen2x2 values)
{
	if (values.success)
	{
		matrix_free(values.V);
		matrix_free(values.D);
		values.success = false;
		values.V = NULL;
		values.D = NULL;
	}
}

//Normalizzazione vettore
Matrix* matrix_column_normalized(const Matrix* in, size_t icolumn)
{
    Matrix* output = matrix_alloc(1, in->h);
    
    double sum = 0.0;
    //calcolo lunghezza del vettore
    for(size_t y=0; y!=in->h; ++y)
        sum += matrix_get(in,icolumn,y) * matrix_get(in,icolumn,y);
    
    double mag = sqrt(sum);
    //normalizzazione degli elementi del vettore
    for(size_t y=0; y!=output->h; ++y)
        matrix_set(output, icolumn, y, matrix_get(in,icolumn,y) / mag );
    
    return output;
}
//normalizza il vettore dato
void matrix_column_normalized_inplace(Matrix* inout, size_t icolumn)
{
    
    double sum = 0.0;
    //calcolo lunghezza del vettore
    for(size_t y=0; y!=inout->h; ++y)
        sum += matrix_get(inout,icolumn,y) * matrix_get(inout,icolumn,y);
    
    double mag = sqrt(sum);
    //normalizzazione degli elementi del vettore
    for(size_t y=0; y!=inout->h; ++y)
        matrix_set(inout, icolumn, y, matrix_get(inout,icolumn,y) / mag );
}
//trasporre matrice
Matrix* matrix_transpose(const Matrix* in)
{
    Matrix* output = matrix_alloc(in->h, in->w);
    
    for(size_t x = 0; x != output->w; ++x)
        for(size_t y = 0; y != output->h; ++y)
        {
            matrix_set(output, x, y, matrix_get(in, y, x));
        }
    return output;
}
//Inversa matrice 2x2 
Matrix* matrix_inv2x2(const Matrix* in)
{
	if (in->w != 2 || in->h != 2) return NULL;

	Matrix* output = matrix_alloc(2, 2);
	
	double det = (matrix_get(in, 0, 0)*matrix_get(in, 1, 1)) - 
		         (matrix_get(in, 1, 0)*matrix_get(in, 0, 1));

    double one_on_det = 1.0/det;
	matrix_set(output, 0, 0, (matrix_get(in, 1, 1)) * one_on_det);
    matrix_set(output, 1, 1, (matrix_get(in, 0, 0)) * one_on_det);
    matrix_set(output, 0, 1,-(matrix_get(in, 0, 1)) * one_on_det);
    matrix_set(output, 1, 0,-(matrix_get(in, 1, 0)) * one_on_det);
    //matrix_set(output, 1, 0, -(matrix_get(in, 1, 0)) / det);
    //matrix_set(output, 0, 1, -(matrix_get(in, 0, 1)) / det);

	return output;
}