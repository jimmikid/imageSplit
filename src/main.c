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
#include <tga_loader.h>

#if 1


int main()
{
    const char* file_img1 = "assets/D.tga";
    const char* file_img2 = "assets/test2.tga";
    //load image
    RGB_Matrix tga = matrix_from_tga_file(file_img1);
    
    //print
    for(short y=0;y!=tga.matrix_r->h;++y)
    {
        for(short x=0;x!=tga.matrix_r->w; ++x)
            printf(matrix_get(tga.matrix_r,x,y) < 0.5 ? "#" : " ");
        //new line
        printf("\n");
    }
    
    matrix_to_tga_file(file_img2,tga);
    
    return 0;
}

#else
double deleteLT128(double value, int x, int y,Matrix* in)
{
    return value < 128 ?  0 : value;
}
double deleteGT128(double value, int x, int y,Matrix* in)
{
    return value > 128 ?  0 : value;
}

int main()
{
    printf("Inserisci nome file: ");
    char filename[64];
    scanf("%s", filename);
    Matrix* mat = matrix_from_gimp_file(filename);
    if(mat)
    {
        matrix_save_only_buffer("gerace_prova.pgm",mat);
        matrix_free(mat);
    }
    else
        printf("errore caricamento file. ciao!\n");

    Matrix *A = matrix_alloc(2, 2);
    matrix_set(A, 0, 0, -5); matrix_set(A, 1, 0, -2);
    matrix_set(A, 0, 1, -2); matrix_set(A, 1, 1, 3);
    
    Eigenvalues2x2 eigen = matrix_eigenvalues2x2(A);
    Matrix *E = matrix_alloc(2, 2);
    matrix_set(E, 0, 0, eigen.lambda1); matrix_set(E, 1, 0, 0);
    matrix_set(E, 0, 1, 0); matrix_set(E, 1, 1, eigen.lambda2);
    
    Eigenvctors2x2 vectors = matrix_eigenvectors2x2(A);
    
    Matrix* V = matrix_alloc(2, 2);
    
    if(vectors.success)
    {
        Matrix* v1n = matrix_column_normalized(vectors.v1,0);
        Matrix* v2n = matrix_column_normalized(vectors.v2,0);
        matrix_set(V, 0, 0, matrix_get(v1n, 0, 0)); matrix_set(V, 1, 0, matrix_get(v2n, 0, 0));
        matrix_set(V, 0, 1, matrix_get(v1n, 0, 1)); matrix_set(V, 1, 1, matrix_get(v2n, 0, 1));
        //dealloc
        matrix_free(v1n);
        matrix_free(v2n);
    }
    else printf("autovettori non trovati");
    
    Matrix* Vt = matrix_transpose(V);
    Matrix* VE = matrix_multiply(V, E);
    Matrix* VEVt = matrix_multiply(VE, Vt);
    
    printf("A:\n");
    matrix_print(A);
    
    printf("VEVt:\n");
    matrix_print(VEVt);
    //dealloc all
    eigenvctors2x2_free(vectors);
    matrix_free(A);
    matrix_free(E);
    matrix_free(V);
    matrix_free(Vt);
    matrix_free(VEVt);
    
    return 0;
}
#endif