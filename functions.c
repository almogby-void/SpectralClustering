#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>

#include "spkmeans.h"

#define epsilon 0.00001
#define max_rotations 100





int main(int argc, char **argv) {
    double **M, **L, **V,**U;
    double *eigenvalues;
    int n,dim;
    M = file_to_matrix("test_input.txt",&n,&dim);
    V = Identity(n);
    L = L_norm(M,n,dim);
    eigenvalues = diag(Jacobi(L,V,n,0),n);
    print_matrix(V,n,n);
    U = eigenvectors(eigenvalues,V,n);
    normalize(U,n,3);
    print_matrix(U,n,3);
    argv[0] = argv[1];
    return argc; 
} 