#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>

double **rand_matrix();
void print_matrix(double ** matrix);
double **adj_matrix(double ** matrix);
double **diag_degree_matrix(double ** matrix);
static double dist(double *X, double *Y, int dim);
double **matrix_mult(double ** A,double ** B);
double **matrix_func(double ** matrix,double (*f)(double));
double **L_norm(double ** matrix);
double divsqrt(double A);


int main(int argc, char *argv[]) {
    double **random;

    random = rand_matrix();
    print_matrix(random);
    print_matrix(L_norm(random));
    argv[0] = argv[1];
    return argc;
} 

double **rand_matrix(){
    int i, j;
    double **matrix;
    matrix = calloc(3,sizeof(double*));
     
    for(i = 0; i < 3; i++) {
        matrix[i] = malloc(sizeof(double*) * 3);
    }
 
    for(i = 0; i < 3; i++){
        for(j = 0; j < 3; j++){
            matrix[i][j] = rand()%10;
        }
    }
    return matrix;
}

double **adj_matrix(double ** matrix){
    int i, j;
    double **adj;
    adj = calloc(3,sizeof(double*));
     
    for(i = 0; i < 3; i++) {
        adj[i] = calloc(3,sizeof(double*) * 3);
    }
 
    for(i = 0; i < 3; i++){
        for(j = 0; j < 3; j++){
            adj[i][j] =  exp(-dist(matrix[i],matrix[j],3)/2);
        }
    }
    return adj;
}

double **diag_degree_matrix(double ** matrix){
    int i, j;
    double **diag;
    diag = calloc(3,sizeof(double*) * 3);
     
    for(i = 0; i < 3; i++) {
        diag[i] = calloc(3,sizeof(double*));
    }
 
    for(i = 0; i < 3; i++){
        for(j = 0; j < 3; j++){
            diag[i][i] +=  exp(-dist(matrix[i],matrix[j],3)/2);
        }
    }
    return diag;
}

double **L_norm(double ** matrix){
    int i, j;
    double **L,**W,**D;
    L = calloc(3,sizeof(double*) * 3);
     
    for(i = 0; i < 3; i++) {
        L[i] = calloc(3,sizeof(double*));
    }
    W = adj_matrix(matrix);
    D = matrix_func(diag_degree_matrix(matrix),divsqrt);
    L = matrix_mult(matrix_mult(D,W),D);
    for(i = 0; i < 3; i++){
        for(j = 0; j < 3; j++){
            L[i][j] = (i==j) - L[i][j];
        }
    }
    return L;
}

double **matrix_func(double ** matrix,double (*f)(double)){
    int i, j;
    double **result;
    result = calloc(3,sizeof(double*) * 3);
     
    for(i = 0; i < 3; i++) {
        result[i] = calloc(3,sizeof(double*));
    }
 
    for(i = 0; i < 3; i++){
        for(j = 0; j < 3; j++){
            result[i][j] = f(matrix[i][j]);
        }
    }
    return result;
}

double divsqrt(double A){return A ? 1/sqrt(A) : 0;}

double **matrix_mult(double ** A,double ** B){
    int i, j,k; 
    double **result;
    result = calloc(3,sizeof(double*) * 3);
     
    for(i = 0; i < 3; i++) {
        result[i] = calloc(3,sizeof(double*));
    }
 
    for(i = 0; i < 3; i++){
        for(j = 0; j < 3; j++){
            for(k = 0; k < 3; k++){
                result[i][j] += A[i][k]*B[k][j];
            }
        }
    }
    return result;
}

static double dist(double *X, double *Y, int dim){
    double result = 0;
    int i = 0;
    for (; i < dim; i++)
        result += (X[i]-Y[i])*(X[i]-Y[i]);
    return sqrt(result);
}

void print_matrix(double ** matrix){
    int i, j;
    for(i = 0; i < 3; i++){
        for(j = 0; j < 3; j++)
            printf("%.3f,",matrix[i][j]);
        printf("%s","\n");
    }
}