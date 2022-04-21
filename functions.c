#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>

#define epsilon 0.00001
#define max_rotations 100

double **rand_matrix(int n);
double **symmetric(double ** matrix,int n);
double **Identity(int n);
void print_matrix(double ** matrix,int n);
double **adj_matrix(double ** matrix);
double **diag_degree_matrix(double ** matrix);
static double dist(double *X, double *Y, int dim);
double **matrix_mult(double ** A,double ** B);
double **matrix_func(double ** matrix,double (*f)(double));
double **L_norm(double ** matrix);
double divsqrt(double A);
double **Jacobi_Rotation_sMatrix(double ** A,int i, int j, int n);
double **Jacobi(double **A, int n, int iter);
double off(double ** matrix,int n);
void normalize(double **A, int n, int k);
double *diag(double ** matrix,int n);


int main(int argc, char *argv[]) {
    double **random;
    random = rand_matrix(3);
    random = symmetric(random,3);
    print_matrix(random,3);
    print_matrix(Jacobi(random,3,0),3);
    argv[0] = argv[1];
    return argc;
} 

double **rand_matrix(int n){
    int i, j;
    double **matrix;
    matrix = calloc(n,sizeof(double*));
     
    for(i = 0; i < n; i++) {
        matrix[i] = malloc(sizeof(double*) * n);
    }
 
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            matrix[i][j] = (rand()*rand())%10;
        }
    }
    return matrix;
}
double **symmetric(double ** matrix,int n){/*TODO*/ 
    int i, j;
    for(i = 0; i < n; i++){
        for(j = 0; j < i+1; j++){
            matrix[i][j] = matrix[j][i];
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
            adj[i][j] =  (i!=j)*exp(-dist(matrix[i],matrix[j],3)/2);
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

double *diag(double ** matrix,int n){
    int i;
    double *diag = calloc(n,sizeof(double));
    for (i = 0; i < n; i++)
    {
        diag[i] = matrix[i][i];
    }
    return diag;
}

double **Identity(int n){ /*identity matrix*/
    int i;
    double **m;
    /* init matrix*/
    m = calloc(3,sizeof(double*) * n); 
    for(i = 0; i < n; i++)
        m[i] = calloc(3,sizeof(double*));
    
    for(i = 0; i < n; i++)
        m[i][i] = 1;
    return m;
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

double **Jacobi_Rotation_Matrix(double ** A,int i, int j, int n){ /* Unused? */
    double **m;
    double c,theta,t;
    m = Identity(n);
    theta = (A[i][i]-A[j][j])/(2*A[i][j]);
    t = ((theta >= 0) - (theta < 0))/(abs(theta) + sqrt(theta*theta + 1));
    c = 1/sqrt(t*t+1);
    m[i][i] = c; m[j][j] = c;
    m[j][i] = ((i > j) - (i < j))*t*c; m[i][j] = -m[j][i];
    return m;
}

double **Jacobi(double **A, int n, int iter){
    int i,j;
    int row,col;
    double **temp; /* A', temp[0]=A'[i], temp[1]=A'[j]*/
    double c,theta,t,s,square_diff,old_off;
    temp = calloc(2,sizeof(double*));
    temp[0] = malloc(sizeof(double*) * n);
    temp[1] = malloc(sizeof(double*) * n);
    i = 1;
    j = 0;
    old_off = off(A,n);
    for (row = 0; row < n; row++){
        for (col = 0; col < n; col++){
            if ((row != col)  & (abs(A[row][col]) > abs(A[i][j]))){ 
                i = row;
                j = col;
            }
        }
    }
    theta = (A[j][j]-A[i][i])/(2*A[i][j]);
    t = ((theta >= 0) - (theta < 0))/(abs(theta) + sqrt(theta*theta + 1));
    c = 1/sqrt(t*t+1);
    s = t*c;
    for (col = 0; col < n; col++){ /* switcheroo'd from the instructions */
        temp[0][col] = c*A[i][col] - s*A[j][col];
        temp[1][col] = c*A[j][col] + s*A[i][col];
    }
    temp[0][i] = c*c*A[i][i] + s*s*A[j][j] - 2*s*c*A[i][j];
    temp[1][j] = s*s*A[i][i] + c*c*A[j][j] + 2*s*c*A[i][j];
    temp[0][j] = 0;/*(c*c-s*s)*A[i][j]+s*c*(A[i][i]-A[j][j]); TODO: Check if true */
    temp[1][i] = temp[0][j];

    for ( col = 0; col < n; col++){
        if ((col != i) & (col != j)){
            square_diff += 2*(A[i][col]*A[i][col]  - temp[0][col]*temp[0][col]); /* I calculate convergence here instead of the offsets individually as there are a lot of duplicates. */
            square_diff += 2*(A[j][col]*A[j][col]  - temp[1][col]*temp[1][col]);
            A[i][col] = temp[0][col];
            A[col][i] = temp[0][col];
            A[j][col] = temp[1][col];
            A[col][j] = temp[1][col];
        }  
    }
    square_diff += A[i][i]*A[i][i] + A[j][j]*A[j][j] + 2*A[i][j]*A[i][j] - temp[0][i]*temp[0][i] - temp[1][j]*temp[1][j] - 2*temp[0][j]*temp[0][j];
    A[i][i] = temp[0][i];
    A[j][j] = temp[1][j];
    A[i][j] = temp[0][j];
    A[j][i] = A[i][j];
    
    print_matrix(A,n);
    square_diff+=1;
    square_diff = old_off - off(A,n);
    printf("theta:%f,t:%f,c:%f,s:%f,old_off:%f,off:%f,diff:%f,i=%d,j=%d\n",theta,t,c,s,old_off,off(A,n),square_diff,i,j);
    if ((square_diff <= epsilon) | (++iter > max_rotations))
        return A;
    return Jacobi(A,n,iter);
}

int cmpfunc (const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
}


int Heuristic (double *list, int n){
    int i,argmax;
    double delta,max_delta = 0;
    n = (int)floor(n/2);
    qsort(list,n,sizeof(double),cmpfunc);
    for (i = 0; i < n; i++){
        delta = (list[i] - list[i+1]);
        if (delta > max_delta){
            max_delta = delta;
            argmax = i;
        }
    }
    return argmax;
}

void normalize(double **A, int n, int k){
    int i,j;
    double sum = 0;
    for (i = 0; i < n; i++){
        sum = 0;
        for (j = 0; j < k; j++){
            sum+=A[i][j]*A[i][j];
        }
        for (j = 0; j < k; j++){
            A[i][j] = A[i][j]/sqrt(sum);
        }
    }
}

double off(double **matrix, int n) {
    double sum = 0;
    int i, j;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            sum += (i != j) * matrix[i][j] * matrix[i][j];
    return sum;
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

void print_matrix(double ** matrix,int n){
    int i, j;
    printf("\n");
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++)
            printf("%.6f,",matrix[i][j]);
        printf("%s","\n");
    }
    printf("-----\n");
}