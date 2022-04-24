#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>

#include "spkmeans.h"

#define epsilon 0.00001
#define max_rotations 100

static double *temp_eigen;

int main(int argc, char **argv) {
    int i,n,dim;
    double **M, **V;
    double *eigenvalues;
    M = file_to_matrix(argv[argc-1],&n,&dim);
    if (!strcmp(argv[argc - 2],"wam"))
        print_matrix(adj_matrix(M,n,dim),n,n);
    else if (!strcmp(argv[argc - 2],"ddg"))
        print_matrix(diag_degree_matrix(M,n,dim),n,n);
    else if (!strcmp(argv[argc - 2],"lnorm"))
        print_matrix(L_norm(M,n,dim),n,n);
    else if (!strcmp(argv[argc - 2],"jacobi"))
    {
        V = Identity(n);
        eigenvalues =  diag(Jacobi(L_norm(M,n,dim),V,n,dim),n);
        for (i = 0; i < n-1; i++) printf("%.4f,",eigenvalues[i]);
        printf("%.4f\n",eigenvalues[n-1]);
        print_matrix(V,n,n);
    }
    else
        input_error();
    return 0;
}


double **rand_matrix(int n) {
    int i, j;
    double **M = matrix(n,n);
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            M[i][j] = (rand()) % 10;
    return M;
}

double **symmetric(double **matrix, int n) { /* TODO */ 
    int i, j;
    for (i = 0; i < n; i++)
        for (j = 0; j < i + 1; j++)
            matrix[i][j] = matrix[j][i];
    return matrix;
}

double **adj_matrix(double ** M, int n, int dim){
    int i, j;   
    double **adj;
    adj = matrix(n,n);
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            adj[i][j] = (i != j) * exp(-dist(M[i], M[j], dim) / 2);
    return adj;
}


double **diag_degree_matrix(double **M,int n, int dim) {
    int i, j;
    double **diag;
    diag = matrix(n,n);
 
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            diag[i][i] +=  (i != j) * exp(-dist(M[i],M[j],dim)/2);
        }
    }
    
    return diag;
}

double *diag(double ** matrix,int n){
    int i;
    double *diag = calloc(n,sizeof(double));
    if (diag == NULL)
        error();
    for (i = 0; i < n; i++)
    {
        diag[i] = matrix[i][i];
    }
    return diag;
}

double **Identity(int n) { /* identity matrix */
    int i;
    double **m;
    /* init matrix*/
    m = matrix(n,n);
    
    for(i = 0; i < n; i++)
        m[i][i] = 1;
    return m;
}

/* Initializes a matrix */
double **matrix(int n,int m){
    int i = 0;
    double **M = calloc(n,sizeof(double*)); 
    for(i = 0; i < n; i++)
        M[i] = calloc(m,sizeof(double));
    if (M == NULL)
        error();
    return M;
}


double **L_norm(double **M,int n,int dim) {
    int i, j;
    double **L, **W, **D;
    L = matrix(n,n);
    W = adj_matrix(M,n,dim);
    D = matrix_func(diag_degree_matrix(M,n,dim), divsqrt,n); /* divsqrt is a function */
    L = matrix_mult(matrix_mult(D, W,n), D,n);
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            L[i][j] = (i == j) - L[i][j];
    return L;
}

double **Jacobi_Rotation_Matrix(double **A, int i, int j, int n) { /* Unused? */
    double **m;
    double c, theta, t;
    m = Identity(n);
    theta = (A[i][i] - A[j][j]) / (2 * A[i][j]);
    t = ((theta >= 0) - (theta < 0)) / (fabs(theta) + sqrt(theta * theta + 1));
    c = 1 / sqrt(t * t + 1);
    m[i][i] = c;
    m[j][j] = c;
    m[j][i] = ((i > j) - (i < j)) * t * c;
    m[i][j] = -m[j][i];
    return m;
}


/* Jacobi algorithm, V should be initialized as the identity matrix and iter should be 0 */
double **Jacobi(double **A, double ** V, int n, int iter) {
    int i, j;
    int row, col;
    double **temp; /* A', temp[0]=A'[i], temp[1]=A'[j] */
    double c, theta, t, s, square_diff, old_off;
    temp = calloc(2, sizeof(double*));
    if (temp == NULL)
        error();
    temp[0] = malloc(sizeof(double) * n);
    temp[1] = malloc(sizeof(double) * n);
    if ((temp[0] == NULL) | (temp[1] == NULL))
        error();
    i = 1;
    j = 0;
    old_off = off(A, n);
    for (row = 0; row < n; row++)
        for (col = 0; col < n; col++)
            if ((row != col) & (fabs(A[row][col]) > fabs(A[i][j]))) { 
                i = row;
                j = col;
            }
    theta = (A[j][j] - A[i][i]) / (2 * A[i][j]);
    t = ((theta >= 0) - (theta < 0)) / (fabs(theta) + sqrt(theta * theta + 1));
    c = 1 / sqrt(t * t + 1);
    /* printf("%f,%f,%f,%d,%d,%f\n",theta,t,c,i,j,A[i][j]); */
    s = t * c;
    for (col = 0; col < n; col++) { /* switcheroo'd from the instructions */
        temp[0][col] = c * A[i][col] - s * A[j][col];
        temp[1][col] = c * A[j][col] + s * A[i][col];
    }
    temp[0][i] = c * c * A[i][i] + s * s * A[j][j] - 2 * s * c * A[i][j];
    temp[1][j] = s * s * A[i][i] + c * c * A[j][j] + 2 * s * c * A[i][j];
    temp[0][j] = (c * c - s * s) * A[i][j] + s * c * (A[i][i] - A[j][j]);
    temp[1][i] = temp[0][j];
    for (col = 0; col < n; col++)
        if ((col != i) & (col != j)) {
            A[i][col] = temp[0][col];
            A[col][i] = temp[0][col];
            A[j][col] = temp[1][col];
            A[col][j] = temp[1][col];
        }
    A[i][i] = temp[0][i];
    A[j][j] = temp[1][j];
    A[i][j] = temp[0][j];
    A[j][i] = A[i][j];
    
    square_diff = old_off - off(A,n);
    /*printf("theta:%f,t:%f,c:%f,s:%f,old_off:%f,off:%f,diff:%f,i=%d,j=%d\n",theta,t,c,s,old_off,off(A,n),square_diff,i,j);*/
    for (row = 0; row < n; row++)
    {
        temp[0][row] = c*V[row][i] - s*V[row][j];
        temp[1][row] = s*V[row][i] + c*V[row][j];
    }
    for (row = 0; row < n; row++)
    {
        V[row][i] = temp[0][row];
        V[row][j] = temp[1][row];
    }
    free_matrix(temp);
    if ((square_diff < epsilon) && (++iter > max_rotations)) /* Check why square_diff is so inaccurate */
        return A;
    return Jacobi(A, V,n, iter);
}

int cmpfunc (const void * a, const void * b) {
    return (*(double*)a > *(double*)b) - (*(double*)a < *(double*)b);
}


int Heuristic (double *list, int n){
    int i,argmax;
    double delta,max_delta = 0;
    qsort(list,n,sizeof(double),cmpfunc);

    for (i = 0; i < (int)floor(n/2); i++){
        
        delta = fabs(list[i] - list[i+1]);
        if (delta > max_delta){
            max_delta = delta;
            argmax = i;
        }
    }
    return argmax;
}



/* Returns eigenvectors, V should be the rotation matrix*/
double** eigenvectors (double *list, double **V, int n){
    int i,j,argmax;
    double delta,max_delta = 0;
    double ** U;
    int * indices = calloc(n,sizeof(int));
    if (indices == NULL)
        error();
    for (i = 0; i < n; i++) indices[i] = i; 
    temp_eigen = list;
    qsort(indices,n,sizeof(int),compare_indexes); /* qsort_s helps us sort by indices and therefore sort the eigenvectors, not just eigenvalues */
    for (i = 0; i < (int)floor(n/2); i++){
        delta = fabs(list[indices[i]] - list[indices[i+1]]);
        /*printf("%f,%f,  %f\n",list[indices[i]],list[indices[i+1]],delta);*/
        if (delta > max_delta){
            max_delta = delta;
            argmax = i;
        }
    }
    U = matrix(n,argmax+1);
    for (i = 0; i < n; i++){
        for (j = 0; j < argmax+1; j++)
            U[i][j] = V[i][indices[j]];
    }
    free(indices);
    return U;
}

/* compare indexes for index sort of the eigenvalues */
int compare_indexes(const void *b, const void * a)
{

    return (temp_eigen[*(int*)a] < temp_eigen[*(int*)b]) - (temp_eigen[*(int*)a] > temp_eigen[*(int*)b]);
}

/* normalize for matrix T */
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
void free_matrix(double **array) {
    free(*array);
    free(array);
}

double off(double **matrix, int n) {
    double sum = 0;
    int i, j;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            sum += (i != j) * matrix[i][j] * matrix[i][j];
    return sum;
}

double **matrix_func(double **M, double (*f)(double), int n) {
    int i, j;
    double **result;
    result = matrix(n,n);
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            result[i][j] = f(M[i][j]);
    return result;
}

double divsqrt(double A) {return A ? 1 / sqrt(A) : 0;} /* Divides and does square root for the laplacian function*/

double **matrix_mult(double **A, double **B,int n) {
    int i, j, k; 
    double **result = matrix(n,n);
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            for (k = 0; k < n; k++)
                result[i][j] += A[i][k] * B[k][j];
    return result;
}


double dist(double *X, double *Y, int dim) {
    double result = 0;
    int i = 0;
    for (; i < dim; i++)
        result += (X[i] - Y[i]) * (X[i] - Y[i]);
    return sqrt(result);
}

void error(){
    printf("An Error Has Occurred");
    exit(0);
}

void input_error(){
    printf("Invalid Input!");
    exit(0);
}

void print_matrix(double **matrix, int m, int n){
    int i, j;
    printf("\n");
    for (i = 0; i < m; i++) {
        for (j = 0; j < n-1; j++){
            printf("%.4f,", matrix[i][j]); 
        }
        printf("%.4f\n", matrix[i][n-1]); 
    }
}



typedef struct list_t {
    double *arr;
    struct list_t *next;
} LIST;

double **file_to_matrix(char *filename, int *m, int *n) { /* m rows and n columns */
    int i, rows = 0, cols = 0, line_size = 10;
    char *line, *token;
    double **mat;
    LIST *head, *curr;
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        printf("Invalid Input!\n");
        exit(1);
    }

    line = malloc(sizeof(char) * line_size);
    head = malloc(sizeof(LIST*));
    curr = head;
    fgets(line, line_size, file);
    if (line == NULL)
        exit(1);
    rows++;
    while (1 + strlen(line) == (unsigned) line_size) {
        line = realloc(line, line_size * 2);
        fgets(line + line_size - 1, line_size + 1, file);
            line_size *= 2;
    }
    curr->arr = malloc(sizeof(double) * 10);
    token = strtok(line, ",");
    for (i = 0; i < 10; i++) {
        if (token == NULL)
            break;
        curr->arr[i] = strtod(token, NULL);
        cols++;
        token = strtok(NULL, ",");
    }
    curr->arr = realloc(curr->arr, sizeof(double) * cols);
    while (fgets(line, line_size, file) != NULL) {
        curr->next = malloc(sizeof(LIST*));
        curr = curr->next;
        curr->arr = malloc(sizeof(double) * cols);
        token = strtok(line, ",");
        for (i = 0; i < cols; i++) {
            if (token == NULL)
                break;
            curr->arr[i] = strtod(token, NULL);
            token = strtok(NULL, ",");
        }
        rows++;
    }
    fclose(file);
    mat = malloc(sizeof(double*) * rows);
    for (i = 0; i < rows; i++) {
        mat[i] = head->arr;
        curr = head->next;
        free(head);
        head = curr;
    }
    *m = rows;
    *n = cols;
    return mat;
}