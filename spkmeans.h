#ifndef SPK_H_
#define SPK_H_

double geo_c(double z, int n);
double **matrix(int n,int m);
double **rand_matrix(int n);
double **symmetric(double **matrix, int n);
double **Identity(int n);
void print_matrix(double **matrix, int m, int n);
double **adj_matrix(double **matrix, int n, int dim);
double **diag_degree_matrix(double **matrix, int n, int dim);
double dist(double *X, double *Y, int dim);
double **matrix_mult(double **A, double **B, int n);
double **matrix_func(double **matrix, double (*f)(double), int n);
double **L_norm(double **matrix, int n, int dim);
double divsqrt(double A);
double **Jacobi_Rotation_sMatrix(double ** A,int i, int j, int n);
double **Jacobi(double **A, double **V, int n, int iter);
double off(double ** matrix,int n);
void normalize(double **A, int n, int k);
double *diag(double ** matrix,int n);
double **file_to_matrix(char *filename, int *m, int *n);
void error();
int Heuristic (double *list, int n);
double** eigenvectors (double *list, double **V, int n);
int compare_indexes(void *context, const void *b, const void * a);
double **eigen(double **M, int n,int dim);
#endif
