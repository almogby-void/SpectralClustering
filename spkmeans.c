#include <stdio.h>
#include <math.h>
<<<<<<< Updated upstream
#include "spkmeans.h"
=======
#include "spkmeans.h"

int main(int argc, char **argv) {
    int i ,n ,dim;
    double **M, **V;
    double *eigenvalues;
    M = file_to_matrix(argv[argc - 1], &n, &dim);
    if (!strcmp(argv[argc - 2], "wam"))
        print_matrix(adj_matrix(M, n, dim), n, dim);
    else if (!strcmp(argv[argc - 2], "ddg"))
        print_matrix(diag_degree_matrix(M, n, dim), n, dim);
    else if (!strcmp(argv[argc - 2], "lnorm"))
        print_matrix(L_norm(M, n, dim), n, dim);
    else if (!strcmp(argv[argc - 2], "lnorm"))
    {
        V = Identity(n);
        eigenvalues =  diag(Jacobi(M, V, n, dim), n);
        for (i = 0; i < n - 1; i++) printf("%.4f,", eigenvalues[i]);
        printf("%.4f\n", eigenvalues[n - 1]);
        print_matrix(V, n, dim);
    }
}
>>>>>>> Stashed changes
