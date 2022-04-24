#define PY_SSIZE_T_CLEAN  /* For all # variants of unit formats (s#, y#, etc.) use Py_ssize_t rather than int. */
#include <Python.h>       /* MUST include <Python.h>, this implies inclusion of the following standard headers:
                             <stdio.h>, <string.h>, <errno.h>, <limits.h>, <assert.h> and <stdlib.h> (if available). */
#include <math.h>         /* include <Python.h> has to be before any standard headers are included */
#include "spkmeans.h"

void input_error();
double **matrix(int n,int m);
static PyObject* geo_capi(PyObject *self, PyObject *args);
/*
 * Helper function that will not be exposed (meaning, should be static)
 */

/*
 * A geometric series up to n. sum_up_to_n(z^n)
 */



/*
 * This actually defines the geo function using a wrapper C API function
 * The wrapping function needs a PyObject* self argument.
 * This is a requirement for all functions and methods in the C API.
 * It has input PyObject *args from Python.
 */
static PyObject* goal(PyObject *self, PyObject *args)
{
    int         i,j,n,dim;
    char        *command;
    PyObject    *pList = NULL;
    PyObject    *py_eigenvalues, *py_T;
    PyObject    *pItem,*c;
    double      **M,results;
    double      *eigenvalues;
    /* This parses the Python arguments into a double (d)  variable named z and int (i) variable named n*/
    if(!PyArg_ParseTuple(args, "siiO", &command,&n,&dim,&pList)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }
    if (PyList_Size(pList) != n*dim)
        input_error();

    M = matrix(n,dim);
    for (i = 0; i < n*dim; i++){
        pItem = PyList_GetItem(pList, i);
        M[n%dim][(int)floor(n/i)] = PyFloat_AsDouble(pItem);
    }
    if (command == "spk"){
        eigenvalues = diag(Jacobi(L_norm(M,n,dim),V,n,dim),n);
        results = eigenvectors(eigenvalues,V,n);
        normalize(results,n,Heuristic(results,n));
        py_T = MatrixPyObject(eigenvectors, n, n);
        free_matrix(eigenvectors);
        return Py_BuildValue("(O)", py_T);
    }
    else if (command == "wam"){
        results = adj_matrix(M,n,dim);
        return Py_BuildValue("(O)", MatrixPyObject(results,n,n));
    }
    else if (command = "ddg")
    {
        results = diag_degree_matrix(M,n,dim);
        return Py_BuildValue("(O)", MatrixPyObject(results,n,n));
    }
    else if (command = "lnorm")
    {
        results = L_norm(M,n,dim);
        return Py_BuildValue("(O)", MatrixPyObject(results,n,n));
    }
    else if (command = "jacobi")
    {
        V = Identity(n);
        eigenvalues =  diag(Jacobi(L_norm(M,n,dim),V,n,dim),n);
        results = diag_degree_matrix(M,n,dim);
        return Py_BuildValue("(OO)", MatrixPyObject(results,n,n));
    }
    
}


static PyObject *eigenvectors_capi(PyObject *self, PyObject *args){
    PyObject *matrix,*py_eigenvectors, *py_eigenvalues;
    double **c_matrix,**eigenvectors;
    int i,size;

    if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &matrix) || (!PyList_Check(matrix))){
        printf("Invalid Input!\n");
        return NULL;
    }

    size = (int)PyList_Size(matrix);
    eigenvectors = calculate_eigenvalues(c_matrix, size);
    py_eigenvectors = MatrixPyObject(eigenvectors, size, size);
    py_eigenvalues = PyList_New(size);
    for (i = 0; i < size; i++)
        PyList_SetItem(py_eigenvalues, i, PyFloat_FromDouble(c_matrix[i][i]));

    free_matrix(c_matrix);
    free_matrix(eigenvectors);
    return Py_BuildValue("(OO)", py_eigenvalues, py_eigenvectors);
}

PyObject *MatrixPyObject(double **matrix, int rows, int columns){
    PyObject *row;
    int i, j;
    PyObject *py_matrix = PyList_New(rows);
    for (i = 0; i < rows; i++){
        row = PyList_New(columns);
        for (j = 0; j < columns; j++)
            PyList_SetItem(row, j, PyFloat_FromDouble(matrix[i][j]));
        PyList_SetItem(py_matrix, i, row);
    }
    return py_matrix;
}


/*
 * This array tells Python what methods this module has.
 * We will use it in the next structure
 */
static PyMethodDef capiMethods[] = {
    {"geo",                   /* the Python method name that will be used */
      (PyCFunction) geo_capi, /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,           /* flags indicating parametersaccepted for this function */
      PyDoc_STR("A geometric series up to n. sum_up_to_n(z^n)")}, /*  The docstring for the function */
    {NULL, NULL, 0, NULL}     /* The last entry must be all NULL as shown to act as a
                                 sentinel. Python looks for this entry to know that all
                                 of the functions for the module have been defined. */
};


/* This initiates the module using the above definitions. */
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "spkmeansmodule ", /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,  /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    capiMethods /* the PyMethodDef array from before containing the methods of the extension */
};


/*
 * The PyModuleDef structure, in turn, must be passed to the interpreter in the module’s initialization function.
 * The initialization function must be named PyInit_name(), where name is the name of the module and should match
 * what we wrote in struct PyModuleDef.
 * This should be the only non-static item defined in the module file
 */
PyMODINIT_FUNC
PyInit_spkmeansmodule(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}

void input_error(){
    printf("Invalid Input!");
    exit(0);
}