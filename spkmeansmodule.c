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
static PyObject* geo_capi(PyObject *self, PyObject *args)
{
    int         i,j,n,m;
    char        *command;
    PyObject    *pList = NULL;
    PyObject    *pResult;
    PyObject    *pItem,*c;
    double      **M;
    /* This parses the Python arguments into a double (d)  variable named z and int (i) variable named n*/
    if(!PyArg_ParseTuple(args, "siiO", &command,&n,&m,&pList)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }
    if (PyList_Size(pList) != n*m)
        input_error();

    M = matrix(n,m);
    for (i = 0; i < n*m; i++)
    {
        pItem = PyList_GetItem(pList, i);
        M[n%m][(int)floor(n/i)] = PyFloat_AsDouble(pItem);
    }
    pResult = PyList_New(0);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            c = Py_BuildValue("f",M[i][j]);
            PyList_Append(pResult,c);
        }
    }
    return pResult;
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
    "capi_project", /* name of module */
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
PyInit_capi_project(void)
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