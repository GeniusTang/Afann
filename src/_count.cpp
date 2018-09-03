#include "Python.h"
#include <numpy/arrayobject.h>
#include "kmer_count_multithreads.h" 
#include <atomic>




static PyObject *kmer_count(PyObject *self, PyObject *args)
{
    char* filename;
    int NumThreads=1;
    bool Reverse=false;
    int K=0;
    if (!PyArg_ParseTuple(args, "siii", &filename, &K, &NumThreads, &Reverse))
        return NULL;
    npy_intp SIZE = pow(4, K);
    static std::vector<std::atomic<int>> count_array;
    count_array = count(filename, K, NumThreads, Reverse);
    return PyArray_SimpleNewFromData(1, &SIZE, NPY_INT32, static_cast<void*>(count_array.data()));
}

static PyObject *kmer_count_m_k(PyObject *self, PyObject *args)
{
    char* filename;
    int K=0;
    int M=0;
    int NumThreads=1;
    bool Reverse=false;
    if (!PyArg_ParseTuple(args, "siiii", &filename, &M, &K, &NumThreads, &Reverse))
        return NULL;
    npy_intp SIZE = pow(4, K) + pow(4, M); 
    static std::vector<std::atomic<int>> count_array;
    count_array = count_M_K(filename, M, K, NumThreads, Reverse);
    return PyArray_SimpleNewFromData(1, &SIZE, NPY_INT32, static_cast<void*>(count_array.data()));
}

static PyMethodDef module_methods[] = {
    {"kmer_count_m_k", kmer_count_m_k, METH_VARARGS, ""},
    {"kmer_count", kmer_count, METH_VARARGS, ""},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef countmodule = {
    PyModuleDef_HEAD_INIT,
    "_count",   /* name of module */
    "", /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    module_methods
};

PyMODINIT_FUNC PyInit__count(void)
{
    import_array();
    return PyModule_Create(&countmodule);
}

