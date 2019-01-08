#include "Python.h"
#include <numpy/arrayobject.h>
#include "kmer_count_multithreads.h" 
#include <atomic>




static PyObject *kmer_count(PyObject *self, PyObject *args)
{
    char* filename;
    int x, K, NumThreads;
    bool Reverse=false;
    if (!PyArg_ParseTuple(args, "siip", &filename, &x, &NumThreads, &Reverse))
        return NULL;
    if (!PyArg_ParseTuple(args, "siip", &filename, &K, &NumThreads, &Reverse))
        return NULL;
    npy_intp SIZE = pow(4, K);
    /*
    std::cout << "K: " << K << std::endl;
    std::cout << "Threads: " << NumThreads << std::endl;
    std::cout << "Reverse: " << Reverse << std::endl;
    */
    static std::vector<std::atomic<int>> count_array;
    count_array = count(filename, K, NumThreads, Reverse);
    return PyArray_SimpleNewFromData(1, &SIZE, NPY_INT32, static_cast<void*>(count_array.data()));
}

static PyObject *kmer_count_seq(PyObject *self, PyObject *args)
{
    char* sequence;
    int x , K, NumThreads;
    bool Reverse=false;
    if (!PyArg_ParseTuple(args, "siip", &sequence, &x, &NumThreads, &Reverse))
        return NULL;
    if (!PyArg_ParseTuple(args, "siip", &sequence, &K, &NumThreads, &Reverse))
        return NULL;
    /*
    std::cout << "K: " << K << std::endl;
    std::cout << "Threads: " << NumThreads << std::endl;
    std::cout << "Reverse: " << Reverse << std::endl;
    */
    npy_intp SIZE = pow(4, K);
    static std::vector<std::atomic<int>> count_array;
    count_array = count_seq(sequence, K, NumThreads, Reverse);
    return PyArray_SimpleNewFromData(1, &SIZE, NPY_INT32, static_cast<void*>(count_array.data()));
}

static PyObject *kmer_count_m_k(PyObject *self, PyObject *args)
{
    char* filename;
    int x, M, K, NumThreads;
    bool Reverse;
    if (!PyArg_ParseTuple(args, "siiip", &filename, &x, &K, &NumThreads, &Reverse))
        return NULL;
    if (!PyArg_ParseTuple(args, "siiip", &filename, &M, &K, &NumThreads, &Reverse))
        return NULL;
    /*
    std::cout << "M: " << M << std::endl;
    std::cout << "K: " << K << std::endl;
    std::cout << "Threads: " << NumThreads << std::endl;
    std::cout << "Reverse: " << Reverse << std::endl;
    */
    npy_intp SIZE = pow(4, K) + pow(4, M); 
    static std::vector<std::atomic<int>> count_array;
    count_array = count_M_K(filename, M, K, NumThreads, Reverse);
    return PyArray_SimpleNewFromData(1, &SIZE, NPY_INT32, static_cast<void*>(count_array.data()));
}

static PyObject *kmer_count_m_k_seq(PyObject *self, PyObject *args)
{
    char* sequence;
    int x, M, K, NumThreads;
    bool Reverse;
    if (!PyArg_ParseTuple(args, "siiip", &sequence, &x, &K, &NumThreads, &Reverse))
        return NULL;
    if (!PyArg_ParseTuple(args, "siiip", &sequence, &M, &K, &NumThreads, &Reverse))
        return NULL;
    /*
    std::cout << "M: " << M << std::endl;
    std::cout << "K: " << K << std::endl;
    std::cout << "Threads: " << NumThreads << std::endl;
    std::cout << "Reverse: " << Reverse << std::endl;
    */
    npy_intp SIZE = pow(4, K) + pow(4, M);
    static std::vector<std::atomic<int>> count_array;
    count_array = count_M_K_seq(sequence, M, K, NumThreads, Reverse);
    return PyArray_SimpleNewFromData(1, &SIZE, NPY_INT32, static_cast<void*>(count_array.data()));
}

static PyMethodDef module_methods[] = {
    {"kmer_count_m_k", kmer_count_m_k, METH_VARARGS, ""},
    {"kmer_count_m_k_seq", kmer_count_m_k_seq, METH_VARARGS, ""},
    {"kmer_count", kmer_count, METH_VARARGS, ""},
    {"kmer_count_seq", kmer_count_seq, METH_VARARGS, ""},
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

