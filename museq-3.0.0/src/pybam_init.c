#include <Python.h>
#include "structmember.h"
#include "fasta.h"
#include "pybam.h"
#include "base.h"

static PyMethodDef pybamMethods[] =
{
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initpybam(void) {
    PyObject *m;
    
    pybam_FastaType.tp_new = PyType_GenericNew;
    pybam_FastaIterType.tp_new = PyType_GenericNew;

    pybam_BamType.tp_new = PyType_GenericNew;
    pybam_BamIterType.tp_new = PyType_GenericNew;

    if (PyType_Ready(&pybam_FastaType) < 0) return;
    if (PyType_Ready(&pybam_FastaIterType) < 0) return;
    if (PyType_Ready(&pybam_BamType) < 0) return;
    if (PyType_Ready(&pybam_BamIterType) < 0) return;

    m = Py_InitModule("pybam", pybamMethods);
    Py_INCREF(&pybam_FastaType);
    Py_INCREF(&pybam_FastaIterType);
    Py_INCREF(&pybam_BamType);
    Py_INCREF(&pybam_BamIterType);

    PyModule_AddObject(m, "Fasta", (PyObject *)&pybam_FastaType);
    PyModule_AddObject(m, "FastaIter", (PyObject *)&pybam_FastaIterType);
    PyModule_AddObject(m, "Bam", (PyObject *)&pybam_BamType);
    PyModule_AddObject(m, "BamIter", (PyObject *)&pybam_BamIterType);

}

