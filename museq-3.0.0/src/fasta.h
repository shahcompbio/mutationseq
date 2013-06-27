#ifndef _fasta_h
#define _fasta_h
#include "debug.h"
#include "base.h"
#include <faidx.h>

// Fasta file object
typedef struct {
    PyObject_HEAD
    uint32_t position;
    PyObject *contig;
    PyObject *contigs;
    char *sequence;
    uint16_t counts[4];
    uint32_t length;
    uint32_t l, r, old_position;
    double entropy, gc;
    PyObject *return_value;
    faidx_t *fd;
} pybam_FastaObject;

// Fasta file object functions
PyObject *Fasta_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
int Fasta_init(pybam_FastaObject *self, PyObject *args, PyObject *kwds);
PyObject *Fasta_load(pybam_FastaObject *self, PyObject *args);
PyObject *Fasta_tuple(pybam_FastaObject *self, PyObject *args);
PyObject *Fasta_nmer(pybam_FastaObject *self, PyObject *args);
void Fasta_dealloc(pybam_FastaObject *self);

PyObject *Fasta_enter(PyObject *self);
PyObject *Fasta_exit(PyObject *self);

// Fasta file object implementation details
extern PyMethodDef Fasta_methods[];
extern PyMemberDef Fasta_members[];

PyTypeObject pybam_FastaType;

// The iterator
typedef struct {
    PyObject_HEAD
    uint32_t position;
    uint32_t length;
    uint16_t gc;
    double entropy;
    char *sequence;
    PyObject *base;
    pybam_FastaObject *fasta;
} pybam_FastaIter;

PyObject *pybam_FastaIter_iter(PyObject *self);
PyObject *pybam_FastaIter_next(PyObject *self);

PyTypeObject pybam_FastaIterType;

uint16_t h_f(char *sequence, uint32_t position);
uint16_t h_b(char *sequence, uint32_t position);
uint16_t nmer(char *sequence, uint32_t position, uint8_t n, int i);
double gc_window(char *sequence, uint16_t counts[4]);
double entropy_window(char *sequence, uint16_t counts[4]);
#endif
