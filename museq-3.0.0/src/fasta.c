#include <Python.h>
#include <faidx.h>
#include "structmember.h"
#include "base.h"
#include "fasta.h"

/* Fasta file iterator functions */
PyMethodDef Fasta_methods[] = {
    {"load", (PyCFunction)Fasta_load, METH_VARARGS, "load chromosome"},
    {"nmer", (PyCFunction)Fasta_nmer, METH_VARARGS, "count tandem repeats"},
    {"vector", (PyCFunction)Fasta_tuple, METH_VARARGS, "fetch tuple"},
    {"__enter__", (PyCFunction)Fasta_enter, METH_VARARGS, "context entry"},
    {"__exit__", (PyCFunction)Fasta_exit, METH_VARARGS, "context exit"},
    {NULL}
};

PyMemberDef Fasta_members[] = {
    {"contig", T_OBJECT_EX, offsetof(pybam_FastaObject, contig), 0, "name"},
    {NULL}
};

PyTypeObject pybam_FastaType = {
    PyObject_HEAD_INIT(NULL)
    0,"pybam.Fasta", sizeof(pybam_FastaObject),
    0,
    (destructor)Fasta_dealloc,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    "Fasta file object",
    0,0,0,0,0,0,
    Fasta_methods,
    Fasta_members,
    0,0,0,0,0,0,
    (initproc)Fasta_init,
    0,
    Fasta_new
};

PyTypeObject pybam_FastaIterType = {
    PyObject_HEAD_INIT(NULL)
    0,"pybam.FastaIter",sizeof(pybam_FastaIter),
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_ITER,
    "Fasta file iterator object",
    0,0,0,0,
    pybam_FastaIter_iter,
    pybam_FastaIter_next,
};

PyObject *
Fasta_enter(PyObject *self)
{
    return self;
}

PyObject *
Fasta_exit(PyObject *self)
{
    return self;
}

PyObject *
pybam_FastaIter_iter(PyObject *self)
{
    Py_INCREF(self);
    return self;
}

PyObject *
pybam_FastaIter_next(PyObject *self)
{
    pybam_FastaIter *s;

    s = (pybam_FastaIter *)self;
    if (!s->base) {
        Py_DECREF(s->base);
    }
    pybam_FastaIter *iter = (pybam_FastaIter *) self;
    if (iter->position < iter->length) {
//        PyObject *t = char_base(iter->sequence[iter->i]);
        iter->position++;
//        return t;
    }
    else {
        free(iter->sequence);
        PyErr_SetNone(PyExc_StopIteration);
    }
    return NULL;
}

void
Fasta_dealloc(pybam_FastaObject *self)
{
    Py_XDECREF(self->contig);
    fai_destroy(self->fd);
    self->ob_type->tp_free((PyObject*)self);
}

PyObject *
Fasta_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    pybam_FastaObject *self;
    self = (pybam_FastaObject *)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->contig = PyString_FromString("");
        if (self->contig == NULL) {
            Py_DECREF(self);
            return NULL;
        }
        self->position = 0;
        self->fd = NULL;
    }
//    self->contigs = Py_
    return (PyObject *)self;
}

int
Fasta_init(pybam_FastaObject *self, PyObject *args, PyObject *kwds)
{
    char *filename = NULL;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return (int)NULL;
    }

    self->return_value = NULL;
    self->sequence = NULL;
    self->fd = fai_load(filename);
    if (!self->fd) {
        PyErr_SetString(PyExc_IOError, "Cannot open fasta file");
        return (int)NULL;
    }
    return 0;
}

PyObject *
Fasta_load(pybam_FastaObject *self, PyObject *args)
{
    int l;
    char *s;

    if (self->sequence != NULL) {
        free(self->sequence);
    }
    if (!PyArg_ParseTuple(args, "s", &s)) return NULL; 

    self->sequence = fai_fetch(self->fd, s, &l);
    self->length = l;
    self->counts[0] = 0;
    self->counts[1] = 0;
    self->counts[2] = 0;
    self->counts[4] = 0;
    self->old_position = 0;
    if (self->sequence == NULL) {
        PyErr_SetString(PyExc_ValueError, "Contig does not exist");
        return NULL;
    }
    
    return Py_None;
}

PyObject *
Fasta_nmer(pybam_FastaObject *self, PyObject *args)
{
    uint32_t position;
    uint16_t f, b;
    uint8_t n;
    int i;
    PyTupleObject *pair;
    if (!PyArg_ParseTuple(args, "ii", &position, &n)) return NULL;
    if (n > 16) {
        return NULL;
    }
    position--;
    pair = (PyTupleObject *)PyTuple_New(2);
    Py_INCREF(pair);
    f = nmer(self->sequence, position, n, 1);
    b = nmer(self->sequence, position, n, -1);
    PyTuple_SET_ITEM(pair, 0, PyInt_FromLong((long)b));
    PyTuple_SET_ITEM(pair, 1, PyInt_FromLong((long)f));
    return pair;
}
PyObject *
Fasta_tuple(pybam_FastaObject *self, PyObject *args)
{
    base2_t b;
    PyTupleObject *tuple; 
    if (self->return_value != NULL) {
        Py_DECREF(self->return_value);
    }
    tuple = (PyTupleObject *)PyTuple_New(5);
    Py_INCREF(tuple);
    uint32_t i, l, r;
    uint32_t pos;
    if (!PyArg_ParseTuple(args, "i", &pos)) return NULL;
    if (pos > self->length)  {
        trace("requesting position outside of contig bounds"); 
        return Py_None;
    }
    pos--; // lower 1-based index to 0-based
    uint16_t x;
    uint16_t y;
    uint16_t range = 500;
    if ((self->old_position == 0) || (abs(self->old_position - pos) > 150)) {
        self->counts[0] = 0;
        self->counts[1] = 0;
        self->counts[2] = 0;
        self->counts[3] = 0;
        if (((int)pos - (range / 2)) < 0) {
            l = 0;
            r = range;
        } else if ((pos + (range / 2)) > self->length) {
            l = self->length - range;
            r = self->length;
        } else {
            l = (pos - (range / 2));
            r = (pos + (range / 2));
        }
         for (i = l; i < r; i++) {
            b = char_base2(self->sequence[i]);
            if (b <= 3) {
                self->counts[b]++;
            }
        }
        self->entropy = entropy_window(self->sequence, self->counts);
        self->gc = gc_window(self->sequence, self->counts);
        self->old_position = pos;
    }

    x = h_f(self->sequence, pos);
    y = h_b(self->sequence, pos);
    PyTuple_SET_ITEM(tuple, 0, PyInt_FromLong(char_base2(self->sequence[pos])));
    PyTuple_SET_ITEM(tuple, 1, PyInt_FromLong(x));
    PyTuple_SET_ITEM(tuple, 2, PyInt_FromLong(y));
    PyTuple_SET_ITEM(tuple, 3, PyFloat_FromDouble(self->gc));
    PyTuple_SET_ITEM(tuple, 4, PyFloat_FromDouble(self->entropy));
    self->return_value = (PyObject *)tuple;
    return (PyObject *)tuple;
}

uint16_t
h_f(char *sequence, uint32_t position) {
    base2_t b;
    base2_t next;
    uint16_t count;
    count = 0;
    b = char_base2(sequence[position + 1]);
    next = char_base2(sequence[position + 2]);
    while ((b == next) && (b <= 3)) {
        count++;
        position++;
        b = char_base2(sequence[position + 1]);
        next = char_base2(sequence[position + 2]);
    }
    return count;

}

uint16_t
h_b(char *sequence, uint32_t position) {
    base2_t b;
    base2_t next;
    uint16_t count;
    count = 0;
    if (position <= 3) {
        return 0;
    }
    b = char_base2(sequence[position - 1]);
    next = char_base2(sequence[position - 2]);
    while ((b == next) && (position >= 2) && (b <= 3)) {
        count++;
        position--;
        b = char_base2(sequence[position - 1]);
        next = char_base2(sequence[position - 2]);
    }
    return count;
}

uint16_t 
nmer(char *sequence, uint32_t position, uint8_t n, int i) {
    ngram32_t b;
    ngram32_t next;
    uint16_t count;
    int offset;
    count = 0;
    b = 0;
    for (offset = 0; offset < n; offset++) {
        ngram32_push(&b, char_base2(sequence[position + (i * offset)]));
//        printf("%c", sequence[position + (i * offset)]);
    }
    next = b;
    while ((b == next) && (position >= n + 1)) {
        count++;
        position += (i * n);

        next = 0;
        for (offset = 0; offset < n; offset++) {
            ngram32_push(&next, char_base2(sequence[position + (i * offset)]));
//            printf("%c", sequence[position + (i * offset)]);
        }
//        printf("^\n");
        if (count > 256) {
            break;
        }
    }
    return count - 1;
}
double gc_window(char *sequence, uint16_t counts[4]) {
    double gc, total;
    gc = counts[1] + counts[2];
    total = counts[0] + counts[1] + counts[2] + counts[3];
    if (total == 0) {
        return 0;
    }
    return gc / total;
}

double entropy_window(char *sequence, uint16_t counts[4]) {
    double entropy, pr, total;
    uint8_t i;
    entropy = 0;
    total = counts[0] + counts[1] + counts[2] + counts[3];
    if (total == 0) {
        return 0;
    }
    for (i = 0; i < 4; i++) {
        pr = (double)counts[i] / total;
        if (pr != 0) {
            entropy += ((log(pr)/log(4)) * pr);
        }
    }
    return -entropy;
}


