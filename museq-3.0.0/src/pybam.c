#include <Python.h>
#include <sam.h>
#include <stdlib.h>
#include "structmember.h"
#include "base.h"
#include "pybam.h"

PyMethodDef Bam_methods[] = {
    {"vector", (PyCFunction)Bam_counts, METH_VARARGS,
                "fixed size tuples"},
    {"xentropy", (PyCFunction)Bam_xentropy, METH_VARARGS,
                "cross entropy"},
    {NULL}
};

PyMemberDef Bam_members[] = {
    {"targets", T_OBJECT_EX, offsetof(pybam_BamObject, targets), 0, "targets"},
    {"tids", T_OBJECT_EX, offsetof(pybam_BamObject, tids), 0, "tids"},
    {NULL}
};

PyTypeObject pybam_BamType = {
    PyObject_HEAD_INIT(NULL)
    0,"pybam.Bam", sizeof(pybam_BamObject),
    0,
    (destructor)Bam_dealloc,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    "Bam file object",
    0,0,0,0,0,0,
    Bam_methods,
    Bam_members,
    0,0,0,0,0,0,
    (initproc)Bam_init,
    0,
    (newfunc)Bam_new
};

PyTypeObject pybam_BamIterType = {
    PyObject_HEAD_INIT(NULL)
    0,"pybam.BamIter",sizeof(pybam_BamIter),
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_ITER,
    "Bam file iterator object",
    0,0,0,0,
    pybam_BamIter_iter,
    pybam_BamIter_next,
};

pybam_BamIter *
pybam_BamIter_iter(pybam_BamIter *self) {
    // Initialize pileup buffer
    Py_INCREF(self);
    self->buffer = queue_init();
    self->return_value = NULL;
    self->position = 0;
    // Check for the rightmost limit of this contig
    return self;
}

PyTupleObject *
pybam_BamIter_next(pybam_BamIter *self) {
    uint32_t start;
    uint8_t i, j;
    int loop = 1;
    uint32_t stop = 0;
    PyTupleObject *nested_tuple;
    PyTupleObject *tuple;
    column_t column;
    
    if (self->position >= self->stop) {
        queue_destroy(self->buffer);
        Py_DECREF(self);
        PyErr_SetNone(PyExc_StopIteration);
        return NULL;
    }

    while (((self->position >= self->buffer->end) || (self->buffer->size == 0)) &&
           (self->buffer->end <= self->stop)) {
      
        if (self->buffer->position == 0) {
            start = self->start;
        } else {
            start = self->buffer->end;
        }

        stop = start + BUFFER_SIZE;

        if (stop > self->stop) {
            stop = self->stop;
        }
        queue_destroy(self->buffer);
        self->buffer = queue_init();
        self->buffer->fetch_start = start;
        self->buffer->fetch_stop = stop;
        self->pileup = bam_plbuf_init(pileup_func, self);
	//self->pileup->flag_mask &= 772; // disable the BAM_FDUP
	if(self->df){
	  bam_plbuf_set_mask(self->pileup, 772);
	  bam_plp_set_maxcnt(self->pileup->iter, 60000); // increase the max depth size 
	}
        bam_fetch(self->bam->fd->x.bam, self->bam->idx, self->tid,
                  start, stop, (void *)self, fetch_f);
     
        // top off the buffer (as per samtools doc)
        bam_plbuf_push(0, self->pileup);
        bam_plbuf_destroy(self->pileup);
         if ((stop >= self->stop) && (self->buffer->size == 0)) {
            queue_destroy(self->buffer);
            Py_DECREF(self);
            PyErr_SetNone(PyExc_StopIteration);
            return NULL;
        } else if (self->buffer->size == 0) {
            self->buffer->position = stop;
            self->buffer->end = stop;
        }

       loop++;
    }

    column = dequeue(self->buffer);
    self->position = column.position;

    if (self->return_value) {
        Py_DECREF(self->return_value);
    }
    tuple = (PyTupleObject *)PyTuple_New(12);
    PyTuple_SET_ITEM(tuple, 0, PyInt_FromLong((long)column.position));
    for (i = 0; i < 5; i++) {
        nested_tuple = (PyTupleObject *)PyTuple_New(6);
        for (j = 0; j < 6; j++) {
            PyTuple_SET_ITEM(nested_tuple, j, PyFloat_FromDouble((double)column.features[i][j]));
        }
        PyTuple_SET_ITEM(tuple, i + 1, (PyObject *)nested_tuple);
    }
    PyTuple_SET_ITEM(tuple, 6, PyInt_FromLong((long)column.major));
    PyTuple_SET_ITEM(tuple, 7, PyInt_FromLong((long)column.minor));
    PyTuple_SET_ITEM(tuple, 8, PyFloat_FromDouble((double)column.ambiguous));
    PyTuple_SET_ITEM(tuple, 9, PyFloat_FromDouble((double)column.indels));
    PyTuple_SET_ITEM(tuple, 10, PyFloat_FromDouble(column.entropy));
    PyTuple_SET_ITEM(tuple, 11, PyInt_FromLong((long)column.is_del));
    Py_INCREF(tuple);
    self->return_value = (PyObject *)tuple; 
    return tuple;

}

void
Bam_dealloc(pybam_BamObject *self) {
    bam_index_destroy(self->idx);
    samclose(self->fd);
    Py_DECREF(self->tids);
    Py_CLEAR(self->contig);
    self->ob_type->tp_free((PyObject*)self);
}

pybam_BamObject *
Bam_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    pybam_BamObject *self;
    self = (pybam_BamObject *)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->contig = PyString_FromString("");
        if (self->contig == NULL) {
            Py_CLEAR(self);
            return NULL;
        }
        self->fd = NULL;
    }
    return self;
}

int
Bam_init(pybam_BamObject *self, PyObject *args, PyObject *kwds) {
    uint8_t i;
    char *filename = NULL;
    self->header = NULL;
    PyObject *target = NULL;
    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return (int)NULL;
    }

    self->tids = (PyDictObject *)PyDict_New();
    Py_INCREF((PyObject *)self->tids);
    self->fd = samopen(filename, "rb", 0);
    self->idx = bam_index_load(filename);
    self->targets = (PyTupleObject *)PyTuple_New(self->fd->header->n_targets);
    Py_INCREF((PyObject *)self->targets);
    for (i = 0; i < self->fd->header->n_targets; i++) {
        target = Py_BuildValue("s", self->fd->header->target_name[i]);
        PyDict_SetItemString((PyObject *)self->tids,
                    self->fd->header->target_name[i], Py_BuildValue("i", i));
        PyTuple_SET_ITEM(self->targets, i, target);
        
    }
   if (!self->fd) {
        PyErr_SetString(PyExc_IOError, "Cannot open bam file");
        return (int)NULL;
    }
    return 0;
    
}


pybam_BamIter *
Bam_counts(pybam_BamObject *self, PyObject *args) {
    pybam_BamIter *iter;
    int tid, start, stop;
    char *s; 
    uint32_t right_bound; 
    int df; // deep_flag
    start = 0;
    stop = 0;
    df = o;
    if (!PyArg_ParseTuple(args, "s|(ii)i", &s, &start, &stop, &df)) return NULL;
    PyObject *index = PyDict_GetItemString((PyObject *)self->tids, s);
    tid = PyInt_AS_LONG(index);
    if (stop == 0) {
        start = 0;
        stop = self->fd->header->target_len[tid];
    } else {
        start--;
        stop--;
    }
    iter = (pybam_BamIter *)PyObject_New(pybam_BamIter,
                                           &pybam_BamIterType);
    Py_INCREF(iter); 
    right_bound = self->fd->header->target_len[tid];
    if (right_bound < stop) {
        stop = right_bound;
    }

    if (stop < 0) {
        stop = right_bound;
    }

    if (start > stop) {
        return NULL;
    }

    iter->bam = self;
    iter->offset = 0;
    iter->position = start;
    iter->return_value = NULL;
    iter->buffer = NULL;
    iter->start = start;
    iter->tid = tid;
    iter->stop = stop;
    iter->deep_flag = df;
    return iter;
}


PyObject *
Bam_single(pybam_BamObject *self, PyObject *args) {
    char *s;
    uint32_t position;
    if(!PyArg_ParseTuple(arg, "si", &s, &position)) {
        return NULL
    }
    

    return NULL;
}

int
fetch_f(const bam1_t *b, void *data) {
    pybam_BamIter *s = (pybam_BamIter *)data;
    bam_plbuf_t *pileup = (bam_plbuf_t *)s->pileup;
    bam_plbuf_push(b, pileup);
    return 0;
}

int
pileup_func(uint32_t tid, uint32_t pos, int n,
            const bam_pileup1_t *pl, void *data) {
    
    if (n <= 4) {
        return 0;
    }
    int offset, distance, length, r;
    uint8_t quality, mapping, reverse;
    uint8_t base2;
    uint8_t i, j;
    bam1_t *b;
    bam_pileup1_t alignment;
    uint16_t base_counts[4] = {0, 0, 0, 0};
    column_t column;
    for (i = 0; i < 5; i++) {
        for (j = 0; j < 6; j++) {
            column.features[i][j] = 0;
        }
    }
    pybam_BamIter *iterator = (pybam_BamIter *)data;
    queue *buffer = iterator->buffer;
    if ((pos < buffer->fetch_start) || (pos > buffer->fetch_stop)) {
        return 0;
    }
    
    column.is_del = 0;
    column.indels = 0;
    column.ambiguous = 0;
        //
    column.depth = n;
    for (r = 0; r < n; r++) {
        // append tuple to list
        length = 100; 
        //ops = b->core.n_cigar;
        alignment = pl[r];
        b = alignment.b;
        column.indels += alignment.indel;
        column.is_del += alignment.is_del;   
        //offset = pos - b->core.pos;
        offset = alignment.qpos; 
        base2 = base4_base2(bam1_seqi(bam1_seq(b), offset));
        if (base2 <= 3)  {
            base_counts[base2]++;
        } else {
            column.ambiguous++;
            continue;
        }
        reverse = 1 && (b->core.flag & BAM_FREVERSE);
        quality = bam1_qual(b)[offset];
        mapping = b->core.qual;
        distance = 0;
        if (reverse) {
            distance = length - offset;
        } else {
            distance = offset;
        }

        column.features[base2][0] += 1;
        column.features[base2][1] += quality;
        column.features[base2][2] += mapping;
        column.features[base2][3] += distance;
        column.features[base2][4] += reverse;
        column.features[base2][5] += 0;

        column.features[4][0] += 1;
        column.features[4][1] += quality;
        column.features[4][2] += mapping;
        column.features[4][3] += distance;
        column.features[4][4] += reverse;
        column.features[4][5] += 0;
    }
    if (column.features[4][0] == 0) {
        return 0;
    }
    column.major = 0;
    uint16_t max = 0;
    for (i = 0; i < 4; i++) {
        if (base_counts[i] > max) {
            column.major = i;
            max = base_counts[i];
        }
    }
    column.minor = column.major;
    max = 0;
    for (i = 0; i < 4; i++) {
        if ((base_counts[i] > max) && (i != column.major)) {
            column.minor = i;
            max = base_counts[i];
        }
    }

    column.position = (pos + 1); // HERE BE DRAGONS
    column.entropy = entropy(base_counts, column.depth);
    enqueue(buffer, column, column.position);
    buffer->end = column.position;
    return 0;
}

int
enqueue(queue *list, column_t content, uint32_t pos) {
    queue_node *node;
    node = (queue_node *)malloc(sizeof(queue_node));
    if (list->size == 0) {
        list->position = pos;
        list->head = node;
    } else {
        list->tail->next = node;
    }
     
    list->tail = node;
    node->content = content;
    node->position = pos;
    node->next = NULL;
    (list->size)++;
    return 0;
}

column_t
dequeue(queue *list) {
    column_t content;

    queue_node *next;
    content = list->head->content;
    list->position = list->head->position;
    list->size--;
    next = list->head->next;
    free(list->head);
    list->head = next;
    return content;
}

queue*
queue_init(void) {
    queue *buffer;
    buffer = (queue *)malloc(sizeof(queue));
    buffer->size = 0;
    buffer->end = 0;
    buffer->tail = NULL;
    buffer->head = NULL;
    buffer->position = 0;
    return buffer;
}

int
queue_destroy(queue *list) {
    column_t item;
    int count = 0;
    while(list->size > 0) {
        item = dequeue(list);
        count++;
    }
    free(list);
    list = NULL;
    return 0;
}

double
entropy(uint16_t bases[4], uint16_t depth) {
    uint8_t i;
    double e = 0;
    float pr;
    for (i = 0; i < 4; i++) {
        pr = (float)bases[i] / (float)depth;
        if (pr != 0) {
            e += (log(pr) * pr); 
        }
    }
    return -e;
}

PyObject*
Bam_xentropy(pybam_BamObject *self, PyObject *args) {
    PyTupleObject *n_tuple;
    PyTupleObject *t_tuple;
    PyObject *r;
    uint8_t i;
    uint16_t n_depth, t_depth;
    double e;
    float pr_n, pr_t;

    if (!PyArg_ParseTuple(args, "OO", &n_tuple, &t_tuple)) return NULL;
   
    n_depth = (uint16_t)PyInt_AsLong(PyTuple_GET_ITEM(n_tuple, 4));
    t_depth = (uint16_t)PyInt_AsLong(PyTuple_GET_ITEM(t_tuple, 4));
    
    e = 0;
    for (i = 0; i < 4; i++) {
        pr_n = (float)PyInt_AsLong(PyTuple_GET_ITEM(n_tuple, i)) / (float)n_depth;
        pr_t = (float)PyInt_AsLong(PyTuple_GET_ITEM(t_tuple, i)) / (float)t_depth;
        if (pr_t != 0) {
            if (pr_n == 0) {
                e += (-7 * pr_t);
            } else {
                e += (log(pr_n) * pr_t);
            }
        }
    }
    r = PyFloat_FromDouble(-e);
    Py_INCREF(r);
    return r;
}
