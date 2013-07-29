#ifndef _pybam_h
#define _pybam_h
#include <sam.h>
#include "debug.h"
#include "base.h"

#define BUFFER_SIZE 1000000 // ~ 200MB

typedef struct {
    PyObject_HEAD
    PyObject *contig;
    PyTupleObject *targets;
    samfile_t *fd;
    bam_index_t *idx;
    bam_header_t *header;
    PyDictObject *tids;
    PyObject *callback;
} pybam_BamObject;

// Bamfile object functions
pybam_BamObject *Bam_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
int Bam_init(pybam_BamObject *self, PyObject *args, PyObject *kwds);
PyObject *Bam_xentropy(pybam_BamObject *self, PyObject *args);
void Bam_dealloc(pybam_BamObject *self);

// Bamfile object implementation details
extern PyMethodDef Bam_methods[];
extern PyMemberDef Bam_members[];
PyTypeObject pybam_BamType;

typedef struct {
    uint32_t position;
// FD: change to resolve overflow bug for deep museq    
//     uint16_t depth;
    uint32_t depth;
// FD: change to resolve overflow bug for deep museq    
//     uint16_t features[5][6];
//     uint16_t features_f[2];
    uint32_t features[5][6];
    uint32_t features_f[2];
    uint8_t major;
    uint8_t minor;
    int indels;
    uint16_t is_del;
    double entropy;
    uint16_t ambiguous;
} column_t;

typedef struct {
    void *next;
    uint32_t position;
    column_t content;
} queue_node;

typedef struct {
    uint32_t size;
    uint32_t fetch_start;
    uint32_t fetch_stop;
    uint32_t end;
    uint32_t position;
    queue_node *head;
    queue_node *tail;
} queue;

// The iterator
typedef struct {
    PyObject_HEAD
    uint32_t position;
    uint16_t offset;
    uint32_t start;
    uint32_t stop;
    int tid;
    pybam_BamObject *bam;
    PyObject *return_value;
    bam_plbuf_t *pileup;
    queue *buffer; // the buffer linked list
    int deep_flag; //J
} pybam_BamIter;


pybam_BamIter *pybam_BamIter_iter(pybam_BamIter *self);
PyTupleObject *pybam_BamIter_next(pybam_BamIter *self);

PyTypeObject pybam_BamIterType;

int fetch_f(const bam1_t *b, void *data);

int pileup_func(uint32_t tid, uint32_t pos, int n,
                       const bam_pileup1_t *pl, void *data);

int enqueue(queue *list, column_t content, uint32_t pos);
column_t dequeue(queue *list);
queue *queue_init(void);
int queue_destroy(queue *list);
pybam_BamIter *Bam_counts(pybam_BamObject *self, PyObject *args);
// FD: change to resolve overflow bug for deep museq
// double entropy(uint16_t bases[4], uint16_t depth);
double entropy(uint32_t bases[4], uint32_t depth);
#endif
