#include <Python.h>
#include "structmember.h"
#include <string.h> /* for NULL pointers */
#include <vector>
#include "pssm_algorithms.hpp"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include    <numpy/arrayobject.h>

#if PY_MAJOR_VERSION >= 3
/* see http://python3porting.com/cextensions.html */
    #define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)
#else
    #define MOD_INIT(name) PyMODINIT_FUNC init##name(void)
#endif

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C API functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

PyDoc_STRVAR(moodsclass__doc__,
              "Do Position Weight Matrix stuff\n");

scoreArray* atoDoubleArray(PyObject *o) {
    if(!PyList_Check(o)) {
        return NULL;
    }

    Py_ssize_t length = PyList_Size(o);
    scoreArray *tp = new scoreArray((int) length);
    scoreArray &t = *tp;
    PyObject *tmp;
    for(int i=0; i< length; i++) {
        tmp = PyList_GET_ITEM(o, i);
        t.push_back(PyFloat_AsDouble(tmp));
    }
    return tp;
}

typedef struct {
    PyObject_HEAD
    int q; 
    std::vector<scoreMatrix> *matrices;
    std::vector<std::vector< OutputListElementMulti> > *output; 
    intArray *window_positions;
    intArray *m; 
    intMatrix *orders; 
    scoreMatrix *L;
    scoreArray *thresholds;
} MOODSSearch;

static void
MOODSSearch_dealloc(MOODSSearch* self) {
    delete self->matrices;
    delete self->output;
    delete self->window_positions;
    delete self->m;
    delete self->orders;
    delete self->L;
    delete self->thresholds;
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject *
MOODSSearch_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    MOODSSearch *self;
    self = (MOODSSearch *)type->tp_alloc(type, 0);
    self->matrices = new std::vector<scoreMatrix>();
    self->output = NULL;
    self->window_positions=NULL;
    self->m = NULL;
    self->orders = NULL;
    self->L = NULL;
    self->thresholds = NULL;
    /* allocate other fields later.
    */
    return (PyObject *)self;
}

static int
MOODSSearch_init(MOODSSearch *self, PyObject *args, PyObject *kwds) {
    PyObject *py_matrices;
    PyObject *py_thresholds;
    const char *algorithm;
    int q;
    PyObject *py_absolute_threshold;
    PyObject *py_bg;
    PyObject *py_combine;
    PyObject *py_both_strands;
    bool absolute_threshold;
    double ps = 0.1;
    bool combine = false;
    bool both_strands;

    std::vector<scoreMatrix> &matrices = self->matrices;

    scoreArray bg;

    if (self == NULL) {
        return -1;
    }

    if (!PyArg_ParseTuple(args, "OOOsiOOO", &py_matrices, &py_thresholds, &py_bg, &algorithm, &q, &py_absolute_threshold, &py_combine, &py_both_strands)) {
        return -1;
    }

    absolute_threshold = (bool) PyObject_IsTrue(py_absolute_threshold);
    combine = (bool) PyObject_IsTrue(py_combine);
    self->thresholds = atoDoubleArray(py_thresholds);
    scoreArray &thresholds = *(self->thresholds);  // Hopefully this works
    
    both_strands = (bool) PyObject_IsTrue(py_both_strands);

    if(py_bg != Py_None) {
        bg = *(atoDoubleArray(py_bg));
    }
    else {
        if(!absolute_threshold) {
            bg = bgFromSequence(c_seq, 4, ps);
        } else {
            bg = flatBG(4);
        }
    }

    if(!PyList_Check(py_matrices)) {
        return NULL;
    }
    int num_matrices = (int) PyList_Size(py_matrices);
    if(num_matrices != thresholds.size()) {
        PyErr_SetString(PyExc_RuntimeError, "Thresholds should be as many as matrices");
        return NULL;
    }

    for(int i=0; i< num_matrices; i++) {
        matrices.push_back(atoDoubleMatrix(PyList_GET_ITEM(py_matrices, i)));
        if(matrices[i].size() != 4) {
            PyErr_SetString(PyExc_RuntimeError, "Matrix size must be 4");
            return NULL;
        }
    }

    //Check if parameter parsing has raised an exception
    if(PyErr_Occurred()) {
        return NULL;
    }

    if(both_strands) {
        for(int i=0; i < num_matrices; i++) {
            matrices.push_back(reverseComplement(matrices[i]));
            thresholds.push_back(thresholds[i]);
        }
    }
    if(!absolute_threshold) {
        for(int i=0; i< matrices.size(); i++) {
            matrices[i] = counts2LogOdds(matrices[i], bg, ps);
            thresholds[i] = tresholdFromP(matrices[i], bg, thresholds[i]);
        }
    }

    if(matrices.size() == 1) {
        combine = false;
    }

    if(strcmp(algorithm, "lf") == 0) {
        const bits_t size = 1 << (BITSHIFT * q); // numA^q
        self->output = new std::vector<std::vector< OutputListElementMulti> >(size);
        self->window_positions = new intArray(matrices.size());
        self->m = new intArray(matrices.size(), 0);
        self->orders = new intMatrix(matrices.size());
        self->L = new scoreMatrix(matrices.size());
        multipleMatrixLookaheadFiltrationDNASetup(q, 
            self->matrices, self->output, self->window_positions, 
            self->m, self->orders, self->L,
            bg, self->thresholds);
    }

    self.q = q;

    return 0;
};

static PyMemberDef MOODSSearch_members[] = {
    {NULL}  /* Sentinel */
};

static PyGetSetDef MOODSSearch_getsetters[] = {
    {NULL}  /* Sentinel */
};

static PyMethodDef MOODSSearch_methods[] = {
    {"checkstuff", (PyCFunction)MOODSSearch_checkstuff, METH_VARARGS,
     "cs\n"
    },
    {NULL}  /* Sentinel */
};

static PyTypeObject MOODSSearchType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "moodsclass.MOODSSearch",        /*tp_name*/
    sizeof(MOODSSearch),            /*tp_basicsize*/
    0,                              /*tp_itemsize*/
    (destructor)MOODSSearch_dealloc,/*tp_dealloc*/
    0,                              /*tp_print*/
    0,                              /*tp_getattr*/
    0,                              /*tp_setattr*/
    0,                              /*tp_compare*/
    0,                              /*tp_repr*/
    0,                              /*tp_as_number*/
    0,                              /*tp_as_sequence*/
    0,                              /*tp_as_mapping*/
    0,                              /*tp_hash */
    0,                              /*tp_call*/
    0,                              /*tp_str*/
    0,                              /*tp_getattro*/
    0,                              /*tp_setattro*/
    0,                              /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "MOODSSearch objects",          /* tp_doc */
    0,                              /* tp_traverse */
    0,                              /* tp_clear */
    0,                              /* tp_richcompare */
    0,                              /* tp_weaklistoffset */
    0,                              /* tp_iter */
    0,                              /* tp_iternext */
    MOODSSearch_methods,            /* tp_methods */
    MOODSSearch_members,            /* tp_members */
    MOODSSearch_getsetters,          /* tp_getset */
    0,                              /* tp_base */
    0,                              /* tp_dict */
    0,                              /* tp_descr_get */
    0,                              /* tp_descr_set */
    0,                              /* tp_dictoffset */
    (initproc)MOODSSearch_init,     /* tp_init */
    0,                              /* tp_alloc */
    MOODSSearch_new,                /* tp_new */
};

static PyMethodDef moodsclass_mod_methods[] = {
    {NULL}
};

extern "C" {
MOD_INIT(moodsclass) {
    if (PyType_Ready(&MOODSSearchType) < 0) {
        return NULL;
    }
    #if PY_MAJOR_VERSION >= 3
        static struct PyModuleDef moduledef = {
            PyModuleDef_HEAD_INIT,
            "moodsclass",            /* m_name */
            moodsclass__doc__,       /* m_doc */
            -1,                     /* m_size */
            moodsclass_mod_methods,  /* m_methods */
            NULL,                   /* m_reload */
            NULL,                   /* m_traverse */
            NULL,                   /* m_clear */
            NULL,                   /* m_free */
        };
        PyObject* m = PyModule_Create(&moduledef);
        import_array();
        if (m == NULL) { return NULL; }

        Py_INCREF(&MOODSSearchType);
        PyModule_AddObject(m, "MOODSSearch", (PyObject *)&MOODSSearchType);
        return m;
    #else
        PyObject* m = Py_InitModule3("moodsclass", moodsclass_methods, moodsclass__doc__);
        if (m == NULL) { return; }
        import_array();
        Py_INCREF(&MOODSSearchType);
        PyModule_AddObject(m, "MOODSSearch", (PyObject *)&MOODSSearchType);
    #endif
};
}