#include <Python.h>
#include "structmember.h"
#include "string.h"             /* for NULL pointers */
#include "mlf.hpp"
#include <new>

// #define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
// #include    <numpy/arrayobject.h>

#if PY_MAJOR_VERSION >= 3
/* see http://python3porting.com/cextensions.html */
    #define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)
#else
    #define MOD_INIT(name) PyMODINIT_FUNC init##name(void)
#endif

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C API functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

PyDoc_STRVAR(moody__doc__,
              "Do Position Weight Matrix stuff\n");


charArray convertSequence(const char *sequence) {
    charArray c_seq;
    int lenght = strlen(sequence);
    for(int i = 0; i < lenght; i++) {
        char toput = 5;
        switch (sequence[i])
        {
            case 'a':
            case 'A': toput = 0; break;
            case 'c':
            case 'C': toput = 1; break;
            case 'g':
            case 'G': toput = 2; break;
            case 't':
            case 'T': toput = 3; break;
            case 'n':
            case 'N': toput = (int) (4 * rand() / (1.0 + RAND_MAX)); break;
            default:
                break;
        }
        if(toput != 5) {
            c_seq.push_back(toput);
        }
    }
    return c_seq;
}



PyObject *atoPyArray(scoreArray a) {
    PyObject *new_row = PyList_New(a.size());
    int il = a.size();
    for(int i=0; i < il; i++) {
        PyList_SET_ITEM(new_row, i, PyFloat_FromDouble(a[i]));
    }
    return new_row;
}

PyObject *atoPyMatrix(scoreMatrix m) {
    PyObject *py_matrix = PyList_New(m.size());
    int il=m.size();
    for(int i = 0; i < il; i++) {
        PyList_SET_ITEM(py_matrix, i, atoPyArray(m[i]));
    }
    return py_matrix;
}

scoreArray atoDoubleArray(PyObject *o) {
    scoreArray t;
    if(!PyList_Check(o)) {
        return t;
    }
    Py_ssize_t length = PyList_Size(o);
    PyObject *tmp;
    for(int i=0; i < length; i++) {
        tmp = PyList_GET_ITEM(o, i);
        t.push_back(PyFloat_AsDouble(tmp));
    }
    return t;
}

scoreMatrix atoDoubleMatrix(PyObject *o) {
    scoreMatrix t;
    if(!PyList_Check(o)) {
        return t;
    }
    Py_ssize_t length = PyList_Size(o);
    PyObject *tmp;
    for(int i=0; i < length; i++) {
        tmp = PyList_GET_ITEM(o, i);
        t.push_back(atoDoubleArray(tmp));
    }
    return t;
}

typedef struct {
    PyObject_HEAD
    MOODS_MLF *mlf;
    int both_strands;
    int num_matrices;
} MOODSSearch;

static void
MOODSSearch_dealloc(MOODSSearch* self) {
    delete self->mlf;
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject *
MOODSSearch_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    MOODSSearch *self;
    self = (MOODSSearch *)type->tp_alloc(type, 0);
    self->mlf = NULL;
    /* allocate other fields later.
    */
    return (PyObject *)self;
}

static int
MOODSSearch_init(MOODSSearch *self, PyObject *args, PyObject *kwds) {
    PyObject *py_matrices;
    PyObject *py_thresholds;
    int q;
    PyObject *py_absolute_threshold;
    PyObject *py_bg;
    PyObject *py_both_strands;
    int absolute_threshold;
    double ps = 0.1;
    int both_strands;

    if (self == NULL) {
        return -1;
    }

    if (!PyArg_ParseTuple(args, "OOOiOO", &py_matrices, &py_thresholds, &py_bg, &q, &py_absolute_threshold, &py_both_strands)) {
        return -1;
    }

    // checkout nothrow here
    // http://www.informit.com/guides/content.aspx?g=cplusplus&seqNum=170
    self->mlf = new(std::nothrow) MOODS_MLF(q);
    if (!self->mlf) {
        return -1;
    }
    MOODS_MLF *mlf_p = self->mlf;

    absolute_threshold = PyObject_IsTrue(py_absolute_threshold);
    mlf_p->thresholds = atoDoubleArray(py_thresholds);
    
    self->both_strands = both_strands = PyObject_IsTrue(py_both_strands);

    mlf_p->bg = atoDoubleArray(py_bg);

    if(!PyList_Check(py_matrices)) {
        return -1;
    }
    unsigned int num_matrices =  PyList_Size(py_matrices);
    if(num_matrices != mlf_p->thresholds.size()) {
        PyErr_SetString(PyExc_RuntimeError, "Thresholds should be as many as matrices");
        return -1;
    }
    self->num_matrices = num_matrices;

    for(unsigned int i=0; i < num_matrices; i++) {
        mlf_p->matrices.push_back(atoDoubleMatrix(PyList_GET_ITEM(py_matrices, i)));
        if(mlf_p->matrices[i].size() != 4) {
            PyErr_SetString(PyExc_RuntimeError, "Matrix size must be 4");
            return -1;
        }
    }

    //Check if parameter parsing has raised an exception
    if(PyErr_Occurred()) {
        return -1;
    }

    if(both_strands) {
        for(unsigned int i=0; i < num_matrices; i++) {
            mlf_p->matrices.push_back(reverseComplement(mlf_p->matrices[i]));
            mlf_p->thresholds.push_back(mlf_p->thresholds[i]);
        }
    }
    if(!absolute_threshold) {
        int il = mlf_p->matrices.size();
        for(int i=0; i < il; i++) {
            mlf_p->matrices[i] = counts2LogOdds(mlf_p->matrices[i], mlf_p->bg, ps);
            mlf_p->thresholds[i] = tresholdFromP(mlf_p->matrices[i], mlf_p->bg, mlf_p->thresholds[i]);
        }
    }

    if (mlf_p->multipleMatrixLookaheadFiltrationDNASetup() < 0) {
        return -1;
    };

    return 0;
};

static PyObject *
MOODSSearch_search(MOODSSearch* self, PyObject *args) {
    const char *sequence;
    std::vector<matchArray> matches;
    charArray c_seq;
    int rc = 0;

    if (!PyArg_ParseTuple(args, "s", &sequence)) {
        return NULL;
    }
    c_seq = convertSequence(sequence);
    matches = self->mlf->doScan(c_seq, &rc);
    if (rc < 0) {
        // cleanup or do goto
        return NULL;
    }

    unsigned int num_matrices = self->num_matrices;
    if(self->both_strands) {
        if(matches.size() != (2 * num_matrices) ) {
            PyErr_SetString(PyExc_RuntimeError, "Unknown error");
            return NULL;
        }
        for(unsigned int i=0; i < num_matrices; i++) {
            while(!matches[num_matrices + i].empty()) {
                matches[num_matrices + i].back().position = -matches[num_matrices + i].back().position;
                matches[i].push_back(matches[num_matrices + i].back());
                matches[num_matrices + i].pop_back();
            }
        }
    }
    PyObject *results = PyList_New(matches.size());
    for(unsigned int i = 0; i < matches.size(); i++) {
        PyObject *new_match_list = PyList_New(matches[i].size());
        for(unsigned int j=0; j < matches[i].size(); j++) {
            PyList_SET_ITEM(new_match_list, j, Py_BuildValue("Ld", matches[i][j].position, matches[i][j].score));
        }
        PyList_SET_ITEM(results, i, new_match_list);
    }

    return results;
}

static PyMemberDef MOODSSearch_members[] = {
    {NULL}  /* Sentinel */
};

static PyGetSetDef MOODSSearch_getsetters[] = {
    {NULL}  /* Sentinel */
};

static PyMethodDef MOODSSearch_methods[] = {
    {"search", (PyCFunction)MOODSSearch_search, METH_VARARGS,
     "do a Lookahead Filtration Search\n"
    },
    {NULL}  /* Sentinel */
};

static PyTypeObject MOODSSearchType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "moody.MOODSSearch",        /*tp_name*/
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

static PyMethodDef moody_mod_methods[] = {
    {NULL}
};

extern "C" {
MOD_INIT(moody) {
    if (PyType_Ready(&MOODSSearchType) < 0) {
        return NULL;
    }
    #if PY_MAJOR_VERSION >= 3
        static struct PyModuleDef moduledef = {
            PyModuleDef_HEAD_INIT,
            "moody",                /* m_name */
            moody__doc__,           /* m_doc */
            -1,                     /* m_size */
            moody_mod_methods,      /* m_methods */
            NULL,                   /* m_reload */
            NULL,                   /* m_traverse */
            NULL,                   /* m_clear */
            NULL,                   /* m_free */
        };
        PyObject* m = PyModule_Create(&moduledef);
        // import_array();
        if (m == NULL) { return NULL; }

        Py_INCREF(&MOODSSearchType);
        PyModule_AddObject(m, "MOODSSearch", (PyObject *)&MOODSSearchType);
        return m;
    #else
        PyObject* m = Py_InitModule3("moody", moody_mod_methods, moody__doc__);
        if (m == NULL) { return; }
        // import_array();
        Py_INCREF(&MOODSSearchType);
        PyModule_AddObject(m, "MOODSSearch", (PyObject *)&MOODSSearchType);
    #endif
};
}