#define PYSPA_MAX_ARGS 64

#include <Python.h>
#include <numpy/arrayobject.h>
#include "spa.h"

/* Docstrings */
static char module_docstring[] =
    "This module provides an interface for the Solar Position Algorithm (SPA).";
static char   calc_docstring[] =
    "Calculate solar zenith, azimuth, and possibly incidence angles according to inputs.";

/* Available functions */
static PyObject *spa_calc(PyObject *self, PyObject *args);

/* Module specification */
static PyMethodDef module_methods[] = {
    {"calc", spa_calc, METH_VARARGS, calc_docstring},
    {NULL, NULL, 0, NULL}
};

/* Initialize the module */
PyMODINIT_FUNC init_spa(void)
{
    PyObject *m = Py_InitModule3("_spa", module_methods, module_docstring);
    if (m == NULL)
        return;

    /* Load `numpy` functionality. */
    import_array();
}

static PyObject *spa_calc(PyObject *self, PyObject *args)
{
    int N_IN, N_DOUBLE_IN, N_ARRAY_IN, N_DOUBLE_OUT, N_ARRAY_OUT, N;
    
    PyObject *input_obj[PYSPA_MAX_ARGS];
    double    input_double[PYSPA_MAX_ARGS];
    
    PyArrayObject *input_arr[PYSPA_MAX_ARGS],     *output_arr[PYSPA_MAX_ARGS];
    double        *input_arr_ptr[PYSPA_MAX_ARGS], *output_arr_ptr[PYSPA_MAX_ARGS];
    int            input_arr_map[PYSPA_MAX_ARGS], input_double_map[PYSPA_MAX_ARGS], input_type[PYSPA_MAX_ARGS];
    
    double    year=0, month=0, day=0, hour=0, minute=0, second=0, latitude=0, longitude=0, elevation=0, slope=0, aspect=0;
    double    zenith, azimuth, incidence;
    double    value;
    
    int       flag;
    
    int       i, j, ndim, dims[NPY_MAXDIMS], rc=0;
    spa_data  spa;
    
    static const double m_air   = 0.02896; // molecular mass of air, kg mol^-1
    static const double R_const = 8.3143;  // gas constant, N M mol^-1 K^-1
    static const double g_const = 9.807;   // gravity constant, m s^-2
    static const double T_const = 288.15;  // "default" air temperature, K
        
    /* Parse the input tuple */
    for (i=0; i<64; i++) input_obj[i] = NULL;
    if (!PyArg_ParseTuple(args, "OOOOOOOO|OOO", &input_obj[0], &input_obj[1], &input_obj[2], &input_obj[3], &input_obj[4], 
                                                &input_obj[5], &input_obj[6], &input_obj[7], &input_obj[8], &input_obj[9], &input_obj[10]))
        return NULL;
    
    if (input_obj[8] == NULL) {
        // fprintf(stderr, "Running SPA_ZA mode ...\n");
        spa.function = SPA_ZA;
        N_IN         = 8;
        N_DOUBLE_OUT = 0;
        N_ARRAY_OUT  = 2;
    }
    else {
        // fprintf(stderr, "Running SPA_ZA_INC mode ...\n");
        spa.function = SPA_ZA_INC;
        N_IN         = 11;
        N_DOUBLE_OUT = 0;
        N_ARRAY_OUT  = 3;
    }
    
    N_DOUBLE_IN = 0;
    N_ARRAY_IN  = 0;
    for (i=0; i<N_IN; i++) {
        /* Interpret the input objects as numpy arrays. */
        input_arr[i] = (PyArrayObject *) PyArray_FROM_OTF(input_obj[i], NPY_DOUBLE, NPY_IN_ARRAY);
        /* Is is really a numpy array? */
        if (PyArray_NDIM(input_arr[i])==0) {
            input_type[i] = 0;
            input_double_map[N_DOUBLE_IN] = i;
            N_DOUBLE_IN++;
            input_double[i] = PyFloat_AsDouble(input_obj[i]);
            Py_XDECREF(input_arr[i]);
        }
        else {
            input_type[i] = 1;
            input_arr_map[N_ARRAY_IN] = i;
            N_ARRAY_IN++;
        }
    }
    
    /* If that didn't work, throw an exception. */
    flag = 0; for (i=0; i<N_ARRAY_IN; i++) if (input_arr[input_arr_map[i]]==NULL) flag = 1;
    if (flag==1) {
        for (i=0; i<N_ARRAY_IN; i++) Py_XDECREF(input_arr[input_arr_map[i]]);
        return NULL;
    }
    
    /* Get pointers to the data as C-types. */
    for (i=0; i<N_ARRAY_IN; i++) input_arr_ptr[input_arr_map[i]] = (double*) PyArray_DATA(input_arr[input_arr_map[i]]);
    
    /* How many data points are there? */
    ndim = (int)PyArray_NDIM(input_arr[input_arr_map[0]]);
    
    N=1;
    for (j=0; j<ndim; j++) {
         dims[j] = PyArray_DIM(input_arr[input_arr_map[0]], j);
         N *= dims[j];
         /* check for dimension size compatibility */
         flag = 0; for (i=1; i<N_ARRAY_IN; i++) if (dims[j] != PyArray_DIM(input_arr[input_arr_map[i]], j)) flag = 1;
         if (flag==1) {
             for (i=0; i<N_ARRAY_IN; i++) Py_XDECREF(input_arr[input_arr_map[i]]);
             PyErr_SetString(PyExc_RuntimeError, "different dimensions of input arrays.");
             return NULL;
         }
    }
        
    for (i=0; i<N_ARRAY_OUT; i++) {
        output_arr[i]     = (PyArrayObject *) PyArray_FromDims(ndim, dims, NPY_DOUBLE);
        output_arr_ptr[i] = (double*) PyArray_DATA(output_arr[i]);
    }
    
    /* Call the external C function to compute the results. */
    for (j=0; j<N; j++) {
        
        for (i=0; i<N_IN; i++) {
            if (input_type[i] == 0) value = input_double[i];
            else                    value = input_arr_ptr[i][j];
            switch (i) {
                case  0: year      = value;
                case  1: month     = value;
                case  2: day       = value;
                case  3: hour      = value;
                case  4: minute    = value;
                case  5: second    = value;
                case  6: latitude  = value;
                case  7: longitude = value;
                case  8: elevation = value;
                case  9: slope     = value;
                case 10: aspect    = value;
            }
        }
        
        /* some checks */
        if (longitude>180) longitude -= 360;
        
        /* SPA's azm_rotation angle is measured from south and most aspect angle is measured from north */
        aspect -= 180; if (aspect<-360) aspect += 360;
        
        spa.year          = (int) year;
        spa.month         = (int) month;
        spa.day           = (int) day;
        spa.hour          = (int) hour;
        spa.minute        = (int) minute;
        spa.second        = (int) second;
        spa.latitude      = latitude;
        spa.longitude     = longitude;
        spa.timezone      = 0.0;
        spa.delta_t       = 0;
        spa.elevation     = elevation;
        spa.pressure      = 1000*exp(-m_air*g_const*elevation/(R_const*T_const));
        spa.temperature   = 0;
        spa.slope         = slope;
        spa.azm_rotation  = aspect;
        spa.atmos_refract = 0.5667;
        
        rc = spa_calculate(&spa);
        
        if (rc == 0) {
            zenith    = spa.zenith;
            azimuth   = spa.azimuth;
            incidence = spa.incidence;
        }
        else {
            zenith    = -9999;
            azimuth   = -9999;
            incidence = -9999;
        }
        
        output_arr_ptr[0][j] = zenith; output_arr_ptr[1][j] = azimuth;
        if (spa.function == SPA_ZA_INC) output_arr_ptr[2][j] = incidence;
        
    }
    
    /* Clean up. */
    for (i=0; i<N_ARRAY_IN; i++) Py_XDECREF(input_arr[input_arr_map[i]]);
    
    /* Build the output tuple */
    if (spa.function == SPA_ZA)
        return Py_BuildValue("OO", output_arr[0], output_arr[1]);
    else
        return Py_BuildValue("OOO", output_arr[0], output_arr[1], output_arr[2]);
}
