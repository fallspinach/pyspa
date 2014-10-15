#include <Python.h>
#include <numpy/arrayobject.h>
#include "spa.h"

/* Docstrings */
static char module_docstring[] =
    "This module provides an interface for the Solar Position Algorithm (SPA).";
static char     za_docstring[] =
    "Calculate solar zenith and azimuth angles.";
static char za_inc_docstring[] =
    "Calculate solar zenith, azimuth, and incidence angles.";

/* Available functions */
static PyObject     *spa_za(PyObject *self, PyObject *args);
static PyObject *spa_za_inc(PyObject *self, PyObject *args);

/* Module specification */
static PyMethodDef module_methods[] = {
    {"za", spa_za, METH_VARARGS, za_docstring},
    {"za_inc", spa_za_inc, METH_VARARGS, za_inc_docstring},
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

static PyObject *spa_za(PyObject *self, PyObject *args)
{
    static const int N_DOUBLE_IN  = 0;
    static const int N_ARRAY_IN   = 8;
    static const int N_DOUBLE_OUT = 0;
    static const int N_ARRAY_OUT  = 2;
    
    PyObject *input_obj[128];
    double    input_double[64];
    
    PyArrayObject *input_arr[64],     *output_arr[64];
    double        *input_arr_ptr[64], *output_arr_ptr[64];
    
    double    year, month, day, hour, minute, second, latitude, longitude;
    double    zenith, azimuth;
    
    int flag;
    
    int      i, j, ndim, dims[NPY_MAXDIMS], rc=0;
    spa_data spa;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OOOOOOOO", &input_obj[0], &input_obj[1], &input_obj[2], &input_obj[3], &input_obj[4], 
                                            &input_obj[5], &input_obj[6], &input_obj[7]))
        return NULL;
    
    /* Interpret the input objects as numpy arrays. */
    for (i=0; i<N_ARRAY_IN; i++)
        input_arr[i] = (PyArrayObject *) PyArray_FROM_OTF(input_obj[i+N_DOUBLE_IN], NPY_DOUBLE, NPY_IN_ARRAY);
    
    /* If that didn't work, throw an exception. */
    flag = 0; for (i=0; i<N_ARRAY_IN; i++) if (input_arr[i]==NULL) flag = 1;
    if (flag==1) {
        for (i=0; i<N_ARRAY_IN; i++) Py_XDECREF(input_arr[i]);
        return NULL;
    }
    
    /* Get pointers to the data as C-types. */
    for (i=0; i<N_ARRAY_IN; i++) input_arr_ptr[i] = (double*) PyArray_DATA(input_arr[i]);
    
    /* How many data points are there? */
    ndim = (int)PyArray_NDIM(input_arr[0]);
    
    int N=1;
    for (j=0; j<ndim; j++) {
         dims[j] = PyArray_DIM(input_arr[0], j);
         N *= dims[j];
         /* check for dimension size compatibility */
         flag = 0; for (i=1; i<N_ARRAY_IN; i++) if (dims[j] != PyArray_DIM(input_arr[i], j)) flag = 1;
         if (flag==1) {
             for (i=0; i<N_ARRAY_IN; i++) Py_XDECREF(input_arr[i]);
             PyErr_SetString(PyExc_RuntimeError, "different dimensions of input arrays.");
             return NULL;
         }
    }
        
    for (i=0; i<N_ARRAY_OUT; i++) {
        output_arr[i]     = (PyArrayObject *) PyArray_FromDims(ndim, dims, NPY_DOUBLE);
        output_arr_ptr[i] = (double*) PyArray_DATA(output_arr[i]);
    }
    
    /* Call the external C function to compute the chi-squared. */
    for (i=0; i<N; i++) {
        
        year   = input_arr_ptr[0][i];   month = input_arr_ptr[1][i];      day = input_arr_ptr[2][i];      hour = input_arr_ptr[3][i];
        minute = input_arr_ptr[4][i];  second = input_arr_ptr[5][i]; latitude = input_arr_ptr[6][i]; longitude = input_arr_ptr[7][i];
        
        /* some checks */
        if (longitude>180) longitude -= 360;
        
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
        spa.elevation     = 0;
        spa.pressure      = 1000;
        spa.temperature   = 0;
        spa.slope         = 0;
        spa.azm_rotation  = 0;
        spa.atmos_refract = 0.5667;
        spa.function      = SPA_ZA;
        
        rc = spa_calculate(&spa);
        
        if (rc == 0) {
            zenith  = spa.zenith;
            azimuth = spa.azimuth;
        }
        else {
            zenith  = -9999;
            azimuth = -9999;
        }
        
        output_arr_ptr[0][i] = zenith; output_arr_ptr[1][i] = azimuth;
        
    }
    
    /* Clean up. */
    for (i=0; i<N_ARRAY_IN; i++) Py_XDECREF(input_arr[i]);
    
    // PyErr_SetString(PyExc_RuntimeError, "SPA here 1."); return NULL;

    /* Build the output tuple */
    return Py_BuildValue("OO", output_arr[0], output_arr[1]);
}

static PyObject *spa_za_inc(PyObject *self, PyObject *args)
{
    static const int N_DOUBLE_IN  = 0;
    static const int N_ARRAY_IN   = 11;
    static const int N_DOUBLE_OUT = 0;
    static const int N_ARRAY_OUT  = 3;
    
    PyObject *input_obj[128];
    double    input_double[64];
    
    PyArrayObject *input_arr[64],     *output_arr[64];
    double        *input_arr_ptr[64], *output_arr_ptr[64];
    
    double    year, month, day, hour, minute, second, latitude, longitude, elevation, slope, aspect;
    double    zenith, azimuth, incidence;
    
    int flag;
    
    int      i, j, ndim, dims[NPY_MAXDIMS], rc=0;
    spa_data spa;
    
    static const double m_air   = 0.02896; // molecular mass of air, kg mol^-1
    static const double R_const = 8.3143;  // gas constant, N M mol^-1 K^-1
    static const double g_const = 9.807;   // gravity constant, m s^-2
    static const double T_const = 288.15;  // "default" air temperature, K

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OOOOOOOOOOO", &input_obj[0], &input_obj[1], &input_obj[2], &input_obj[3], &input_obj[4], 
                                               &input_obj[5], &input_obj[6], &input_obj[7], &input_obj[8], &input_obj[9], &input_obj[10]))
        return NULL;
    
    /* Interpret the input objects as numpy arrays. */
    for (i=0; i<N_ARRAY_IN; i++)
        input_arr[i] = (PyArrayObject *) PyArray_FROM_OTF(input_obj[i+N_DOUBLE_IN], NPY_DOUBLE, NPY_IN_ARRAY);
    
    /* If that didn't work, throw an exception. */
    flag = 0; for (i=0; i<N_ARRAY_IN; i++) if (input_arr[i]==NULL) flag = 1;
    if (flag==1) {
        for (i=0; i<N_ARRAY_IN; i++) Py_XDECREF(input_arr[i]);
        return NULL;
    }
    
    /* Get pointers to the data as C-types. */
    for (i=0; i<N_ARRAY_IN; i++) input_arr_ptr[i] = (double*) PyArray_DATA(input_arr[i]);
    
    /* How many data points are there? */
    ndim = (int)PyArray_NDIM(input_arr[0]);
    
    int N=1;
    for (j=0; j<ndim; j++) {
         dims[j] = PyArray_DIM(input_arr[0], j);
         N *= dims[j];
         /* check for dimension size compatibility */
         flag = 0; for (i=1; i<N_ARRAY_IN; i++) if (dims[j] != PyArray_DIM(input_arr[i], j)) flag = 1;
         if (flag==1) {
             for (i=0; i<N_ARRAY_IN; i++) Py_XDECREF(input_arr[i]);
             PyErr_SetString(PyExc_RuntimeError, "different dimensions of input arrays.");
             return NULL;
         }
    }
        
    for (i=0; i<N_ARRAY_OUT; i++) {
        output_arr[i]     = (PyArrayObject *) PyArray_FromDims(ndim, dims, NPY_DOUBLE);
        output_arr_ptr[i] = (double*) PyArray_DATA(output_arr[i]);
    }
    
    /* Call the external C function to compute the chi-squared. */
    for (i=0; i<N; i++) {
        
        year      = input_arr_ptr[0][i];   month = input_arr_ptr[1][i];      day = input_arr_ptr[2][i];      hour = input_arr_ptr[3][i];
        minute    = input_arr_ptr[4][i];  second = input_arr_ptr[5][i]; latitude = input_arr_ptr[6][i]; longitude = input_arr_ptr[7][i];
        elevation = input_arr_ptr[8][i];   slope = input_arr_ptr[9][i];   aspect = input_arr_ptr[10][i]; 
        
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
        spa.function      = SPA_ZA_INC;
        
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
        
        output_arr_ptr[0][i] = zenith; output_arr_ptr[1][i] = azimuth; output_arr_ptr[2][i] = incidence;
        
    }
    
    /* Clean up. */
    for (i=0; i<N_ARRAY_IN; i++) Py_XDECREF(input_arr[i]);
    
    // PyErr_SetString(PyExc_RuntimeError, "SPA here 1."); return NULL;

    /* Build the output tuple */
    return Py_BuildValue("OOO", output_arr[0], output_arr[1], output_arr[2]);
}
