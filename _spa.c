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
    PyObject  *year_obj, *month_obj, *day_obj, *hour_obj, *minute_obj, *second_obj, *latitude_obj, *longitude_obj;
    
    int      i, ndim, dims[NPY_MAXDIMS], rc=0;
    spa_data spa;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OOOOOOOO", &year_obj, &month_obj, &day_obj, &hour_obj, &minute_obj, &second_obj,
                                            &latitude_obj, &longitude_obj))
        return NULL;

    
    /* Interpret the input objects as numpy arrays. */
    PyObject      *year_arr = PyArray_FROM_OTF(     year_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject     *month_arr = PyArray_FROM_OTF(    month_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject       *day_arr = PyArray_FROM_OTF(      day_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject      *hour_arr = PyArray_FROM_OTF(     hour_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject    *minute_arr = PyArray_FROM_OTF(   minute_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject    *second_arr = PyArray_FROM_OTF(   second_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject  *latitude_arr = PyArray_FROM_OTF( latitude_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *longitude_arr = PyArray_FROM_OTF(longitude_obj, NPY_DOUBLE, NPY_IN_ARRAY);

    /* If that didn't work, throw an exception. */
    if (year_arr == NULL || month_arr == NULL || day_arr == NULL  || hour_arr == NULL
        || minute_arr == NULL|| second_arr == NULL || latitude_arr == NULL || longitude_arr == NULL) {
        Py_XDECREF(     year_arr);
        Py_XDECREF(    month_arr);
        Py_XDECREF(      day_arr);
        Py_XDECREF(     hour_arr);
        Py_XDECREF(   minute_arr);
        Py_XDECREF(   second_arr);
        Py_XDECREF( latitude_arr);
        Py_XDECREF(longitude_arr);
        return NULL;
    }

    /* Get pointers to the data as C-types. */
    double      *year = (double*)PyArray_DATA(     year_arr);
    double     *month = (double*)PyArray_DATA(    month_arr);
    double       *day = (double*)PyArray_DATA(      day_arr);
    double      *hour = (double*)PyArray_DATA(     hour_arr);
    double    *minute = (double*)PyArray_DATA(   minute_arr);
    double    *second = (double*)PyArray_DATA(   second_arr);
    double  *latitude = (double*)PyArray_DATA( latitude_arr);
    double *longitude = (double*)PyArray_DATA(longitude_arr);
    
    /* How many data points are there? */
    ndim =   (int)PyArray_NDIM(year_arr);
    
    int N=1;
    for (i=0; i<ndim; i++) {
         dims[i] = PyArray_DIM(year_arr, i);
         N *= dims[i];
         /* check for dimension size compatibility */
         if ( dims[i] != PyArray_DIM(month_arr, i) || dims[i] != PyArray_DIM(day_arr, i) || dims[i] != PyArray_DIM(hour_arr, i) || dims[i] != PyArray_DIM(minute_arr, i) ||
              dims[i] != PyArray_DIM(second_arr, i) || dims[i] != PyArray_DIM(latitude_arr, i) || dims[i] != PyArray_DIM(longitude_arr, i) ) {
             PyErr_SetString(PyExc_RuntimeError, "different dimensions of input arrays.");
             return NULL;
         }
    }
    
    // fprintf(stderr, "ndim = %d, N=%d, dims[0]=%d, dims[1]=%d\n", ndim, N, dims[0], dims[1]);
    // PyErr_SetString(PyExc_RuntimeError, "SPA");
   
    PyArrayObject  *zenith_arr = (PyArrayObject *) PyArray_FromDims(ndim, dims, NPY_DOUBLE);
    PyArrayObject *azimuth_arr = (PyArrayObject *) PyArray_FromDims(ndim, dims, NPY_DOUBLE);    
    
    double  *zenith = (double*)PyArray_DATA( zenith_arr);
    double *azimuth = (double*)PyArray_DATA(azimuth_arr);
    
    /* Call the external C function to compute the chi-squared. */
    for (i=0; i<N; i++) {
        
        /* some checks */
        if (longitude[i]>180) longitude[i] -= 360;
        
        spa.year          = (int) year[i];
        spa.month         = (int) month[i];
        spa.day           = (int) day[i];
        spa.hour          = (int) hour[i];
        spa.minute        = (int) minute[i];
        spa.second        = (int) second[i];
        spa.latitude      = latitude[i];
        spa.longitude     = longitude[i];
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
            zenith[i]  = spa.zenith;
            azimuth[i] = spa.azimuth;
        }
        else {
            zenith[i]  = -9999;
            azimuth[i] = -9999;
        }
        
    }
    
    /* Clean up. */
    Py_DECREF(year_arr);
    Py_DECREF(month_arr);
    Py_DECREF(day_arr);
    Py_DECREF(hour_arr);
    Py_DECREF(minute_arr);
    Py_DECREF(second_arr);
    Py_DECREF(latitude_arr);
    Py_DECREF(longitude_arr);
    
    // PyErr_SetString(PyExc_RuntimeError, "SPA here 1."); return NULL;

    /* Build the output tuple */
    return Py_BuildValue("OO", zenith_arr, azimuth_arr);
}

static PyObject *spa_za_inc(PyObject *self, PyObject *args)
{
    PyObject  *year_obj, *month_obj, *day_obj, *hour_obj, *minute_obj, *second_obj, *latitude_obj, *longitude_obj, *elevation_obj, *slope_obj, *aspect_obj;
    
    int      i, ndim, dims[NPY_MAXDIMS], rc=0;
    spa_data spa;
    static const double m_air   = 0.02896; // molecular mass of air, kg mol^-1
    static const double R_const = 8.3143;  // gas constant, N M mol^-1 K^-1
    static const double g_const = 9.807;   // gravity constant, m s^-2
    static const double T_const = 288.15;  // "default" air temperature, K

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OOOOOOOOOOO", &year_obj, &month_obj, &day_obj, &hour_obj, &minute_obj, &second_obj,
                                               &latitude_obj, &longitude_obj, &elevation_obj, &slope_obj, &aspect_obj))
        return NULL;

    
    /* Interpret the input objects as numpy arrays. */
    PyObject      *year_arr = PyArray_FROM_OTF(     year_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject     *month_arr = PyArray_FROM_OTF(    month_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject       *day_arr = PyArray_FROM_OTF(      day_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject      *hour_arr = PyArray_FROM_OTF(     hour_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject    *minute_arr = PyArray_FROM_OTF(   minute_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject    *second_arr = PyArray_FROM_OTF(   second_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject  *latitude_arr = PyArray_FROM_OTF( latitude_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *longitude_arr = PyArray_FROM_OTF(longitude_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *elevation_arr = PyArray_FROM_OTF(elevation_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject     *slope_arr = PyArray_FROM_OTF(    slope_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject    *aspect_arr = PyArray_FROM_OTF(   aspect_obj, NPY_DOUBLE, NPY_IN_ARRAY);

    /* If that didn't work, throw an exception. */
    if (year_arr == NULL || month_arr == NULL || day_arr == NULL  || hour_arr == NULL
        || minute_arr == NULL|| second_arr == NULL || latitude_arr == NULL || longitude_arr == NULL
        || elevation_arr == NULL || slope_arr == NULL || aspect_arr == NULL) {
        Py_XDECREF(     year_arr);
        Py_XDECREF(    month_arr);
        Py_XDECREF(      day_arr);
        Py_XDECREF(     hour_arr);
        Py_XDECREF(   minute_arr);
        Py_XDECREF(   second_arr);
        Py_XDECREF( latitude_arr);
        Py_XDECREF(longitude_arr);
        Py_XDECREF(elevation_arr);
        Py_XDECREF(    slope_arr);
        Py_XDECREF(   aspect_arr);
        return NULL;
    }

    /* Get pointers to the data as C-types. */
    double      *year = (double*)PyArray_DATA(     year_arr);
    double     *month = (double*)PyArray_DATA(    month_arr);
    double       *day = (double*)PyArray_DATA(      day_arr);
    double      *hour = (double*)PyArray_DATA(     hour_arr);
    double    *minute = (double*)PyArray_DATA(   minute_arr);
    double    *second = (double*)PyArray_DATA(   second_arr);
    double  *latitude = (double*)PyArray_DATA( latitude_arr);
    double *longitude = (double*)PyArray_DATA(longitude_arr);
    double *elevation = (double*)PyArray_DATA(elevation_arr);
    double     *slope = (double*)PyArray_DATA(    slope_arr);
    double    *aspect = (double*)PyArray_DATA(   aspect_arr);
    
    /* How many data points are there? */
    ndim =   (int)PyArray_NDIM(year_arr);
    
    int N=1;
    for (i=0; i<ndim; i++) {
         dims[i] = PyArray_DIM(year_arr, i);
         N *= dims[i];
         /* check for dimension size compatibility */
         if ( dims[i] != PyArray_DIM(month_arr, i) || dims[i] != PyArray_DIM(day_arr, i) || dims[i] != PyArray_DIM(hour_arr, i) || dims[i] != PyArray_DIM(minute_arr, i) ||
              dims[i] != PyArray_DIM(second_arr, i) || dims[i] != PyArray_DIM(latitude_arr, i) || dims[i] != PyArray_DIM(longitude_arr, i) ||
              dims[i] != PyArray_DIM(elevation_arr, i) || dims[i] != PyArray_DIM(slope_arr, i) || dims[i] != PyArray_DIM(aspect_arr, i) ) {
             PyErr_SetString(PyExc_RuntimeError, "different dimensions of input arrays.");
             return NULL;
         }
    }
    
    // fprintf(stderr, "ndim = %d, N=%d, dims[0]=%d, dims[1]=%d\n", ndim, N, dims[0], dims[1]);
    // PyErr_SetString(PyExc_RuntimeError, "SPA");
   
    PyArrayObject    *zenith_arr = (PyArrayObject *) PyArray_FromDims(ndim, dims, NPY_DOUBLE);
    PyArrayObject   *azimuth_arr = (PyArrayObject *) PyArray_FromDims(ndim, dims, NPY_DOUBLE);    
    PyArrayObject *incidence_arr = (PyArrayObject *) PyArray_FromDims(ndim, dims, NPY_DOUBLE);    
    
    double  *zenith   = (double*)PyArray_DATA(   zenith_arr);
    double *azimuth   = (double*)PyArray_DATA(  azimuth_arr);
    double *incidence = (double*)PyArray_DATA(incidence_arr);
    
    /* Call the external C function to compute the chi-squared. */
    for (i=0; i<N; i++) {
        
        /* some checks */
        if (longitude[i]>180) longitude[i] -= 360;
        
        /* SPA's azm_rotation angle is measured from south and most aspect angle is measured from north */
        aspect[i] -= 180;
        if (aspect[i]<-360) aspect[i] += 360;
        
        spa.year          = (int) year[i];
        spa.month         = (int) month[i];
        spa.day           = (int) day[i];
        spa.hour          = (int) hour[i];
        spa.minute        = (int) minute[i];
        spa.second        = (int) second[i];
        spa.latitude      = latitude[i];
        spa.longitude     = longitude[i];
        spa.timezone      = 0.0;
        spa.delta_t       = 0;
        spa.elevation     = elevation[i];
        spa.pressure      = 1000*exp(-m_air*g_const*elevation[i]/(R_const*T_const));
        spa.temperature   = 0;
        spa.slope         = slope[i];
        spa.azm_rotation  = aspect[i];
        spa.atmos_refract = 0.5667;
        spa.function      = SPA_ZA_INC;
        
        rc = spa_calculate(&spa);
        
        if (rc == 0) {
            zenith[i]    = spa.zenith;
            azimuth[i]   = spa.azimuth;
            incidence[i] = spa.incidence;
        }
        else {
            zenith[i]    = -9999;
            azimuth[i]   = -9999;
            incidence[i] = -9999;
        }
        
    }
    
    /* Clean up. */
    Py_XDECREF(     year_arr);
    Py_XDECREF(    month_arr);
    Py_XDECREF(      day_arr);
    Py_XDECREF(     hour_arr);
    Py_XDECREF(   minute_arr);
    Py_XDECREF(   second_arr);
    Py_XDECREF( latitude_arr);
    Py_XDECREF(longitude_arr);
    Py_XDECREF(elevation_arr);
    Py_XDECREF(    slope_arr);
    Py_XDECREF(   aspect_arr);
    
    // PyErr_SetString(PyExc_RuntimeError, "SPA here 1."); return NULL;

    /* Build the output tuple */
    return Py_BuildValue("OOO", zenith_arr, azimuth_arr, incidence_arr);
}
