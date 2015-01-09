/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_SIGNAL_LIGHT
#define DEF_HOA_SIGNAL_LIGHT

#include "Defs.hpp"

namespace hoa
{
    //! The signal class.
    /** The signal class owns a set of static methods to perform operations over vectors and matrices.
     */
    template <typename T> class Signal;
    
    template<> class Signal<float>
    {
    public:
        
        static inline void matrix_vector_mul(const ulong inputsize, const ulong outputsize, const float* vector, const float* matrix, float* outputs)
        {
            cblas_sgemv(CblasRowMajor, CblasNoTrans, (const int)outputsize, (const int)inputsize, 1.f, matrix, (const int)inputsize, vector, 1, 0.f, outputs, 1);
        }
        
        static inline float vector_max(const ulong vectorsize, const float* vector)
        {
#ifdef __APPLE__
            float result;
            vDSP_maxmgv(vector, 1, &result, vectorsize);
            return result;
#else
            return vector[cblas_isamax((const int)vectorsize, vector, 1)];
#endif
        }
        
        static inline float vector_sum(const ulong vectorsize, const float* vector)
        {
            return cblas_sasum(vectorsize, vector, 1);
        }
        
        static inline void vector_scale(const ulong vectorsize, const float value, float* vector)
        {
            cblas_sscal((const int)vectorsize, value, vector, 1.);
        }
        
        static inline void vector_clear(const ulong vectorsize, float* vector)
        {
#ifdef __APPLE__
            vDSP_vclr(vector, 1, vectorsize);
#else
            memset(vector, 0, vectorsize * sizeof(float));
#endif
        }
        
        static inline void vector_copy(const ulong vectorsize, const float* source, float* dest)
        {
            cblas_scopy(vectorsize, source, 1, dest, 1);
        }
        
        static inline void vector_add(const ulong vectorsize, const float* source, float* dest)
        {
            cblas_saxpy(vectorsize, 1., source, 1, dest, 1);
        }
        
        static inline float vectors_dot_product(const ulong vectorsize, const float* vector1, const float* vector2)
        {
            return cblas_sdot(vectorsize, vector1, 1, vector2, 1);
        }
    };
    
    template<> class Signal<double>
    {
    public:
        
        static inline void matrix_vector_mul(const ulong inputsize, const ulong outputsize, const double* vector, const double* matrix, double* outputs)
        {
            cblas_dgemv(CblasRowMajor, CblasNoTrans, (const int)outputsize, (const int)inputsize, 1., matrix, (const int)inputsize, vector, 1, 0., outputs, 1);
        }
        
        static inline double vector_max(const ulong vectorsize, const double* vector)
        {
#ifdef __APPLE__
            double result;
            vDSP_maxmgvD(vector, 1, &result, vectorsize);
            return result;
#else
            return vector[cblas_idamax((const int)vectorsize, vector, 1)];
#endif
        }
        
        static inline double vector_sum(const ulong vectorsize, const double* vector)
        {
            return cblas_dasum(vectorsize, vector, 1);
        }
        
        static inline void vector_scale(const ulong vectorsize, const double value, double* vector)
        {
            cblas_dscal((const int)vectorsize, value, vector, 1.);
        }
        
        static inline void vector_clear(const ulong vectorsize, double* vector)
        {
#ifdef __APPLE__
            vDSP_vclrD(vector, 1, vectorsize);
#else
            memset(vector, 0, vectorsize * sizeof(float));
#endif
        }
        
        static inline void vector_copy(const ulong vectorsize, const double* source, double* dest)
        {
            cblas_dcopy(vectorsize, source, 1, dest, 1);
        }
        
        static inline void vector_add(const ulong vectorsize, const double* source, double* dest)
        {
            cblas_daxpy(vectorsize, 1., source, 1, dest, 1);
        }
        
        static inline double vectors_dot_product(const ulong vectorsize, const double* vector1, const double* vector2)
        {
            return cblas_ddot(vectorsize, vector1, 1, vector2, 1);
        }
    };
}

#endif
