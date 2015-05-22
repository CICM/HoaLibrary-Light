/*
// Copyright (c) 2012-2015 Eliott Paris, Julien Colafrancesco, Thomas Le Meur & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_SIGNAL_LIGHT
#define DEF_HOA_SIGNAL_LIGHT

#include "Defs.hpp"

namespace hoa
{
    //! The signal class perform all the signal operations with matrix and vector.
    /** The signal class has to be used for every hoa signal operation.
     */
    template <typename T> class Signal
    {
    public:
        //! Multiply a matrix by a vector.
        /** Multiply a matrix by a vector.
        @param inputSize   The number of samples of the input.
        @param outputSize  The number of samples of the input.
        @param vector      The vector.
        @param matrix      The matrix.
        @param outputs     The outputs.
         */
        virtual inline void matrix_vector_mul(const ulong inputsize, const ulong outputsize, const T* vector, const T* matrix, T* outputs) = 0;

        //! Multiply a matrix by a matrix.
        /** Multiply a matrix by a matrix.
        @param nbOfInRow   The matrix number of in row.
        @param nbOfOutCol  The matrix number of out column.
        @param nbOfInCol   The matrix number of in column.
        @param in1         The first matrix.
        @param in2         The second matrix.
        @param out         The final matrix.
         */
        virtual inline void matrix_matrix_mul(const ulong innrow, const ulong outcolumn, const ulong incolumn, const T* in1, const T* in2, T* out) = 0;

        //! Get the max value of a vector.
        /** Get the max value of a vector.
        @param   vectorSize   The size of the vector.
        @param   vector       The vector.
        @return  The max value of the vector
         */
        virtual inline T vector_max(const ulong vectorsize, const T* vector) = 0;

        //! Get the sum value of a vector.
        /** Get the sum value of a vector.
        @param   vectorSize   The size of the vector.
        @param   vector       The vector.
        @return  The sum value of the vector
         */
        virtual inline T vector_sum(const ulong vectorsize, const T* vector) = 0;

        //! Multiply each element of a vector by a factor.
        /** Multiply each element of a vector by a factor.
        @param   vectorSize   The size of the vector.
        @param   factor       The factor of the scale.
        @param   vector       The vector.
         */
        virtual inline void vector_scale(const ulong vectorsize, const T value, T* vector) = 0;

        //! Clear a vector.
        /** Clear a vector.
        @param   vectorSize   The size of the vector.
        @param   vector       The vector.
         */
        virtual inline void vector_clear(const ulong vectorsize, T* vector) = 0;

        //! Copy a vector into an other.
        /** Copy a vector into an other.
        @param   vectorSize   The size of the vector.
        @param   srcVector    The source vector.
        @param   destVector   The destination vector.
         */
        virtual inline void vector_copy(const ulong vectorsize, const T* source, T* dest) = 0;

        //! Add a vector to an other.
        /** Add a vector to an other value by value.
        @param   vectorSize   The size of the vector.
        @param   vectortoAdd  The vector to add.
        @param   destVector   The destination vector.
         */
        virtual inline void vector_add(const ulong vectorsize, const T* source, T* dest) = 0;

        //! Add a vector to an other.
        /** Add a vector to an other with a fixed step.
        @param   vectorSize      The size of the vector.
        @param   vectorToAdd     The vector to add.
        @param   vectorToAddStep The step of the vector to add to select the value to add.
        @param   destVector      The destination vector.
        @param   destVectorStep  The step of the destination vector to select the value to add.
         */
        virtual inline void vector_add(const ulong vectorsize, const T* source, const ulong incs, T* dest, const ulong incd) = 0;

        //! Compute the dot product of two vectors.
        /** Compute the dot product of two vectors.
        @param   vectorSize  The size of the vector.
        @param   vector1     The first vector.
        @param   vector2     The second vector.
        @return The dot product of the two vectors.
         */
        virtual inline T vectors_dot_product(const ulong vectorsize, const T* vector1, const T* vector2) = 0;
    };

    template<> class Signal<float>
    {
    public:
        //! Multiply a matrix by a vector.
        /** Multiply a matrix by a vector.
        @param inputSize   The number of samples of the input.
        @param outputSize  The number of samples of the input.
        @param vector      The vector.
        @param matrix      The matrix.
        @param outputs     The outputs.
         */
        static inline void matrix_vector_mul(const ulong inputsize, const ulong outputsize, const float* vector, const float* matrix, float* outputs)
        {
            cblas_sgemv(CblasRowMajor, CblasNoTrans, (const int)outputsize, (const int)inputsize, 1.f, matrix, (const int)inputsize, vector, 1, 0.f, outputs, 1);
        }

        //! Multiply a matrix by a matrix.
        /** Multiply a matrix by a matrix.
        @param nbOfInRow   The matrix number of in row.
        @param nbOfOutCol  The matrix number of out column.
        @param nbOfInCol   The matrix number of in column.
        @param in1         The first matrix.
        @param in2         The second matrix.
        @param out         The final matrix.
         */
        static inline void matrix_matrix_mul(const ulong innrow, const ulong outcolumn, const ulong incolumn, const float* in1, const float* in2, float* out)
        {
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, (const int)innrow, (const int)outcolumn, (const int)incolumn, 1., in1, incolumn, in2, outcolumn, 0., out,  outcolumn);
        }

        //! Get the max value of a vector.
        /** Get the max value of a vector.
        @param   vectorSize   The size of the vector.
        @param   vector       The vector.
        @return  The max value of the vector
         */
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

        //! Get the sum value of a vector.
        /** Get the sum value of a vector.
        @param   vectorSize   The size of the vector.
        @param   vector       The vector.
        @return  The sum value of the vector
         */
        static inline float vector_sum(const ulong vectorsize, const float* vector)
        {
            return cblas_sasum(vectorsize, vector, 1);
        }

        //! Multiply each element of a vector by a factor.
        /** Multiply each element of a vector by a factor.
        @param   vectorSize   The size of the vector.
        @param   factor       The factor of the scale.
        @param   vector       The vector.
         */
        static inline void vector_scale(const ulong vectorsize, const float value, float* vector)
        {
            cblas_sscal((const int)vectorsize, value, vector, 1.);
        }

        //! Clear a vector.
        /** Clear a vector.
        @param   vectorSize   The size of the vector.
        @param   vector       The vector.
         */
        static inline void vector_clear(const ulong vectorsize, float* vector)
        {
#ifdef __APPLE__
            vDSP_vclr(vector, 1, vectorsize);
#else
            memset(vector, 0, vectorsize * sizeof(float));
#endif
        }

        //! Copy a vector into an other.
        /** Copy a vector into an other.
        @param   vectorSize   The size of the vector.
        @param   srcVector    The source vector.
        @param   destVector   The destination vector.
         */
        static inline void vector_copy(const ulong vectorsize, const float* source, float* dest)
        {
            cblas_scopy(vectorsize, source, 1, dest, 1);
        }

        //! Add a vector to an other.
        /** Add a vector to an other value by value.
        @param   vectorSize   The size of the vector.
        @param   vectortoAdd  The vector to add.
        @param   destVector   The destination vector.
         */
        static inline void vector_add(const ulong vectorsize, const float* source, float* dest)
        {
            cblas_saxpy(vectorsize, 1., source, 1, dest, 1);
        }

        //! Add a vector to an other.
        /** Add a vector to an other with a fixed step.
        @param   vectorSize      The size of the vector.
        @param   vectorToAdd     The vector to add.
        @param   vectorToAddStep The step of the vector to add to select the value to add.
        @param   destVector      The destination vector.
        @param   destVectorStep  The step of the destination vector to select the value to add.
         */
        static inline void vector_add(const ulong vectorsize, const float* source, const ulong incs, float* dest, const ulong incd)
        {
            cblas_saxpy(vectorsize, 1., source, incs, dest, incd);
        }

        //! Compute the dot product of two vectors.
        /** Compute the dot product of two vectors.
        @param   vectorSize  The size of the vector.
        @param   vector1     The first vector.
        @param   vector2     The second vector.
        @return The dot product of the two vectors.
         */
        static inline float vectors_dot_product(const ulong vectorsize, const float* vector1, const float* vector2)
        {
            return cblas_sdot(vectorsize, vector1, 1, vector2, 1);
        }
    };

    template<> class Signal<double>
    {
    public:
        //! Multiply a matrix by a vector.
        /** Multiply a matrix by a vector.
        @param inputSize   The number of samples of the input.
        @param outputSize  The number of samples of the input.
        @param vector      The vector.
        @param matrix      The matrix.
        @param outputs     The outputs.
         */
        static inline void matrix_vector_mul(const ulong inputsize, const ulong outputsize, const double* vector, const double* matrix, double* outputs)
        {
            cblas_dgemv(CblasRowMajor, CblasNoTrans, (const int)outputsize, (const int)inputsize, 1., matrix, (const int)inputsize, vector, 1, 0., outputs, 1);
        }

        //! Multiply a matrix by a matrix.
        /** Multiply a matrix by a matrix.
        @param nbOfInRow   The matrix number of in row.
        @param nbOfOutCol  The matrix number of out column.
        @param nbOfInCol   The matrix number of in column.
        @param in1         The first matrix.
        @param in2         The second matrix.
        @param out         The final matrix.
         */
        static inline void matrix_matrix_mul(const ulong innrow, const ulong outcolumn, const ulong incolumn, const double* in1, const double* in2, double* out)
        {
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, (const int)innrow, (const int)outcolumn, (const int)incolumn, 1., in1, incolumn, in2, outcolumn, 0., out,  outcolumn);
        }

        //! Get the max value of a vector.
        /** Get the max value of a vector.
        @param   vectorSize   The size of the vector.
        @param   vector       The vector.
        @return  The max value of the vector
         */
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

        //! Get the sum value of a vector.
        /** Get the sum value of a vector.
        @param   vectorSize   The size of the vector.
        @param   vector       The vector.
        @return  The sum value of the vector
         */
        static inline double vector_sum(const ulong vectorsize, const double* vector)
        {
            return cblas_dasum(vectorsize, vector, 1);
        }

        //! Multiply each element of a vector by a factor.
        /** Multiply each element of a vector by a factor.
        @param   vectorSize   The size of the vector.
        @param   factor       The factor of the scale.
        @param   vector       The vector.
         */
        static inline void vector_scale(const ulong vectorsize, const double value, double* vector)
        {
            cblas_dscal((const int)vectorsize, value, vector, 1.);
        }

        //! Clear a vector.
        /** Clear a vector.
        @param   vectorSize   The size of the vector.
        @param   vector       The vector.
         */
        static inline void vector_clear(const ulong vectorsize, double* vector)
        {
#ifdef __APPLE__
            vDSP_vclrD(vector, 1, vectorsize);
#else
            memset(vector, 0, vectorsize * sizeof(float));
#endif
        }

        //! Copy a vector into an other.
        /** Copy a vector into an other.
        @param   vectorSize   The size of the vector.
        @param   srcVector    The source vector.
        @param   destVector   The destination vector.
         */
        static inline void vector_copy(const ulong vectorsize, const double* source, double* dest)
        {
            cblas_dcopy(vectorsize, source, 1, dest, 1);
        }

        //! Add a vector to an other.
        /** Add a vector to an other value by value.
        @param   vectorSize   The size of the vector.
        @param   vectortoAdd  The vector to add.
        @param   destVector   The destination vector.
         */
        static inline void vector_add(const ulong vectorsize, const double* source, double* dest)
        {
            cblas_daxpy(vectorsize, 1., source, 1, dest, 1);
        }

        //! Add a vector to an other.
        /** Add a vector to an other with a fixed step.
        @param   vectorSize      The size of the vector.
        @param   vectorToAdd     The vector to add.
        @param   vectorToAddStep The step of the vector to add to select the value to add.
        @param   destVector      The destination vector.
        @param   destVectorStep  The step of the destination vector to select the value to add.
         */
        static inline void vector_add(const ulong vectorsize, const double* source, const ulong incs, double* dest, const ulong incd)
        {
            cblas_daxpy(vectorsize, 1., source, incs, dest, incd);
        }

        //! Compute the dot product of two vectors.
        /** Compute the dot product of two vectors.
        @param   vectorSize  The size of the vector.
        @param   vector1     The first vector.
        @param   vector2     The second vector.
        @return The dot product of the two vectors.
         */
        static inline double vectors_dot_product(const ulong vectorsize, const double* vector1, const double* vector2)
        {
            return cblas_ddot(vectorsize, vector1, 1, vector2, 1);
        }
    };
}

#endif
