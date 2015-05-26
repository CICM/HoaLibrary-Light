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
        @param innrow     The matrix number of in row.
        @param outcolumn  The matrix number of out column.
        @param incolumn   The matrix number of in column.
        @param in1        The first matrix.
        @param in2        The second matrix.
        @param out        The final matrix.
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
        @param   vectorSize  The size of the vector.
        @param   value       The factor of the scale.
        @param   vector      The vector.
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
        @param   source       The source vector.
        @param   dest         The destination vector.
         */
        virtual inline void vector_copy(const ulong vectorsize, const T* source, T* dest) = 0;

        //! Copy a vector into an other.
        /** Copy a vector into an other.
         @param   vectorSize   The size of the vector.
         @param   source       The source vector.
         @param   incs         Number of columns of the source vector.
         @param   dest         Number of columns of the source vector.
         @param   incd         Number of columns of the destination vector.
         */
        virtual inline void vector_copy(const ulong vectorsize, const double* source, const ulong incs, double* dest, const ulong incd) = 0;

        //! Add a vector to an other.
        /** Add a vector to an other value by value.
        @param   vectorSize   The size of the vector.
        @param   source       The vector to add.
        @param   dest         The destination vector.
         */
        virtual inline void vector_add(const ulong vectorsize, const T* source, T* dest) = 0;

        //! Add a vector to an other.
        /** Add a vector to an other with a fixed step.
        @param   vectorSize  The size of the vector.
        @param   source      The vector to add.
        @param   incs        Number of columns of the source vector.
        @param   dest        Number of columns of the source vector.
        @param   incd        Number of columns of the destination vector.
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
#ifdef HOA_USE_CBLAS
            cblas_sgemv(CblasRowMajor, CblasNoTrans, (const int)outputsize, (const int)inputsize, 1.f, matrix, (const int)inputsize, vector, 1, 0.f, outputs, 1);
#else
            outputs[0] = vectors_dot_product(inputsize, vector, matrix);
            for(ulong i = 1ul; i < outputsize; i++)
            {
                outputs[i] = vectors_dot_product(inputsize, vector, matrix+i*inputsize);
            }
#endif
        }

        //! Multiply a matrix by a matrix.
        /** Multiply a matrix by a matrix.
        @param innrow     The matrix number of in row.
        @param outcolumn  The matrix number of out column.
        @param incolumn   The matrix number of in column.
        @param in1        The first matrix.
        @param in2        The second matrix.
        @param out        The final matrix.
         */
        static inline void matrix_matrix_mul(const ulong innrow, const ulong outcolumn, const ulong incolumn, const float* in1, const float* in2, float* out)
        {
#ifdef HOA_USE_CBLAS
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, (const int)innrow, (const int)outcolumn, (const int)incolumn, 1., in1, incolumn, in2, outcolumn, 0., out,  outcolumn);
#else
            const float* a = in1;
            for(ulong c = 0; c < innrow; c++)
            {
                const float* b = in2;
                for(ulong d = 0; d < outcolumn; d++)
                {
                    *out = 0.f;
                    const float* a2 = a;
                    for(ulong k = 0; k < incolumn; k++)
                    {
                        *out += *(a2++) * *(b + k * outcolumn);
                    }
                    ++out;
                    ++b;
                }
                a += incolumn;
            }
#endif
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
#elif HOA_USE_CBLAS
            return vector[cblas_isamax((const int)vectorsize, vector, 1)];
#else
            float max = fabs(vector[0]);
            for(ulong i = 1ul; i < vectorsize; i++)
            {
                const double temp = fabs(vector[1]);
                if(temp > max)
                {
                    max = temp;
                }
            }
            return max;
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
#ifdef HOA_USE_CBLAS
            return cblas_sasum(vectorsize, vector, 1);
#else
            float sum = fabs(vector[0]);
            for(ulong i = 1ul; i < vectorsize; i++)
            {
                sum += fabs(vector[i]);
            }
            return sum;
#endif
        }

        //! Multiply each element of a vector by a factor.
        /** Multiply each element of a vector by a factor.
        @param   vectorSize  The size of the vector.
        @param   value       The factor of the scale.
        @param   vector      The vector.
         */
        static inline void vector_scale(const ulong vectorsize, const float value, float* vector)
        {
#ifdef HOA_USE_CBLAS
            cblas_sscal((const int)vectorsize, value, vector, 1.);
#else
            for(ulong i = 0ul; i < vectorsize; i++)
            {
                vector[i] += value;
            }
#endif
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
        @param   source       The source vector.
        @param   dest         The destination vector.
         */
        static inline void vector_copy(const ulong vectorsize, const float* source, float* dest)
        {
#ifdef HOA_USE_CBLAS
            cblas_scopy(vectorsize, source, 1, dest, 1);
#else
            memcpy(dest, source, vectorsize * sizeof(float));
#endif
        }

        //! Copy a vector into an other.
        /** Copy a vector into an other.
         @param   vectorSize   The size of the vector.
         @param   source       The source vector.
         @param   incs         Number of columns of the source vector.
         @param   dest         Number of columns of the source vector.
         @param   incd         Number of columns of the destination vector.
         */
        static inline void vector_copy(const ulong vectorsize, const float* source, const ulong incs, float* dest, const ulong incd)
        {
#ifdef HOA_USE_CBLAS
            cblas_scopy(vectorsize, source, incs, dest, incd);
#else
            ulong is = incs;
            ulong id = incd;
            dest[0] = source[0];
            for(ulong i = 1; i < vectorsize; i++)
            {
                dest[id] = source[is];
                is += incs;
                id += incd;
            }
#endif
        }

        //! Add a vector to an other.
        /** Add a vector to an other value by value.
        @param   vectorSize   The size of the vector.
        @param   source       The vector to add.
        @param   dest         The destination vector.
         */
        static inline void vector_add(const ulong vectorsize, const float* source, float* dest)
        {
#ifdef HOA_USE_CBLAS
            cblas_saxpy(vectorsize, 1., source, 1, dest, 1);
#else
            for(ulong i = 0ul; i < vectorsize; i++)
            {
                dest[i] += source[i];
            }
#endif
        }

        //! Add a vector to an other.
        /** Add a vector to an other with a fixed step.
        @param   vectorSize  The size of the vector.
        @param   source      The vector to add.
        @param   incs        Number of columns of the source vector.
        @param   dest        Number of columns of the source vector.
        @param   incd        Number of columns of the destination vector.
         */
        static inline void vector_add(const ulong vectorsize, const float* source, const ulong incs, float* dest, const ulong incd)
        {
#ifdef HOA_USE_CBLAS
            cblas_saxpy(vectorsize, 1., source, incs, dest, incd);
#else
            ulong is = incs;
            ulong id = incd;
            dest[0] += source[0];
            for(ulong i = 1; i < vectorsize; i++)
            {
                dest[id] += source[is];
                is += incs;
                id += incd;
            }
#endif

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
#ifdef HOA_USE_CBLAS
            return cblas_sdot(vectorsize, vector1, 1, vector2, 1);
#else
            float sum = vector1[0] * vector2[0];
            for(ulong i = 1ul; i < vectorsize; i++)
            {
                sum += vector1[i] * vector2[i];
            }
            return sum;
#endif

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
#ifdef HOA_USE_CBLAS
            cblas_dgemv(CblasRowMajor, CblasNoTrans, (const int)outputsize, (const int)inputsize, 1., matrix, (const int)inputsize, vector, 1, 0., outputs, 1);
#else
            outputs[0] = vectors_dot_product(inputsize, vector, matrix);
            for(ulong i = 1ul; i < outputsize; i++)
            {
                outputs[i] = vectors_dot_product(inputsize, vector, matrix+i*inputsize);
            }
#endif
        }

        //! Multiply a matrix by a matrix.
        /** Multiply a matrix by a matrix.
        @param innrow     The matrix number of in row.
        @param outcolumn  The matrix number of out column.
        @param incolumn   The matrix number of in column.
        @param in1        The first matrix.
        @param in2        The second matrix.
        @param out        The final matrix.
         */
        static inline void matrix_matrix_mul(const ulong innrow, const ulong outcolumn, const ulong incolumn, const double* in1, const double* in2, double* out)
        {
#ifdef HOA_USE_CBLAS
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, (const int)innrow, (const int)outcolumn, (const int)incolumn, 1., in1, incolumn, in2, outcolumn, 0., out,  outcolumn);
#else
            const double* a = in1;
            for(ulong c = 0; c < innrow; c++)
            {
                const double* b = in2;
                for(ulong d = 0; d < outcolumn; d++)
                {
                    *out = 0.f;
                    const double* a2 = a;
                    for(ulong k = 0; k < incolumn; k++)
                    {
                        *out += *(a2++) * *(b + k * outcolumn);
                    }
                    ++out;
                    ++b;
                }
                a += incolumn;
            }
#endif
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
#elif HOA_USE_CBLAS
            return vector[cblas_idamax((const int)vectorsize, vector, 1)];
#else
            double max = fabs(vector[0]);
            for(ulong i = 1ul; i < vectorsize; i++)
            {
                const double temp = fabs(vector[1]);
                if(temp > max)
                {
                    max = temp;
                }
            }
            return max;
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
#ifdef HOA_USE_CBLAS
            return cblas_dasum(vectorsize, vector, 1);
#else
            double sum = fabs(vector[0]);
            for(ulong i = 1ul; i < vectorsize; i++)
            {
                sum += fabs(vector[i]);
            }
            return sum;
#endif
        }

        //! Multiply each element of a vector by a factor.
        /** Multiply each element of a vector by a factor.
        @param   vectorSize  The size of the vector.
        @param   value       The factor of the scale.
        @param   vector      The vector.
         */
        static inline void vector_scale(const ulong vectorsize, const double value, double* vector)
        {
#ifdef HOA_USE_CBLAS
            cblas_dscal((const int)vectorsize, value, vector, 1.);
#else
            for(ulong i = 0ul; i < vectorsize; i++)
            {
                vector[i] += value;
            }
#endif
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
            memset(vector, 0, vectorsize * sizeof(double));
#endif
        }

        //! Copy a vector into an other.
        /** Copy a vector into an other.
        @param   vectorSize   The size of the vector.
        @param   source       The source vector.
        @param   dest         The destination vector.
         */
        static inline void vector_copy(const ulong vectorsize, const double* source, double* dest)
        {
#ifdef HOA_USE_CBLAS
            cblas_dcopy(vectorsize, source, 1, dest, 1);
#else
            memcpy(dest, source, vectorsize * sizeof(double));
#endif
        }

        //! Copy a vector into an other.
        /** Copy a vector into an other.
         @param   vectorSize   The size of the vector.
         @param   source       The source vector.
         @param   incs         Number of columns of the source vector.
         @param   dest         Number of columns of the source vector.
         @param   incd         Number of columns of the destination vector.
         */
        static inline void vector_copy(const ulong vectorsize, const double* source, const ulong incs, double* dest, const ulong incd)
        {
#ifdef HOA_USE_CBLAS
            cblas_dcopy(vectorsize, source, incs, dest, incd);
#else
            ulong is = incs;
            ulong id = incd;
            dest[0] = source[0];
            for(ulong i = 1; i < vectorsize; i++)
            {
                dest[id] = source[is];
                is += incs;
                id += incd;
            }
#endif
        }

        //! Add a vector to an other.
        /** Add a vector to an other value by value.
        @param   vectorSize   The size of the vector.
        @param   source       The vector to add.
        @param   dest         The destination vector.
         */
        static inline void vector_add(const ulong vectorsize, const double* source, double* dest)
        {
#ifdef HOA_USE_CBLAS
            cblas_daxpy(vectorsize, 1., source, 1, dest, 1);
#else
            for(ulong i = 0ul; i < vectorsize; i++)
            {
                dest[i] += source[i];
            }
#endif
        }

        //! Add a vector to an other.
        /** Add a vector to an other with a fixed step.
        @param   vectorSize  The size of the vector.
        @param   source      The vector to add.
        @param   incs        Number of columns of the source vector.
        @param   dest        Number of columns of the source vector.
        @param   incd        Number of columns of the destination vector.
         */
        static inline void vector_add(const ulong vectorsize, const double* source, const ulong incs, double* dest, const ulong incd)
        {
#ifdef HOA_USE_CBLAS
            cblas_daxpy(vectorsize, 1., source, incs, dest, incd);
#else
            ulong is = incs;
            ulong id = incd;
            dest[0] += source[0];
            for(ulong i = 1; i < vectorsize; i++)
            {
                dest[id] += source[is];
                is += incs;
                id += incd;
            }
#endif
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
#ifdef HOA_USE_CBLAS
            return cblas_ddot(vectorsize, vector1, 1, vector2, 1);
#else
            double sum = vector1[0] * vector2[0];
            for(ulong i = 1ul; i < vectorsize; i++)
            {
                sum += vector1[i] * vector2[i];
            }
            return sum;
#endif
        }
    };
}

#endif
