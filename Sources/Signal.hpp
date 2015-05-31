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
    template<typename T> class Signal
    {
    public:
        //! Multiply a matrix by a vector.
        /** Multiply a matrix by a vector.
        @param inputsize   The number of samples of the input.
        @param outputsize  The number of samples of the input.
        @param vector      The vector.
        @param matrix      The matrix.
        @param outputs     The outputs.
         */
        static inline void matrix_vector_mul(const ulong inputsize, const ulong outputsize, const T* vector, const T* matrix, T* outputs)
        {
            for(ulong i = 0ul; i < outputsize; i++)
            {
                T sum = 0.f;
                for(ulong j = 0ul; j < inputsize; j++)
                {
                    sum += vector[j] * *(matrix++);
                }
                outputs[i] = sum;
            }
        }

        //! Multiply a matrix by a matrix.
        /** Multiply a matrix by a matrix.
        @param m          The number of in row in the first matrix and the number of column in the second matrix.
        @param n          The matrix number of out column.
        @param l          The matrix number of in column.
        @param in1        The first matrix.
        @param in2        The second matrix.
        @param out        The final matrix.
         */
        static inline void matrix_matrix_mul(const ulong m, const ulong n, const ulong l, const T* in1, const T* in2, T* out)
        {
            ulong i, j, k;
            for(i = 0; i < m; i++)
            {
                for(j = 0; j < n; j++)
                {
                    out[n * i + j]      = 0.f;
                }
            }
            for(k = 0; k < l; k++)
            {
                for(i = 0; i < m; i++)
                {
                    const T temp = in1[l * i + k];
                    if(temp != 0.f)
                    {
                        for(j = 0; j < n; j++)
                        {
                            out[n * i + j] += temp * in2[n * k + j];
                        }
                    }
                }
            }
        }

        //! Get the max value of a vector.
        /** Get the max value of a vector.
        @param   vectorsize   The size of the vector.
        @param   vector       The vector.
        @return  The max value of the vector
         */
        static inline T vector_max(const ulong vectorsize, const T* vector)
        {
            T max = fabs(vector[0]);
            for(ulong i = 1ul; i < vectorsize; i++)
            {
                const T temp = fabs(vector[1]);
                if(temp > max)
                {
                    max = temp;
                }
            }
            return max;
        }

        //! Get the sum value of a vector.
        /** Get the sum value of a vector.
        @param   vectorsize   The size of the vector.
        @param   vector       The vector.
        @return  The sum value of the vector
         */
        static inline T vector_sum(const ulong vectorsize, const T* vector)
        {
            T sum = 0.f;
            for(ulong i = 0ul; i < vectorsize; i++)
            {
                sum += fabs(vector[i]);
            }
            return sum;
        }

        //! Multiply each element of a vector by a factor.
        /** Multiply each element of a vector by a factor.
        @param   vectorsize  The size of the vector.
        @param   value       The factor of the scale.
        @param   vector      The vector.
         */
        static inline void vector_scale(const ulong vectorsize, const T value, T* vector)
        {
            for(ulong i = 0ul; i < vectorsize; i++)
            {
                vector[i]   *= value;
            }
        }

        //! Clear a vector.
        /** Clear a vector.
        @param   vectorsize   The size of the vector.
        @param   vector       The vector.
         */
        static inline void vector_clear(const ulong vectorsize, T* vector)
        {
            for(ulong i = 0ul; i < vectorsize; i++)
            {
                vector[i]   = 0.f;
            }
        }

        //! Copy a vector into an other.
        /** Copy a vector into an other.
        @param   vectorsize   The size of the vector.
        @param   source       The source vector.
        @param   dest         The destination vector.
         */
        static inline void vector_copy(const ulong vectorsize, const T* source, T* dest)
        {
            for(ulong i = 0ul; i < vectorsize; i++)
            {
                dest[i] = source[i];
            }
        }

        //! Copy a vector into an other.
        /** Copy a vector into an other.
         @param   vectorsize   The size of the vector.
         @param   source       The source vector.
         @param   incs         Number of columns of the source vector.
         @param   dest         Number of columns of the source vector.
         @param   incd         Number of columns of the destination vector.
         */
        static inline void vector_copy(const ulong vectorsize, const T* source, const ulong incs, T* dest, const ulong incd)
        {
            ulong is = 0ul;
            ulong id = 0ul;
            for(ulong i = 0ul; i < vectorsize; i++)
            {
                dest[id] = source[is];
                is += incs;
                id += incd;
            }
        }

        //! Add a vector to an other.
        /** Add a vector to an other value by value.
        @param   vectorsize   The size of the vector.
        @param   source       The vector to add.
        @param   dest         The destination vector.
         */
        static inline void vector_add(const ulong vectorsize, const T* source, T* dest)
        {
            for(ulong i = 0ul; i < vectorsize; i++)
            {
                dest[i]     += source[i];
            }
        }

        //! Add a vector to an other.
        /** Add a vector to an other with a fixed step.
        @param   vectorsize  The size of the vector.
        @param   source      The vector to add.
        @param   incs        Number of columns of the source vector.
        @param   dest        Number of columns of the source vector.
        @param   incd        Number of columns of the destination vector.
         */
        static inline void vector_add(const ulong vectorsize, const T* source, const ulong incs, T* dest, const ulong incd)
        {
            ulong is = 0ul;
            ulong id = 0ul;
            for(ulong i = 0ul; i < vectorsize; i++)
            {
                dest[id] += source[is];
                is += incs;
                id += incd;
            }
        }

        //! Compute the dot product of two vectors.
        /** Compute the dot product of two vectors.
        @param   vectorsize  The size of the vector.
        @param   in1     The first vector.
        @param   in2     The second vector.
        @return The dot product of the two vectors.
         */
        static inline T vectors_dot_product(const ulong vectorsize, const T* in1, const T* in2)
        {
            T sum = 0.f;
            for(ulong i = 0ul; i < vectorsize; i++)
            {
                sum += in1[i]   * in2[i];
            }
            return sum;
        }
    };
}

#endif
