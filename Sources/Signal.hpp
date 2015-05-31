/*
// Copiesright (c) 2012-2015 Eliott Paris, Julien Colafrancesco, Thomas Le Meur & Pierre Guillot, CICM, Universite Paris 8.
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
        
        //! Multiplies a matrix by a vector.
        /** Multiplies a matrix by a vector.
        @param colsize  The size of the input vector and the number of columns.
        @param rowsize  The size of the output vector and the number of rows.
        @param vector      The input vector.
        @param matrix      The input matrix.
        @param output      The output vector.
         */
        static inline void mul(const ulong colsize, const ulong rowsize, const T* input, const T* matrix, T* output) noexcept
        {
            for(ulong i = 0ul; i < rowsize; i++)
            {
                T sum = 0.f;
                for(ulong j = 0ul; j < colsize; j++)
                {
                    sum += input[j] * *(matrix++);
                }
                output[i] = sum;
            }
        }

        //! Multiplies a matrix by a matrix.
        /** Multiplies a matrix by a matrix.
        @param m        The number of rows in the first matrix and the number of columns in the second matrix.
        @param n        The number of rows in the second matrix and the number of column in the output matrix.
        @param l        The number of columns in the first matrix and the number of rows in the output matrix.
        @param in1      The first matrix.
        @param in2      The second matrix.
        @param output   The output matrix.
         */
        static inline void mul(const ulong m, const ulong n, const ulong l, const T* in1, const T* in2, T* output) noexcept
        {
            /*
            {
                ulong i, j, k;
                for(i = 0; i < m; i++)
                {
                    for(j = 0; j < n; j++)
                    {
                        output[n * i + j]      = 0.f;
                    }
                }
                
                for(i = 0; i < m; i++)
                {
                    for(j = 0; j < n; j++)
                    {
                        for(k = 0; k < l; k++)
                        {
                            output[n * i + j] += in1[l * i + k] * in2[n * j + k];
                        }
                    }
                }
            }*/
            /*
            const float* a = in1;
            for(ulong i = 0; i < n; i++)
            {
                const float* b = in2;
                for(ulong j = 0; j < l; j++)
                {
                    __m128 c = _mm_setzero_ps();
                    //const float* a2 = a;
                    for(ulong k = 0; k < m; k += 4)
                    {
                        _mm_add_ps(_mm_load_ps(&a[k]), _mm_load_ps(b));
                        b+=4;
                    }
                    c = _mm_hadd_ps(c, c);
                    c = _mm_hadd_ps(c, c);
                    _mm_store_ss(output, c);
                    output++;
                }
                a += m;
            }
            */
            ulong i, j, k;
            T* out = output;
            for(i = 0; i < m; i++)
            {
                for(j = n; j; j -= 8, out += 8)
                {
                    out[0] = 0;
                    out[1] = 0;
                    out[2] = 0;
                    out[3] = 0;
                    out[4] = 0;
                    out[5] = 0;
                    out[6] = 0;
                    out[7] = 0;
                }
            }
            for(k = 0; k < l; k++)
            {
                out = output;
                for(i = 0; i < m; i++)
                {
                    const T g = in1[l * i + k];
                    if(g != 0)
                    {
                        const T* in = in2+n*k;
                        for(j = n; j; j -= 8, out += 8, in += 8)
                        {
                            const T f0 = in[0] * g, f1 = in[1] * g, f2 = in[2] * g, f3 = in[3] * g;
                            const T f4 = in[4] * g, f5 = in[5] * g, f6 = in[6] * g, f7 = in[7] * g;
                            out[0] += f0; out[1] += f1; out[2] += f2; out[3] += f3;
                            out[4] += f4; out[5] += f5; out[6] += f6; out[7] += f7;
                        }
                    }
                }
            }
        }

        //! Gets the maximum of the absolute values of a vector.
        /** Gets the maximum of the absolute values of a vector.
        @param   vectorsize   The size of the vector.
        @param   vector       The vector.
        @return  The maximum of the absolute values of the vector
         */
        static inline T max(const ulong vectorsize, const T* vector) noexcept
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

        //! Computes the sum of each element of a vector.
        /** Computes the sum of each element of a vector.
        @param   size   The size of the vector.
        @param   vector The vector.
        @return  The sum of each element of the vector
         */
        static inline T sum(const ulong size, const T* vector) noexcept
        {
            T sum = 0;
            for(ulong i = 0ul; i < size; i++)
            {
                sum += fabs(vector[i]);
            }
            return sum;
        }

        //! Multiplies each element of a vector by a factor.
        /** Multiplies each element of a vector by a factor.
        @param   size   The size of the vector.
        @param   factor The factor of the scale.
        @param   vector The vector.
         */
        static inline void scale(const ulong size, const T factor, T* vector) noexcept
        {
            for(ulong i = 0ul; i < size; i++)
            {
                vector[i]   *= factor;
            }
        }

        //! Clears a vector.
        /** Clears a vector.
        @param   size   The size of the vector.
        @param   vector The vector.
         */
        static inline void clear(const ulong size, T* vector) noexcept
        {
            for(ulong i = 0ul; i < size; i++)
            {
                vector[i]   = 0.f;
            }
        }

        //! Copies a vector into an other.
        /** Copies a vector into an other.
        @param   size   The size of the vectors.
        @param   source The source vector.
        @param   dest   The destination vector.
         */
        static inline void copy(const ulong size, const T* source, T* dest) noexcept
        {
            for(ulong i = 0ul; i < size; i++)
            {
                dest[i] = source[i];
            }
        }

        //! Copies a vector into an other.
        /** Copies a vector into an other.
         @param   size   The size of the vectors.
         @param   source The source vector.
         @param   incs   The increment of the source vector.
         @param   dest   The destination vector.
         @param   incd   The increment of the destination vector.
         */
        static inline void copy(const ulong size, const T* source, const ulong incs, T* dest, const ulong incd) noexcept
        {
            ulong is = 0ul;
            ulong id = 0ul;
            for(ulong i = 0ul; i < size; i++)
            {
                dest[id] = source[is];
                is += incs;
                id += incd;
            }
        }

        //! Adds a vector to an other.
        /** Adds a vector to an other value by value.
        @param   size   The size of the vectors.
        @param   source The source vector.
        @param   dest   The destination vector.
         */
        static inline void add(const ulong size, const T* source, T* dest) noexcept
        {
            for(ulong i = 0ul; i < size; i++)
            {
                dest[i] += source[i];
            }
        }

        //! Adds a vector to an other.
        /** Adds a vector to an other.
         @param   size   The size of the vectors.
         @param   source The source vector.
         @param   incs   The increment of the source vector.
         @param   dest   The destination vector.
         @param   incd   The increment of the destination vector.
         */
        static inline void add(const ulong size, const T* source, const ulong incs, T* dest, const ulong incd) noexcept
        {
            ulong is = 0ul;
            ulong id = 0ul;
            for(ulong i = 0ul; i < size; i++)
            {
                dest[id] += source[is];
                is += incs;
                id += incd;
            }
        }

        //! Computes the dot product of two vectors.
        /** Computes the dot product of two vectors.
        @param   size   The size of the vectors.
        @param   in1    The first vector.
        @param   in2    The second vector.
        @return The dot product of the two vectors.
         */
        static inline T dot(const ulong size, const T* in1, const T* in2) noexcept
        {
            T sum = 0;
            for(ulong i = 0ul; i < size; i++)
            {
                sum += in1[i] * in2[i];
            }
            return sum;
        }
    };
}

#endif
