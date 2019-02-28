/*
 // Copyright (c) 2012-2017 CICM - Universite Paris 8 - Labex Arts H2H.
 // Authors :
 // 2012: Pierre Guillot, Eliott Paris & Julien Colafrancesco.
 // 2012-2015: Pierre Guillot & Eliott Paris.
 // 2015: Pierre Guillot & Eliott Paris & Thomas Le Meur (Light version)
 // 2016-2017: Pierre Guillot.
 // For information on usage and redistribution, and for a DISCLAIMER OF ALL
 // WARRANTIES, see the file, "LICENSE.txt," in this distribution.
 */

#pragma once

#include "Hoa_Defs.hpp"

namespace hoa
{
    // ================================================================================ //
    // SIGNAL //
    // ================================================================================ //
    
    //! @brief The signal class perform all the signal operations with matrix and vector.
    //! @details The signal class has to be used for every hoa signal operation.
    template<typename T>
    class Signal
    {
    public:
        
        //! @brief Allocates a vector.
        //! @param size The size of the vector.
        //! @return A pointer to a vector.
        static inline T* alloc(const size_t size) noexcept
        {
            #ifdef __APPLE__
            T* vec = static_cast<T*>(malloc(size * sizeof(T)));
            if(vec) {clear(size, vec);}
            return vec;
            
            #elif _WINDOWS
            T* vec = (T *)_aligned_malloc(size * sizeof(T), (size_t)pow(float(2), int(sizeof(T))));
            if(vec) {clear(size, vec);}
            return vec;
            
            #else
            T* vec = (T *)memalign((size_t)pow(2, sizeof(T)), size * sizeof(T));
            if(vec) {clear(size, vec);}
            return vec;
            #endif
        }
        
        //! @brief Frees a vector.
        //! @param vec A pointer to a vector.
        //! @return A pointer to a vector.
        static inline T* free(T* vec) noexcept
        {
            #ifdef __APPLE__
            if(vec) {std::free(vec);}
            return nullptr;
            
            #elif _WINDOWS
            if(vec) {_aligned_free(vec);}
            return nullptr;
            
            #else
            if(vec) {std::free(vec);}
            return nullptr;
            #endif
        }
        
        //! @brief Multiplies a matrix by a vector.
        //! @param colsize  The size of the input vector and the number of columns.
        //! @param rowsize  The size of the output vector and the number of rows.
        //! @param in      The input vector.
        //! @param in2      The input matrix.
        //! @param output      The output vector.
        static inline void mul(const size_t colsize,
                               const size_t rowsize,
                               const T* in, const T* in2, T* output) noexcept
        {
            for(size_t i = 0ul; i < rowsize; i++)
            {
                T result = 0;
                const T* in1 = in;
                for(size_t j = colsize>>3; j; --j, in1 += 8, in2 += 8)
                {
                    result += in1[0] * in2[0]; result += in1[1] * in2[1];
                    result += in1[2] * in2[2]; result += in1[3] * in2[3];
                    result += in1[4] * in2[4]; result += in1[5] * in2[5];
                    result += in1[6] * in2[6]; result += in1[7] * in2[7];
                }
                for(size_t j = colsize&7; j; --j, in1++, in2++)
                {
                    result += in1[0] * in2[0];
                }
                output[i] = result;
            }
        }
        
        //! @brief Multiplies a matrix by a matrix.
        //! @param m        The number of rows in the first matrix and the number of columns in the second matrix.
        //! @param n        The number of rows in the second matrix and the number of column in the output matrix.
        //! @param l        The number of columns in the first matrix and the number of rows in the output matrix.
        //! @param in1      The first matrix.
        //! @param in2      The second matrix.
        //! @param output   The output matrix.
        static inline void mul(const size_t m, const size_t n, const size_t l,
                               const T* in1, const T* in2, T* output) noexcept
        {
            size_t i, j, k;
            memset(output, 0, m * n * sizeof(T));
            T* out = output;
            for(k = 0; k < l; k++)
            {
                out = output;
                for(i = 0; i < m; i++)
                {
                    const T g0 = in1[l * i + k];
                    if(g0 != 0)
                    {
                        const T* in = in2+n*k;
                        for(j = n; j; j -= 8, out += 8, in += 8)
                        {
                            const T f0 = in[0] * g0, f1 = in[1] * g0, f2 = in[2] * g0, f3 = in[3] * g0;
                            const T f4 = in[4] * g0, f5 = in[5] * g0, f6 = in[6] * g0, f7 = in[7] * g0;
                            out[0] += f0; out[1] += f1; out[2] += f2; out[3] += f3;
                            out[4] += f4; out[5] += f5; out[6] += f6; out[7] += f7;
                        }
                    }
                }
            }
        }
        
        //! @brief Gets the maximum of the absolute values of a vector.
        //! @param   vectorsize   The size of the vector.
        //! @param   vector       The vector.
        //! @return  The maximum of the absolute values of the vector
        static inline T max(const size_t vectorsize, const T* vector) noexcept
        {
            T max = fabs(vector[0]);
            for(size_t i = 1ul; i < vectorsize; i++)
            {
                const T temp = fabs(vector[i]);
                if(temp > max)
                {
                    max = temp;
                }
            }
            return max;
        }
        
        //! @brief Computes the sum of each element of a vector.
        //! @param   size   The size of the vector.
        //! @param   vector The vector.
        //! @return  The sum of each element of the vector
        static inline T sum(const size_t size, const T* vector) noexcept
        {
            T sum = 0;
            for(size_t i = 0ul; i < size; i++)
            {
                sum += fabs(vector[i]);
            }
            return sum;
        }
        
        //! @brief Multiplies each element of a vector by a factor.
        //! @param   size   The size of the vector.
        //! @param   factor The factor of the scale.
        //! @param   vector The vector.
        static inline void scale(const size_t size, const T factor, T* vector) noexcept
        {
            for(size_t i = 0ul; i < size; i++)
            {
                vector[i] *= factor;
            }
        }
        
        //! @brief Clears a vector.
        //! @param   size   The size of the vector.
        //! @param   vector The vector.
        static inline void clear(const size_t size, T* vector) noexcept
        {
            memset(vector, 0, size * sizeof(T));
        }
        
        //! @brief Copies a vector into an other.
        //! @param   size   The size of the vectors.
        //! @param   source The source vector.
        //! @param   dest   The destination vector.
        static inline void copy(const size_t size, const T* source, T* dest) noexcept
        {
            memcpy(dest, source, size * sizeof(T));
        }
        
        //! @brief Copies a vector into an other.
        //! @param   size   The size of the vectors.
        //! @param   source The source vector.
        //! @param   incs   The increment of the source vector.
        //! @param   dest   The destination vector.
        //! @param   incd   The increment of the destination vector.
        static inline void copy(const size_t size, const T* source, const size_t incs, T* dest, const size_t incd) noexcept
        {
            size_t is = 0ul;
            size_t id = 0ul;
            for(size_t i = 0ul; i < size; i++)
            {
                dest[id] = source[is];
                is += incs;
                id += incd;
            }
        }
        
        //! @brief Adds a vector to an other value by value.
        //! @param   size The size of the vectors.
        //! @param   in The source vector.
        //! @param   out The destination vector.
        static inline void add(const size_t size, const T* in, T* out) noexcept
        {
            for(size_t i = size>>3; i; --i, in += 8, out += 8)
            {
                out[0] += in[0]; out[1] += in[1]; out[2] += in[2]; out[3] += in[3];
                out[4] += in[4]; out[5] += in[5]; out[6] += in[6]; out[7] += in[7];
            }
            for(size_t i = size&7; i; --i, in++, out++)
            {
                out[0] += in[0];
            }
        }
        
        //! @brief Adds a vector to an other.
        //! @details Adds a vector to an other.
        //! @param   size   The size of the vectors.
        //! @param   source The source vector.
        //! @param   incs   The increment of the source vector.
        //! @param   dest   The destination vector.
        //! @param   incd   The increment of the destination vector.
        static inline void add(const size_t size, const T* source, const size_t incs, T* dest, const size_t incd) noexcept
        {
            size_t is = 0ul;
            size_t id = 0ul;
            for(size_t i = 0ul; i < size; i++)
            {
                dest[id] += source[is];
                is += incs;
                id += incd;
            }
        }
        
        //! @brief Computes the dot product of two vectors.
        //! @details Computes the dot product of two vectors.
        //! @param   size   The size of the vectors.
        //! @param   in1    The first vector.
        //! @param   in2    The second vector.
        //! @return The dot product of the two vectors.
        static inline T dot(const size_t size, const T* in1, const T* in2) noexcept
        {
            T result = 0;
            for(size_t i = size>>3; i; --i, in1 += 8, in2 += 8)
            {
                result += in1[0] * in2[0]; result += in1[1] * in2[1];
                result += in1[2] * in2[2]; result += in1[3] * in2[3];
                result += in1[4] * in2[4]; result += in1[5] * in2[5];
                result += in1[6] * in2[6]; result += in1[7] * in2[7];
            }
            
            for(size_t i = size&7; i; --i, in1++, in2++)
            {
                result += in1[0] * in2[0];
            }
            
            return result;
        }
    };
}

