/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_3D_VECTOR
#define DEF_HOA_3D_VECTOR

#include "Ambisonic_2D.hpp"

namespace hoa
{
    //! The ambisonic vector.
    /** The vector class compute the energy and the velocity vector of a soudfield for a set of channels. It is an useful tool to characterize the quality of the sound field resitution. For futher information : Michael A. Gerzon, General metatheorie of auditory localisation. Audio Engineering Society Preprint, 3306, 1992. This class retreive the cartesian coordinates of the vectors, the abscissa and the ordinate.
     */
    template <typename T> class Vector : public Planewaves<T>
    {
    private:
        T* m_channels_square;
    public:
        
        //! The vector constructor.
        /**	The optimization constructor allocates and initialize the member values to computes vectors. The number of channels must be at least 1.
            @param     numberOfChannels	The number of channels.
         */
        Vector(unsigned long numberOfChannels) noexcept : Planewaves<T>(numberOfChannels)
        {
            m_channels_square = new T[Planewaves<T>::m_number_of_channels];
        }
        
        //! The optimization destructor.
        /**	The optimization destructor free the memory.
         */
        ~Vector()
        {
            delete [] m_channels_square;
        }
        
        //! This method compute the energy and velocity vectors.
        /**	You should use this method for in-place or not-in-place processing and compute the vectors sample by sample. The inputs array contains the channels samples and the minimum size must be the number of channels. The outputs array contains the vectors cartesian coordinates and the minimum size must be 4. The coordinates arrangement in the outputs array is velocity abscissa, velocity ordinate, energy abscissa, energy ordinate.
         
            @param     inputs   The inputs array.
            @param     outputs  The outputs array.
         */
        inline void process(const T* inputs, T* outputs) noexcept
        {
            processVelocity(inputs, outputs);
            processEnergy(inputs, outputs+2);
        }
        
        //! This method compute the velocity vector.
        /**	You should use this method for in-place or not-in-place processing and compute the vectors sample by sample. The inputs array and contains the channels samples and the minimum size must be the number of channels. The outputs array contains the vectors cartesian coordinates and the minimum size must be 2. The coordinates arrangement in the outputs array is velocity abscissa, velocity ordinate.
         
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        inline void processVelocity(const T* inputs, T* outputs) noexcept
        {
            T veclocitySum = (*inputs++);
            for(unsigned long i = 1; i < Planewaves<T>::m_number_of_channels; i++)
                veclocitySum += (*inputs++);
            
            const T velocityAbscissa = vectors_dot_product(Planewaves<T>::m_number_of_channels, inputs, Planewaves<T>::m_channels_abscissa);
            const T velocityOrdinate = vectors_dot_product(Planewaves<T>::m_number_of_channels, inputs, Planewaves<T>::m_channels_ordinate);
            if(veclocitySum)
            {
                (*outputs++) = velocityAbscissa / veclocitySum;
                (*outputs)  = velocityOrdinate / veclocitySum;
            }
            else
            {
                (*outputs++) = 0.;
                (*outputs) = 0.;
            }
        }
        
        //! This method compute the energy vector.
        /**	You should use this method for in-place or not-in-place processing and compute the vectors sample by sample. The inputs array and contains the channels samples and the minimum size must be the number of harmonics. The outputs array contains the vectors cartesian coordinates and the minimum size must be 2. The coordinates arrangement in the outputs array is energy abscissa, energy ordinate.
         
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        inline void processEnergy(const T* inputs, T* outputs) noexcept
        {
            (*m_channels_square) = (*inputs) * (*inputs);
            for(unsigned long i = 1; i < Planewaves<T>::m_number_of_channels; i++)
                m_channels_square[i] = inputs[i] * inputs[i];
            
            T energySum = vector_sum(Planewaves<T>::m_number_of_channels, m_channels_square);
            const T energyAbscissa = vectors_dot_product(Planewaves<T>::m_number_of_channels, m_channels_square, Planewaves<T>::m_channels_abscissa);
            const T energyOrdinate = vectors_dot_product(Planewaves<T>::m_number_of_channels, m_channels_square, Planewaves<T>::m_channels_ordinate);
            
            if(energySum)
            {
                (*outputs++) = energyAbscissa / energySum;
                (*outputs) = energyOrdinate / energySum;
            }
            else
            {
                (*outputs++) = 0.;
                (*outputs) = 0.;
            }
        }
    };
}

#endif



