/*
// Copyright (c) 2012-2015 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_VECTOR_LIGHT
#define DEF_HOA_VECTOR_LIGHT

#include "Planewaves.hpp"

namespace hoa
{
    //! The vector class computes the energy and the velocity vectors for a set of loudspeakers.
    /** The vector class compute the energy and the velocity vectors of a soud field for a set of channels. It is an useful tool to characterize the quality of the sound field resitution. For futher information : Michael A. Gerzon, General metatheorie of auditory localisation. Audio Engineering Society Preprint, 3306, 1992. This class retreive the cartesian coordinates of the vectors.
     */
    template <Dimension D, typename T> class Vector : public Processor< Planewave<D, T> >
    {
        //! The vector constructor.
        /**	The vector constructor allocates and initialize the member values to computes vectors. The number of channels must be at least 1.
         @param     numberOfChannels	The number of channels.
         */
        Vector(const ulong numberOfChannels) noexcept = 0;
        
        //! The vector destructor.
        /**	The vector destructor free the memory.
         */
        virtual ~Vector() noexcept = 0;
        
        //! This method pre-computes the necessary values to process.
        /**	You should use this method before calling the process methods and after changing the azimuth, the elevation or the offset of the channels.
         */
        virtual void computeRendering() noexcept = 0;
        
        //! This method computes the energy and velocity vectors.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array contains the channels samples and the minimum size must be the number of channels. The outputs array contains the vectors cartesian coordinates and the minimum size must be 4 for 2d and 6 for 3d. The coordinates arrangement in the outputs array is velocity abscissa, velocity ordinate, (velocity height), energy abscissa and energy ordinate (and energy height).
         
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        virtual void process(const T* inputs, T* outputs) noexcept = 0;
        
        //! This method computes the velocity vector.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and contains the channels samples and the minimum size must be the number of channels. The outputs array contains the vectors cartesian coordinates and the minimum size must be 2 for 2d and 3 for 3d. The coordinates arrangement in the outputs array is velocity abscissa and velocity ordinate (and velocity height).
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        virtual void processVelocity(const T* inputs, T* outputs) noexcept = 0;
        
        //! This method computes the energy vector.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and contains the channels samples and the minimum size must be the number of harmonics. The outputs array contains the vectors cartesian coordinates and the minimum size must be 2 for 2d and 3 for 3d. The coordinates arrangement in the outputs array is energy abscissa and energy ordinate (and energy height).
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        virtual void processEnergy(const T* inputs, T* outputs) noexcept = 0;
    };
    
#ifndef DOXYGEN_SHOULD_SKIP_THIS
   
    template <typename T> class Vector<Hoa2d, T> : public Processor< Planewave<Hoa2d, T> >
    {
    private:
        T* m_channels_square;
        T* m_channels_abscissa;
        T* m_channels_ordinate;
    public:
        
        //! The vector constructor.
        /**	The vector constructor allocates and initialize the member values to computes vectors. The number of channels must be at least 1.
         @param     numberOfChannels	The number of channels.
         */
        Vector(const ulong numberOfChannels) noexcept : Processor< Planewave<Hoa2d, T> >(numberOfChannels)
        {
            m_channels_square   = new T[Processor< Planewave<Hoa2d, T> >::getNumberOfPlanewaves()];
            m_channels_abscissa = new T[Processor< Planewave<Hoa2d, T> >::getNumberOfPlanewaves()];
            m_channels_ordinate = new T[Processor< Planewave<Hoa2d, T> >::getNumberOfPlanewaves()];
        }
        
        //! The vector destructor.
        /**	The vector destructor free the memory.
         */
        ~Vector() noexcept
        {
            delete [] m_channels_square;
            delete [] m_channels_abscissa;
            delete [] m_channels_ordinate;
        }
        
        inline void computeRendering() noexcept
        {
            for(ulong i = 0; i < Processor< Planewave<Hoa2d, T> >::getNumberOfPlanewaves(); i++)
            {
                m_channels_abscissa[i] = Processor< Planewave<Hoa2d, T> >::getPlanewaveAbscissa(i);
                m_channels_ordinate[i] = Processor< Planewave<Hoa2d, T> >::getPlanewaveOrdinate(i);
            }
        }
        
        //! This method computes the energy and velocity vectors.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array contains the channels samples and the minimum size must be the number of channels. The outputs array contains the vectors cartesian coordinates and the minimum size must be 4. The coordinates arrangement in the outputs array is velocity abscissa, velocity ordinate, energy abscissa and energy ordinate.
            @param     inputs   The inputs array.
            @param     outputs  The outputs array.
         */
        inline void process(const T* inputs, T* outputs) noexcept
        {
            processVelocity(inputs, outputs);
            processEnergy(inputs, outputs+2);
        }
        
        //! This method computes the velocity vector.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and contains the channels samples and the minimum size must be the number of channels. The outputs array contains the vectors cartesian coordinates and the minimum size must be 2. The coordinates arrangement in the outputs array is velocity abscissa and velocity ordinate.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        inline void processVelocity(const T* inputs, T* outputs) noexcept
        {
            T veclocitySum = (*inputs++);
            for(ulong i = 1; i < Processor< Planewave<Hoa2d, T> >::getNumberOfPlanewaves(); i++)
                veclocitySum += (*inputs++);
            
            const T velocityAbscissa = Signal<T>::vectors_dot_product(Processor< Planewave<Hoa2d, T> >::getNumberOfPlanewaves(), inputs, m_channels_abscissa);
            const T velocityOrdinate = Signal<T>::vectors_dot_product(Processor< Planewave<Hoa2d, T> >::getNumberOfPlanewaves(), inputs, m_channels_ordinate);
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
        
        //! This method computes the energy vector.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and contains the channels samples and the minimum size must be the number of harmonics. The outputs array contains the vectors cartesian coordinates and the minimum size must be 2. The coordinates arrangement in the outputs array is energy abscissa and energy ordinate.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        inline void processEnergy(const T* inputs, T* outputs) noexcept
        {
            (*m_channels_square) = (*inputs) * (*inputs);
            for(ulong i = 1; i < Processor< Planewave<Hoa2d, T> >::getNumberOfPlanewaves(); i++)
                m_channels_square[i] = inputs[i] * inputs[i];
            
            T energySum = Signal<T>::vector_sum(Processor< Planewave<Hoa2d, T> >::getNumberOfPlanewaves(), m_channels_square);
            const T energyAbscissa = Signal<T>::vectors_dot_product(Processor< Planewave<Hoa2d, T> >::getNumberOfPlanewaves(), m_channels_square, m_channels_abscissa);
            const T energyOrdinate = Signal<T>::vectors_dot_product(Processor< Planewave<Hoa2d, T> >::getNumberOfPlanewaves(), m_channels_square, m_channels_ordinate);
            
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
    
    template <typename T> class Vector<Hoa3d, T> : public Processor< Planewave<Hoa3d, T> >
    {
    private:
        T* m_channels_square;
        T* m_channels_abscissa;
        T* m_channels_ordinate;
        T* m_channels_height;
    public:
        
        //! The vector constructor.
        /**	The vector constructor allocates and initialize the member values to computes vectors. The number of channels must be at least 1.
         @param     numberOfChannels	The number of channels.
         */
        Vector(const ulong numberOfChannels) noexcept : Processor< Planewave<Hoa3d, T> >(numberOfChannels)
        {
            m_channels_square   = new T[Processor< Planewave<Hoa3d, T> >::getNumberOfPlanewaves()];
            m_channels_abscissa = new T[Processor< Planewave<Hoa3d, T> >::getNumberOfPlanewaves()];
            m_channels_ordinate = new T[Processor< Planewave<Hoa3d, T> >::getNumberOfPlanewaves()];
            m_channels_height   = new T[Processor< Planewave<Hoa3d, T> >::getNumberOfPlanewaves()];
        }
        
        //! The vector destructor.
        /**	The vector destructor free the memory.
         */
        ~Vector() noexcept
        {
            delete [] m_channels_square;
            delete [] m_channels_abscissa;
            delete [] m_channels_ordinate;
            delete [] m_channels_height;
        }
        
        inline void computeRendering() noexcept
        {
            for(ulong i = 0; i < Processor< Planewave<Hoa3d, T> >::getNumberOfPlanewaves(); i++)
            {
                m_channels_abscissa[i] = Processor< Planewave<Hoa3d, T> >::getPlanewaveAbscissa(i);
                m_channels_ordinate[i] = Processor< Planewave<Hoa3d, T> >::getPlanewaveOrdinate(i);
                m_channels_height[i]   = Processor< Planewave<Hoa3d, T> >::getPlanewaveHeight(i);
            }
        }
        
        //! This method compute the energy and velocity vectors.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array contains the channels samples and the minimum size must be the number of channels. The outputs array contains the vectors cartesian coordinates and the minimum size must be 6. The coordinates arrangement in the outputs array is velocity abscissa, velocity ordinate, velocity height, energy abscissa, energy ordinate and energy height.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        inline void process(const T* inputs, T* outputs) noexcept
        {
            processVelocity(inputs, outputs);
            processEnergy(inputs, outputs+3);
        }
        
        //! This method compute the velocity vector.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and contains the channels samples and the minimum size must be the number of channels. The outputs array contains the vectors cartesian coordinates and the minimum size must be 3. The coordinates arrangement in the outputs array is velocity abscissa, velocity ordinate and velocity height.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        inline void processVelocity(const T* inputs, T* outputs) noexcept
        {
            T veclocitySum = (*inputs++);
            for(ulong i = 1; i < Processor< Planewave<Hoa3d, T> >::getNumberOfPlanewaves(); i++)
                veclocitySum += (*inputs++);
                
                const T velocityAbscissa    = Signal<T>::vectors_dot_product(Processor< Planewave<Hoa3d, T> >::getNumberOfPlanewaves(), inputs, m_channels_abscissa);
                const T velocityOrdinate    = Signal<T>::vectors_dot_product(Processor< Planewave<Hoa3d, T> >::getNumberOfPlanewaves(), inputs, m_channels_ordinate);
                const T velocityHeight      = Signal<T>::vectors_dot_product(Processor< Planewave<Hoa3d, T> >::getNumberOfPlanewaves(), inputs, m_channels_height);
                
                if(veclocitySum)
                {
                    (*outputs++) = velocityAbscissa / veclocitySum;
                    (*outputs++) = velocityOrdinate / veclocitySum;
                    (*outputs)   = velocityHeight / veclocitySum;
                }
                else
                {
                    (*outputs++) = 0.;
                    (*outputs++) = 0.;
                    (*outputs)  = 0.;
                }
        }
        
        //! This method compute the energy vector.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and contains the channels samples and the minimum size must be the number of harmonics. The outputs array contains the vectors cartesian coordinates and the minimum size must be 3. The coordinates arrangement in the outputs array is energy abscissa, energy ordinate and energy height.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        inline void processEnergy(const T* inputs, T* outputs) noexcept
        {
            (*m_channels_square) = (*inputs) * (*inputs);
            for(ulong i = 1; i < Processor< Planewave<Hoa3d, T> >::getNumberOfPlanewaves(); i++)
                m_channels_square[i] = inputs[i] * inputs[i];
                
                T energySum = Signal<T>::vector_sum(Processor< Planewave<Hoa3d, T> >::getNumberOfPlanewaves(), m_channels_square);
                const T energyAbscissa  = Signal<T>::vectors_dot_product(Processor< Planewave<Hoa3d, T> >::getNumberOfPlanewaves(), m_channels_square, m_channels_abscissa);
                const T energyOrdinate  = Signal<T>::vectors_dot_product(Processor< Planewave<Hoa3d, T> >::getNumberOfPlanewaves(), m_channels_square, m_channels_ordinate);
                const T energyHeight    = Signal<T>::vectors_dot_product(Processor< Planewave<Hoa3d, T> >::getNumberOfPlanewaves(), m_channels_square, m_channels_height);
                
                if(energySum)
                {
                    (*outputs++) = energyAbscissa / energySum;
                    (*outputs++) = energyOrdinate / energySum;
                    (*outputs)   = energyHeight / energySum;
                }
                else
                {
                    (*outputs++) = 0.;
                    (*outputs++) = 0.;
                    (*outputs)   = 0.;
                }
        }
    };
    
#endif
}

#endif



