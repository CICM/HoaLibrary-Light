/*
// Copyright (c) 2012-2015 Eliott Paris, Julien Colafrancesco, Thomas Le Meur & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_VECTOR_LIGHT
#define DEF_HOA_VECTOR_LIGHT

#include "Planewaves.hpp"

namespace hoa
{
    //! The vector class computes the energy and the velocity vectors for a set of loudspeakers.
    /** The vector class compute the energy and the velocity vectors of a sound field for a set of channels. It is an useful tool to characterize the quality of the sound field restitution. For further information : Michael A. Gerzon, General metatheorie of auditory localization. Audio Engineering Society Preprint, 3306, 1992. This class retrieve the cartesian coordinates of the vectors.
     */
    template <Dimension D, typename T> class Vector : public Processor<D, T>::Planewaves
    {
    public:
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

    template <typename T> class Vector<Hoa2d, T> : public Processor<Hoa2d, T>::Planewaves
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
        Vector(const ulong numberOfChannels) noexcept : Processor<Hoa2d, T>::Planewaves(numberOfChannels)
        {
            m_channels_square   = new T[Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves()];
            m_channels_abscissa = new T[Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves()];
            m_channels_ordinate = new T[Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves()];
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

        //! This method pre-computes the necessary values to process.
        /**	You should use this method before calling the process methods and after changing the azimuth, the elevation or the offset of the channels.
         */
        inline void computeRendering() noexcept
        {
            for(ulong i = 0; i < Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves(); i++)
            {
                m_channels_abscissa[i] = Processor<Hoa2d, T>::Planewaves::getPlanewaveAbscissa(i);
                m_channels_ordinate[i] = Processor<Hoa2d, T>::Planewaves::getPlanewaveOrdinate(i);
            }
        }

        //! This method computes the energy and velocity vectors.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array contains the channels samples and the minimum size must be the number of channels. The outputs array contains the vectors cartesian coordinates and the minimum size must be 4. The coordinates arrangement in the outputs array is velocity abscissa, velocity ordinate, energy abscissa and energy ordinate.
            @param     inputs   The inputs array.
            @param     outputs  The outputs array.
         */
        inline void process(const T* inputs, T* outputs) noexcept override
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
            T veclocitySum = inputs[0];
            for(ulong i = 1; i < Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves(); i++)
                veclocitySum += inputs[i];

            const T velocityAbscissa = Signal<T>::dot(Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves(), inputs, m_channels_abscissa);
            const T velocityOrdinate = Signal<T>::dot(Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves(), inputs, m_channels_ordinate);
            if(veclocitySum)
            {
                outputs[0] = velocityAbscissa / veclocitySum;
                outputs[1] = velocityOrdinate / veclocitySum;
            }
            else
            {
                outputs[0] = outputs[1] = 0.;
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
            for(ulong i = 1; i < Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves(); i++)
                m_channels_square[i] = inputs[i] * inputs[i];

            T energySum = Signal<T>::sum(Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves(), m_channels_square);
            const T energyAbscissa = Signal<T>::dot(Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves(), m_channels_square, m_channels_abscissa);
            const T energyOrdinate = Signal<T>::dot(Processor<Hoa2d, T>::Planewaves::getNumberOfPlanewaves(), m_channels_square, m_channels_ordinate);

            if(energySum)
            {
                outputs[0] = energyAbscissa / energySum;
                outputs[1] = energyOrdinate / energySum;
            }
            else
            {
                outputs[0] = outputs[1] = 0.;
            }
        }
    };

    template <typename T> class Vector<Hoa3d, T> : public Processor<Hoa3d, T>::Planewaves
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
        Vector(const ulong numberOfChannels) noexcept : Processor<Hoa3d, T>::Planewaves(numberOfChannels)
        {
            m_channels_square   = new T[Processor<Hoa3d, T>::Planewaves::getNumberOfPlanewaves()];
            m_channels_abscissa = new T[Processor<Hoa3d, T>::Planewaves::getNumberOfPlanewaves()];
            m_channels_ordinate = new T[Processor<Hoa3d, T>::Planewaves::getNumberOfPlanewaves()];
            m_channels_height   = new T[Processor<Hoa3d, T>::Planewaves::getNumberOfPlanewaves()];
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

        //! This method pre-computes the necessary values to process.
        /**	You should use this method before calling the process methods and after changing the azimuth, the elevation or the offset of the channels.
         */
        inline void computeRendering() noexcept
        {
            for(ulong i = 0; i < Processor<Hoa3d, T>::Planewaves::getNumberOfPlanewaves(); i++)
            {
                m_channels_abscissa[i] = Processor<Hoa3d, T>::Planewaves::getPlanewaveAbscissa(i);
                m_channels_ordinate[i] = Processor<Hoa3d, T>::Planewaves::getPlanewaveOrdinate(i);
                m_channels_height[i]   = Processor<Hoa3d, T>::Planewaves::getPlanewaveHeight(i);
            }
        }

        //! This method compute the energy and velocity vectors.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array contains the channels samples and the minimum size must be the number of channels. The outputs array contains the vectors cartesian coordinates and the minimum size must be 6. The coordinates arrangement in the outputs array is velocity abscissa, velocity ordinate, velocity height, energy abscissa, energy ordinate and energy height.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        inline void process(const T* inputs, T* outputs) noexcept override
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
            T veclocitySum = inputs[0];
            for(ulong i = 1; i < Processor<Hoa3d, T>::Planewaves::getNumberOfPlanewaves(); i++)
                veclocitySum += inputs[i];

                const T velocityAbscissa    = Signal<T>::dot(Processor<Hoa3d, T>::Planewaves::getNumberOfPlanewaves(), inputs, m_channels_abscissa);
                const T velocityOrdinate    = Signal<T>::dot(Processor<Hoa3d, T>::Planewaves::getNumberOfPlanewaves(), inputs, m_channels_ordinate);
                const T velocityHeight      = Signal<T>::dot(Processor<Hoa3d, T>::Planewaves::getNumberOfPlanewaves(), inputs, m_channels_height);

                if(veclocitySum)
                {
                    outputs[0] = velocityAbscissa / veclocitySum;
                    outputs[1] = velocityOrdinate / veclocitySum;
                    outputs[2] = velocityHeight   / veclocitySum;
                }
                else
                {
                    outputs[0] = outputs[1] = outputs[2] = 0.;
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
            for(ulong i = 1; i < Processor<Hoa3d, T>::Planewaves::getNumberOfPlanewaves(); i++)
                m_channels_square[i] = inputs[i] * inputs[i];

                T energySum = Signal<T>::sum(Processor<Hoa3d, T>::Planewaves::getNumberOfPlanewaves(), m_channels_square);
                const T energyAbscissa  = Signal<T>::dot(Processor<Hoa3d, T>::Planewaves::getNumberOfPlanewaves(), m_channels_square, m_channels_abscissa);
                const T energyOrdinate  = Signal<T>::dot(Processor<Hoa3d, T>::Planewaves::getNumberOfPlanewaves(), m_channels_square, m_channels_ordinate);
                const T energyHeight    = Signal<T>::dot(Processor<Hoa3d, T>::Planewaves::getNumberOfPlanewaves(), m_channels_square, m_channels_height);

                if(energySum)
                {
                    outputs[0] = energyAbscissa / energySum;
                    outputs[1] = energyOrdinate / energySum;
                    outputs[2] = energyHeight   / energySum;
                }
                else
                {
                    outputs[0] = outputs[1] = outputs[2] = 0.;
                }
        }
    };

#endif
}

#endif



