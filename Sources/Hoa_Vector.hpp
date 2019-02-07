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

#include "Hoa_Planewaves.hpp"

namespace hoa
{
    // ================================================================================ //
    // VECTOR //
    // ================================================================================ //
    
    //! @brief Computes the energy and the velocity vectors for a set of loudspeakers.
    //! @details Computes the energy and the velocity vectors of a sound field for a set of channels.
    //! It is a useful tool to characterize the quality of the sound field restitution.
    //! For further information see :
    //! Michael A. Gerzon, General metatheorie of auditory localization.
    //! Audio Engineering Society Preprint, 3306, 1992.
    //! This class retrieve the cartesian coordinates of the vectors.
    template <Dimension D, typename T>
    class Vector
    : public ProcessorPlanewaves<D, T>
    {
    public:
        
        //! @brief Constructor.
        //! @param channels	The number of channels (minimum 1).
        Vector(const size_t channels) noexcept;

        //! @brief Destructor.
        ~Vector() = default;

        //! @brief This method pre-computes the necessary values to process.
        //! @details You should use this method before calling the process methods
        //! and after changing the azimuth, the elevation or the offset of the channels.
        virtual void computeRendering() noexcept = 0;

        //! @brief This method computes the energy and velocity vectors.
        //! @details You should use this method for in-place or not-in-place processing and sample by sample.
        //! The inputs array contains the channels samples and the minimum size must be the number of channels.
        //! The outputs array contains the vectors cartesian coordinates
        //! and the minimum size must be 4 for 2d and 6 for 3d.
        //! The coordinates arrangement in the outputs array is :
        //! velocity abscissa, velocity ordinate, (velocity height),
        //! energy abscissa and energy ordinate (and energy height).
        //! @param inputs   The inputs array.
        //! @param outputs  The outputs array.
        virtual void process(const T* inputs, T* outputs) noexcept = 0;

        //! @brief This method computes the velocity vector.
        //! @details You should use this method for in-place or not-in-place processing and sample by sample.
        //! The inputs array and contains the channels samples and the minimum size must be the number of channels.
        //! The outputs array contains the vectors cartesian coordinates
        //! and the minimum size must be 2 for 2d and 3 for 3d.
        //! The coordinates arrangement in the outputs array is :
        //! velocity abscissa and velocity ordinate (and velocity height).
        //! @param inputs   The inputs array.
        //! @param outputs  The outputs array.
        virtual void processVelocity(const T* inputs, T* outputs) noexcept = 0;

        //! @brief This method computes the energy vector.
        //! @details You should use this method for in-place or not-in-place processing and sample by sample.
        //! The inputs array and contains the channels samples
        //! and the minimum size must be the number of harmonics.
        //! The outputs array contains the vectors cartesian coordinates
        //! and the minimum size must be 2 for 2d and 3 for 3d.
        //! The coordinates arrangement in the outputs array is :
        //! energy abscissa and energy ordinate (and energy height).
        //! @param inputs   The inputs array.
        //! @param outputs  The outputs array.
        virtual void processEnergy(const T* inputs, T* outputs) noexcept = 0;
    };

#ifndef DOXYGEN_SHOULD_SKIP_THIS

    // ================================================================================ //
    // VECTOR 2D //
    // ================================================================================ //
    
    template <typename T>
    class Vector<Hoa2d, T>
    : public ProcessorPlanewaves<Hoa2d, T>
    {
    public:

        //! @brief The vector constructor.
        //! @details The vector constructor allocates and initialize the member values to computes vectors. The number of channels must be at least 1.
        //! @param channels	The number of channels.
        Vector(const size_t channels)
        : ProcessorPlanewaves<Hoa2d, T>(channels)
        {
            const auto size = ProcessorPlanewaves<Hoa2d, T>::getNumberOfPlanewaves();
            m_channels_square   = Signal<T>::alloc(size);
            m_channels_abscissa = Signal<T>::alloc(size);
            m_channels_ordinate = Signal<T>::alloc(size);
        }

        //! @brief Destructor.
        ~Vector()
        {
            Signal<T>::free(m_channels_square);
            Signal<T>::free(m_channels_abscissa);
            Signal<T>::free(m_channels_ordinate);
        }

        //! @brief This method pre-computes the necessary values to process.
        //! @details You should use this method before calling the process methods
        //! and after changing the azimuth, the elevation or the offset of the channels.
        inline void computeRendering() noexcept
        {
            for(size_t i = 0; i < ProcessorPlanewaves<Hoa2d, T>::getNumberOfPlanewaves(); i++)
            {
                m_channels_abscissa[i] = ProcessorPlanewaves<Hoa2d, T>::getPlanewaveAbscissa(i);
                m_channels_ordinate[i] = ProcessorPlanewaves<Hoa2d, T>::getPlanewaveOrdinate(i);
            }
        }

        //! @brief This method computes the energy and velocity vectors.
        //! @details You should use this method for in-place or not-in-place processing and sample by sample.
        //! The inputs array contains the channels samples
        //! and the minimum size must be the number of channels.
        //! The outputs array contains the vectors cartesian coordinates and the minimum size must be 4.
        //! The coordinates arrangement in the outputs array is:
        //! velocity abscissa, velocity ordinate, energy abscissa and energy ordinate.
        //! @param inputs   The inputs array.
        //! @param outputs  The outputs array.
        inline void process(const T* inputs, T* outputs) noexcept override
        {
            processVelocity(inputs, outputs);
            processEnergy(inputs, outputs+2);
        }

        //! @brief This method computes the velocity vector.
        //! @details You should use this method for in-place or not-in-place processing and sample by sample.
        //! The inputs array and contains the channels samples and the minimum size must be the number of channels.
        //! The outputs array contains the vectors cartesian coordinates and the minimum size must be 2.
        //! The coordinates arrangement in the outputs array is: velocity abscissa and velocity ordinate.
        //! @param     inputs   The inputs array.
        //! @param     outputs  The outputs array.
        inline void processVelocity(const T* inputs, T* outputs) noexcept
        {
            const auto planewaves = ProcessorPlanewaves<Hoa2d, T>::getNumberOfPlanewaves();
            
            auto sum = inputs[0];
            for(size_t i = 1; i < planewaves; i++)
            {
                sum += inputs[i];
            }

            if(sum)
            {
                const auto abscissa = Signal<T>::dot(planewaves, inputs, m_channels_abscissa);
                const auto ordinate = Signal<T>::dot(planewaves, inputs, m_channels_ordinate);
                
                outputs[0] = abscissa / sum;
                outputs[1] = ordinate / sum;
            }
            else
            {
                outputs[0] = outputs[1] = 0.;
            }
        }

        //! @brief This method computes the energy vector.
        //! @details You should use this method for in-place or not-in-place processing and sample by sample.
        //! The inputs array and contains the channels samples
        //! and the minimum size must be the number of harmonics.
        //! The outputs array contains the vectors cartesian coordinates and the minimum size must be 2.
        //! The coordinates arrangement in the outputs array is: energy abscissa and energy ordinate.
        //! @param inputs   The inputs array.
        //! @param outputs  The outputs array.
        inline void processEnergy(const T* inputs, T* outputs) noexcept
        {
            const auto planewaves = ProcessorPlanewaves<Hoa2d, T>::getNumberOfPlanewaves();
            
            (*m_channels_square) = (*inputs) * (*inputs);
            
            for(size_t i = 1; i < planewaves; i++)
            {
                m_channels_square[i] = inputs[i] * inputs[i];
            }

            const auto sum = Signal<T>::sum(planewaves, m_channels_square);

            if(sum)
            {
                const auto abscissa = Signal<T>::dot(planewaves, m_channels_square, m_channels_abscissa);
                const auto ordinate = Signal<T>::dot(planewaves, m_channels_square, m_channels_ordinate);
                
                outputs[0] = abscissa / sum;
                outputs[1] = ordinate / sum;
            }
            else
            {
                outputs[0] = outputs[1] = 0.;
            }
        }
        
    private:
        
        T* m_channels_square = nullptr;
        T* m_channels_abscissa = nullptr;
        T* m_channels_ordinate = nullptr;
    };
    
    // ================================================================================ //
    // VECTOR 3D //
    // ================================================================================ //

    template <typename T>
    class Vector<Hoa3d, T>
    : public ProcessorPlanewaves<Hoa3d, T>
    {
    public:

        //! @brief Constructor.
        //! @param channels	The number of channels.
        Vector(const size_t channels)
        : ProcessorPlanewaves<Hoa3d, T>(channels)
        {
            const auto size = ProcessorPlanewaves<Hoa3d, T>::getNumberOfPlanewaves();
            m_channels_square   = Signal<T>::alloc(size);
            m_channels_abscissa = Signal<T>::alloc(size);
            m_channels_ordinate = Signal<T>::alloc(size);
            m_channels_height   = Signal<T>::alloc(size);
        }

        //! @brief Destructor.
        ~Vector()
        {
            Signal<T>::free(m_channels_square);
            Signal<T>::free(m_channels_abscissa);
            Signal<T>::free(m_channels_ordinate);
            Signal<T>::free(m_channels_height);
        }

        //! @brief This method pre-computes the necessary values to process.
        //! @details You should use this method before calling the process methods
        //! and after changing the azimuth, the elevation or the offset of the channels.
        inline void computeRendering()
        {
            for(size_t i = 0; i < ProcessorPlanewaves<Hoa3d, T>::getNumberOfPlanewaves(); i++)
            {
                m_channels_abscissa[i] = ProcessorPlanewaves<Hoa3d, T>::getPlanewaveAbscissa(i);
                m_channels_ordinate[i] = ProcessorPlanewaves<Hoa3d, T>::getPlanewaveOrdinate(i);
                m_channels_height[i]   = ProcessorPlanewaves<Hoa3d, T>::getPlanewaveHeight(i);
            }
        }

        //! @brief This method compute the energy and velocity vectors.
        //! @details You should use this method for in-place or not-in-place processing and sample by sample.
        //! The inputs array contains the channels samples and the minimum size must be the number of channels.
        //! The outputs array contains the vectors cartesian coordinates and the minimum size must be 6.
        //! The coordinates arrangement in the outputs array is:
        //! velocity abscissa, velocity ordinate, velocity height,
        //! energy abscissa, energy ordinate and energy height.
        //! @param inputs   The inputs array.
        //! @param outputs  The outputs array.
        inline void process(const T* inputs, T* outputs) noexcept override
        {
            processVelocity(inputs, outputs);
            processEnergy(inputs, outputs+3);
        }

        //! @brief This method compute the velocity vector.
        //! @details You should use this method for in-place or not-in-place processing and sample by sample.
        //! The inputs array and contains the channels samples and the minimum size must be the number of channels.
        //! The outputs array contains the vectors cartesian coordinates and the minimum size must be 3.
        //! The coordinates arrangement in the outputs array is:
        //! velocity abscissa, velocity ordinate and velocity height.
        //! @param inputs   The inputs array.
        //! @param outputs  The outputs array.
        inline void processVelocity(const T* inputs, T* outputs) noexcept
        {
            const auto planewaves = ProcessorPlanewaves<Hoa3d, T>::getNumberOfPlanewaves();
            
            auto sum = inputs[0];
            for(size_t i = 1; i < planewaves; i++)
            {
                sum += inputs[i];
            }
            
            if(sum)
            {
                const auto abscissa = Signal<T>::dot(planewaves, inputs, m_channels_abscissa);
                const auto ordinate = Signal<T>::dot(planewaves, inputs, m_channels_ordinate);
                const auto height   = Signal<T>::dot(planewaves, inputs, m_channels_height);
                
                outputs[0] = abscissa / sum;
                outputs[1] = ordinate / sum;
                outputs[2] = height / sum;
            }
            else
            {
                outputs[0] = outputs[1] = outputs[2] = 0.;
            }
        }

        //! @brief This method compute the energy vector.
        //! @details You should use this method for in-place or not-in-place processing and sample by sample.
        //! The inputs array and contains the channels samples
        //! and the minimum size must be the number of harmonics.
        //! The outputs array contains the vectors cartesian coordinates and the minimum size must be 3.
        //! The coordinates arrangement in the outputs array is:
        //! energy abscissa, energy ordinate and energy height.
        //! @param inputs   The inputs array.
        //! @param outputs  The outputs array.
        inline void processEnergy(const T* inputs, T* outputs) noexcept
        {
            const auto planewaves = ProcessorPlanewaves<Hoa3d, T>::getNumberOfPlanewaves();
            
            (*m_channels_square) = (*inputs) * (*inputs);
            for(size_t i = 1; i < planewaves; i++)
            {
                m_channels_square[i] = inputs[i] * inputs[i];
            }
            
            auto sum = Signal<T>::sum(planewaves, m_channels_square);
            
            if(sum)
            {
                const T abscissa  = Signal<T>::dot(planewaves, m_channels_square, m_channels_abscissa);
                const T ordinate  = Signal<T>::dot(planewaves, m_channels_square, m_channels_ordinate);
                const T height    = Signal<T>::dot(planewaves, m_channels_square, m_channels_height);
                
                outputs[0] = abscissa / sum;
                outputs[1] = ordinate / sum;
                outputs[2] = height   / sum;
            }
            else
            {
                outputs[0] = outputs[1] = outputs[2] = 0.;
            }
        }
        
    private:
        
        T* m_channels_square = nullptr;
        T* m_channels_abscissa = nullptr;
        T* m_channels_ordinate = nullptr;
        T* m_channels_height = nullptr;
    };

#endif
    
}
