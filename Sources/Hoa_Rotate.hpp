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

#include "Hoa_Processor.hpp"

namespace hoa
{
    // ================================================================================ //
    // ROTATE //
    // ================================================================================ //
    
    //! @brief The rotate class rotates a sound field in the harmonics domain (2d available only).
    //! @details Rotate a sound field by weighting the harmonics depending on the rotation.
    template <Dimension D, typename T>
    class Rotate
    : public ProcessorHarmonics<D, T>
    {
    public:
        
        //! @brief Constructor.
        //! @param order The order (minimum 1).
        Rotate(const size_t order) noexcept;
        
        //! @brief Destructor.
        virtual ~Rotate() noexcept = 0;
        
        //! @brief This method sets the angle of the rotation around the z axis, the yaw value.
        //! @details The yaw is equivalent to a rotation around the z axis,
        //! the yaw value \f$\theta\f$ is in radian and should be between \f$0\f$ and \f$2\pi\f$.
        //! @param yaw The yaw value.
        virtual void setYaw(const T yaw) noexcept;
        
        //! @brief Get the angle of the rotation around the z axis, the yaw value.
        //! @details The method returns the angle of the rotation around the z axis, the yaw value \f$\theta\f$, in radian between \f$0\f$ and \f$2\pi\f$.
        //! @return The yaw value.
        virtual T getYaw() const noexcept;
        
        //! @brief This method performs the rotation.
        //! @details You should use this method for in-place or not-in-place processing and sample by sample.
        //! The inputs array and outputs array contains the spherical harmonics samples
        //! and the minimum size must be the number of harmonics.
        //! If \f$l = 0\f$
        //! \f[Y^{rotated}_{0,0}(\theta) = Y_{0,0}\f]
        //! else
        //! \f[Y^{rotated}_{l,-l}(\theta) = \sin{(\theta l)} \times Y_{l,l} + \cos{(\theta l)} \times Y_{l,-l}\f]
        //! and
        //! \f[Y^{rotated}_{l,l}(\theta) = \cos{(\theta l)} \times Y_{l,l} - \sin{(\theta l)} \times Y_{l,-l}\f]
        //! with \f$\theta\f$ the rotation in radian, \f$l\f$ the degree and \f$m\f$ the order.
        //! @param inputs The input array.
        //! @param outputs The output array.
        virtual void process(const T* inputs, T* outputs) noexcept;
    };
    
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    
    // ================================================================================ //
    // ROTATE 2D //
    // ================================================================================ //
    
    //! @brief 2d specialisation.
    template <typename T>
    class Rotate<Hoa2d, T>
    : public ProcessorHarmonics<Hoa2d, T>
    {
    public:
        
        //! @brief Constructor.
        //! @param order The order (minimum 1).
        Rotate(const size_t order)
        : ProcessorHarmonics<Hoa2d, T>(order)
        {}
        
        //! @brief Destructor.
        ~Rotate() = default;
        
        //! @brief This method sets the angle of the rotation around the z axis.
        //! @details The yaw is equivalent to a rotation around the z axis,
        //! the value is in radian and should be between 0 and 2π.
        //! @param yaw The yaw value.
        inline void setYaw(const T yaw)
        {
            m_yaw     =	yaw;
            m_cosx    = std::cos(m_yaw);
            m_sinx    = std::sin(m_yaw);
        }
        
        //! @brief Get the angle of the rotation around the z axis, the yaw value.
        //! @details Returns the angle of the rotation around the z axis.
        //! The yaw value is in radian between 0 and 2π.
        //! @return The yaw value.
        inline T getYaw() const noexcept
        {
            return math<T>::wrap_two_pi(m_yaw);
        }
        
        //! @brief This method performs the rotation.
        //! @details You should use this method for in-place or not-in-place processing and sample by sample.
        //! The inputs array and outputs array contains the spherical harmonics samples.
        //! The minimum size must be the number of harmonics.
        //! @param inputs   The input array.
        //! @param outputs  The output array.
        inline void process(const T* inputs, T* outputs) noexcept override
        {
            T cos_x = m_cosx;
            T sin_x = m_sinx;
            T tcos_x = cos_x;
            
            (*outputs++) = (*inputs++);
            T sig = (*inputs++);
            (*outputs++) = sin_x * (*inputs) + cos_x * sig;
            (*outputs++) = cos_x * (*inputs++) - sin_x * sig;
            for(size_t i = 2; i <= ProcessorHarmonics<Hoa2d, T>::getDecompositionOrder(); i++)
            {
                cos_x = tcos_x * m_cosx - sin_x * m_sinx;
                sin_x = tcos_x * m_sinx + sin_x * m_cosx;
                tcos_x = cos_x;
                sig = (*inputs++);
                (*outputs++) = sin_x * (*inputs) + cos_x * sig;
                (*outputs++) = cos_x * (*inputs++) - sin_x * sig;
            }
        }
        
    private:
        
        T m_yaw = 0.;
        T m_cosx = 0.;
        T m_sinx = 0.;
    };
    
#endif // DOXYGEN_SHOULD_SKIP_THIS
    
    // ================================================================================ //
    // ROTATE 3D //
    // ================================================================================ //
    
    //! @todo to be implemented..
}
