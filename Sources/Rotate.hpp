/*
// Copyright (c) 2012-2015 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_ROTATE_LIGHT
#define DEF_HOA_ROTATE_LIGHT

#include "Harmonics.hpp"

namespace hoa
{
    //! The rotate class rotates a sound field in the harmonics domain (2d available only).
    /** The rotate should be used to rotate a sound field by wheighting the harmonics depending on the rotation.
     */
    template <Dimension D, typename T> class Rotate : public  Harmonic<D, T>::Processor
    {
    public:
        
        //! The rotate constructor.
        /**	The rotate constructor allocates and initialize the member values. The order must be at least 1.
         @param     order	The order.
         */
        Rotate(const ulong order) noexcept  = 0;
        
        //! The Rotate destructor.
        /**	The Rotate destructor free the memory.
         */
        virtual ~Rotate() noexcept = 0;
        
        //! This method sets the angle of the rotation around the z axis, the yaw value,
        /** The yaw is equivalent to a rotation around the z axis, the value is in radian and should be between 0 and 2π.
         @param     yaw The yaw value.
         */
        virtual void setYaw(const T yaw) noexcept = 0;
        
        //! Get the angle of the rotation around the z axis, the yaw value.
        /** The method returns the angle of the rotation around the z axis, the yaw value, in radian between 0 and 2π.
         @return     The yaw value.
         */
        virtual T getYaw() const noexcept = 0;
        
        //! This method performs the rotation.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The input array.
         @param     outputs  The output array.
         */
        virtual void process(const T* inputs, T* outputs) const noexcept = 0;
    };
    
    template <typename T> class Rotate<Hoa2d, T> : public Harmonic<Hoa2d, T>::Processor
    {
    private:
        T   m_yaw;
        T   m_cosx;
        T   m_sinx;
    public:
        
        //! The rotate constructor.
        /**	The rotate constructor allocates and initialize the member values to computes spherical harmonics rotation depending on a decomposition order. The order must be at least 1.
         @param     order	The order.
         */
        Rotate(const ulong order) noexcept :
        Harmonic<Hoa2d, T>::Processor(order)
        {
            ;
        }
        
        //! The Rotate destructor.
        /**	The Rotate destructor free the memory.
         */
        ~Rotate() noexcept
        {
            ;
        }
        
        //! This method sets the angle of the rotation around the z axis, the yaw value,
        /** The yaw is equivalent to a rotation around the z axis, the value is in radian and should be between 0 and 2π.
         @param     yaw The yaw value.
         */
        inline void setYaw(const T yaw) noexcept
        {
            m_yaw     =	yaw;
            m_cosx    = std::cos(m_yaw);
            m_sinx    = std::sin(m_yaw);
        }
        
        //! Get the angle of the rotation around the z axis, the yaw value.
        /** The method returns the angle of the rotation around the z axis, the yaw value, in radian between 0 and 2π.
         @return     The yaw value.
         */
        inline T getYaw() const noexcept
        {
            return Math<T>::wrap_twopi(m_yaw);
        };
        
        //! This method performs the rotation.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The input array.
         @param     outputs  The output array.
         */
        inline void process(const T* inputs, T* outputs) const noexcept
        {
            T cos_x = m_cosx;
            T sin_x = m_sinx;
            T tcos_x = cos_x;
            
            (*outputs++) = (*inputs++);
            T sig = (*inputs++);
            (*outputs++) = sin_x * (*inputs) + cos_x * sig;
            (*outputs++) = cos_x * (*inputs++) - sin_x * sig;
            for(ulong i = 2; i <= Harmonic<Hoa2d, T>::Processor::getDecompositionOrder(); i++)
            {
                cos_x = tcos_x * m_cosx - sin_x * m_sinx;
                sin_x = tcos_x * m_sinx + sin_x * m_cosx;
                tcos_x = cos_x;
                sig = (*inputs++);
                (*outputs++) = sin_x * (*inputs) + cos_x * sig;
                (*outputs++) = cos_x * (*inputs++) - sin_x * sig;
            }
        }
    };
}

#endif



