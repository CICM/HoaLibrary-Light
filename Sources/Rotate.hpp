/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_ROTATE
#define DEF_HOA_ROTATE

#include "Harmonics.hpp"

namespace hoa
{
    //! The ambisonic encoder.
    /** The encoder should be used to encode a source in the spherical harmonics domain depending of an order of decomposition. It allows to control the azimuth of the source.
     */
    template <Dimension D, typename T> class Rotate : public Harmonic<D>::Processor
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
        Harmonic<D>::Processor(order)
        {
            ;
        }
        
        //! The Rotate destructor.
        /**	The Rotate destructor free the memory.
         */
        ~Rotate()
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
            return wrap_twopi(m_yaw);
        };
        
        //! This method performs the rotation.
        /**	You should use this method for in-place or not-in-place processing and performs the rotation sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
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
            for(ulong i = 2; i <= Harmonic<D>::Processor::getDecompositionOrder(); i++)
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



