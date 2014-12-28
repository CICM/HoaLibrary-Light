/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_2D_ROTATE
#define DEF_HOA_2D_ROTATE

#include "Ambisonic_2D.hpp"

namespace hoa
{
    //! The ambisonic rotation.
    /** The rotate should be used to rotate the ambisonic soundfield.
     */
    template <typename Type> class Rotate : public Ambisonic2D<Type>
    {
    private:
		Type    m_yaw;
        Type    m_cosx;
        Type    m_sinx;
        
    public:
        
        //! The rotate constructor.
        /**	The rotate constructor allocates and initialize the member values to computes spherical harmonics rotation depending on a decomposition order. The order must be at least 1.
            @param     order	The order.
         */
        Rotate(unsigned long order) noexcept : Ambisonic2D<Type>(order)
        {
            setYaw(0.);
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
		 
            @param     value The yaw value.
         */
        inline void setYaw(const Type value) noexcept
        {
            m_yaw     =	wrap_twopi(value);
            m_cosx    = std::cos(m_yaw);
            m_sinx    = std::sin(m_yaw);
        }
		
		//! Get the angle of the rotation around the z axis, the yaw value.
        /** The method returns the angle of the rotation around the z axis, the yaw value, in radian between 0 and 2π.
            @return     The yaw value.
         */
        inline Type getYaw() const noexcept
        {
            return m_yaw;
        };
        
        //! This method performs the rotation.
        /**	You should use this method for in-place or not-in-place processing and performs the rotation sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
            @param     inputs   The input array.
            @param     outputs  The output array.
         */
        inline void process(const Type* inputs, Type* outputs) const noexcept 
        {
            Type cos_x = m_cosx;
            Type sin_x = m_sinx;
            Type tcos_x = cos_x;
            Type sig;
            (*outputs++) = (*inputs++);
            for(unsigned long i = 2; i < Ambisonic<Type>::m_number_of_harmonics; i += 2)
            {
                sig = (*inputs++);
                (*outputs++) = sin_x * (*inputs) + cos_x * sig;
                (*outputs++) = cos_x * (*inputs++) - sin_x * sig;
                cos_x = tcos_x * m_cosx - sin_x * m_sinx;
                sin_x = tcos_x * m_sinx + sin_x * m_cosx;
                tcos_x = cos_x;
            }
        }
    };
	
}

#endif



