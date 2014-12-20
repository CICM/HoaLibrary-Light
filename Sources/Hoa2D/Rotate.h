/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_2D_ROTATE
#define DEF_HOA_2D_ROTATE

#include "Ambisonic.h"

namespace Hoa2D
{
    //! The ambisonic rotation.
    /** The rotate should be used to rotate the ambisonic soundfield.
     */
    class Rotate : public Ambisonic
    {
    private:
		
		double  m_yaw;
        double  m_cosx;
        double  m_sinx;
        
    public:
        
        //! The rotate constructor.
        /**	The rotate constructor allocates and initialize the member values to computes spherical harmonics rotation depending on a decomposition order. The order must be at least 1.
         
            @param     order	The order.
         */
        Rotate(unsigned int order);
        
        //! The Rotate destructor.
        /**	The Rotate destructor free the memory.
         */
        ~Rotate();
		
		//! This method sets the angle of the rotation around the z axis, the yaw value,
        /** The yaw is equivalent to a rotation around the z axis, the value is in radian and should be between 0 and 2π.
		 
            @param     value The yaw value.
         */
        void setYaw(const double value);
		
		//! Get the angle of the rotation around the z axis, the yaw value.
        /** The method returns the angle of the rotation around the z axis, the yaw value, in radian between 0 and 2π.
		 
            @return     The yaw value.
         */
        inline double getYaw() const {return m_yaw;};
        
        //! This method performs the rotation with single precision.
        /**	You should use this method for in-place or not-in-place processing and performs the rotation sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         
            @param     inputs   The input array.
            @param     outputs  The output array.
         */
        void process(const float* inputs, float* outputs);
        
        //! This method performs the rotation with double precision.
        /**	You should use this method for in-place or not-in-place processing and performs the rotation sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         
            @param     inputs   The input array.
            @param     outputs  The output array.
         */
        void process(const double* inputs, double* outputs);
    };
	
}

#endif



