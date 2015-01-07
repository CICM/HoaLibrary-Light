/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_2D_ENCODER
#define DEF_HOA_2D_ENCODER

#include "Ambisonic.h"

namespace Hoa2D
{
    //! The ambisonic encoder.
    /** The encoder should be used to encode a signal in the spherical harmonics domain depending of an order of decomposition. It allows to control the azimuth of the encoding.
     */
    class Encoder : public Ambisonic
    {
        
    private:
        
        double  m_azimuth;
        double  m_cosx;
        double  m_sinx;
        
    public:
        
        //! The encoder constructor.
        /**	The encoder constructor allocates and initialize the member values to computes spherical harmonics coefficients depending of a decomposition order. The order must be at least 1.
         
            @param     order	The order.
         */
        Encoder(unsigned int order);
        
        //! The encoder destructor.
        /**	The encoder destructor free the memory.
         */
        ~Encoder();
        
        //! This method set the angle of azimuth.
        /**	The angle of azimuth in radian and you should prefer to use it between 0 and 2 Pi to avoid recursive wrapping of the value. The direction of rotation is counterclockwise. The 0 radian is Pi/2 phase shifted relative to a mathematical representation of a circle, then the 0 radian is at the "front" of the soundfield.
         
            @param     azimuth	The azimuth.
         */
        void setAzimuth(const double azimuth);
        
        //! This method performs the encoding with single precision.
        /**	You should use this method for in-place or not-in-place processing and performs the encoding sample by sample. The outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         
            @param     input	The input sample.
            @param     outputs The output array.
         */
        
        //! Get the azimuth angle
        /** The method returns the last angle of encoding between 0 and 2π.
		 
            @return     The azimuth.
         */
        inline double getAzimuth() const {return m_azimuth;};
        
        //! This method performs the encoding with single precision.
        /**	You should use this method for in-place or not-in-place processing and performs the encoding sample by sample. The outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         
         @param     input	The input sample.
         @param     outputs The output array.
         */
        void process(const float input, float* outputs);
        
        //! This method performs the encoding with double precision.
        /**	You should use this method for in-place or not-in-place processing and performs the encoding sample by sample. The outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         
            @param     input	The input sample.
            @param     outputs The output array.
         */
        void process(const double input, double* outputs);
    };
}

#endif



