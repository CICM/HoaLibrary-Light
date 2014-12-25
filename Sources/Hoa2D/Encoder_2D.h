/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_2D_ENCODER
#define DEF_HOA_2D_ENCODER

#include "Ambisonic_2D.h"

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
        Encoder(unsigned long order) : Ambisonic(order)
        {
            setAzimuth(0.);
        }
        
        //! The encoder destructor.
        /**	The encoder destructor free the memory.
         */
        ~Encoder()
        {
            ;
        }
        
        //! This method set the angle of azimuth.
        /**	The angle of azimuth in radian and you should prefer to use it between 0 and 2 π to avoid recursive wrapping of the value. The direction of rotation is counterclockwise. The 0 radian is π/2 phase shifted relative to a mathematical representation of a circle, then the 0 radian is at the "front" of the soundfield.
            @param     azimuth	The azimuth.
         */
        inline void setAzimuth(const double azimuth) noexcept
        {
            m_azimuth = wrap_twopi(azimuth);
            m_cosx    = std::cos(m_azimuth);
            m_sinx    = std::sin(m_azimuth);
        }
        
        //! Get the azimuth angle.
        /** The method returns the azimuth between 0 and 2π.
            @return     The azimuth.
         */
        inline double getAzimuth() const noexcept
        {
            return m_azimuth;
        }
        
        //! This method performs the encoding with single precision.
        /**	You should use this method for in-place or not-in-place processing and performs the encoding sample by sample. The outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         
         @param     input	The input sample.
         @param     outputs The output array.
         */
        inline void process(const float input, float* outputs) const noexcept
        {
            float cos_x = m_cosx, tcos_x = m_cosx, sin_x = m_sinx;
            outputs[0]  = input;                            // Hamonic [0,0]
            outputs[1]  = input * cos_x;                    // Hamonic [1,-1]
            outputs[2]  = input * sin_x;                    // Hamonic [1,1]
            for(unsigned long i = 3; i < m_number_of_harmonics; i += 2)
            {
                cos_x   = tcos_x * m_cosx - sin_x * m_sinx; // cos(x + b) = cos(x) * cos(b) - sin(x) * sin(b)
                sin_x   = tcos_x * m_sinx + sin_x * m_cosx; // sin(x + b) = cos(x) * sin(b) + sin(x) * cos(b)
                tcos_x  = cos_x;
                outputs[i]  = input * sin_x;                // Hamonic [l,-l]
                outputs[i+1]= input * cos_x;                // Hamonic [l,l]
            }
        }
        
        //! This method performs the encoding with double precision.
        /**	You should use this method for in-place or not-in-place processing and performs the encoding sample by sample. The outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         
            @param     input	The input sample.
            @param     outputs The output array.
         */
        inline void process(const double input, double* outputs) const noexcept
        {
            double cos_x = m_cosx, tcos_x = m_cosx, sin_x = m_sinx;
            outputs[0]  = input;                            // Hamonic [0,0]
            outputs[1]  = input * cos_x;                    // Hamonic [1,-1]
            outputs[2]  = input * sin_x;                    // Hamonic [1,1]
            for(unsigned long i = 3; i < m_number_of_harmonics; i += 2)
            {
                cos_x   = tcos_x * m_cosx - sin_x * m_sinx; // cos(x + b) = cos(x) * cos(b) - sin(x) * sin(b)
                sin_x   = tcos_x * m_sinx + sin_x * m_cosx; // sin(x + b) = cos(x) * sin(b) + sin(x) * cos(b)
                tcos_x  = cos_x;
                outputs[i]  = input * sin_x;                // Hamonic [l,-l]
                outputs[i+1]= input * cos_x;                // Hamonic [l,l]
            }
        }
    };
}

#endif



