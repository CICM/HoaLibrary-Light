/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef __DEF_HOA_3D_ENCODER__
#define __DEF_HOA_3D_ENCODER__

#include "Ambisonic_3D.h"

namespace Hoa3D
{
    //! The ambisonic encoder.
    /** The encoder should be used to encode a signal in the spherical harmonics domain depending of an order of decomposition. It allows to control the azimuth and the elevation of the encoding. If you want to spatialize with distance compensation, you should use the Map class.
        
        @see Map
     */
    class Encoder : public Ambisonic
    {
        
    private:
        double  m_azimuth;
        double  m_elevation;
        double  m_cos_phi;
        double  m_sin_phi;
        double  m_cos_theta;
        double  m_sqrt_rmin;
    public:
        
        //! The encoder constructor.
        /**	The encoder constructor allocates and initialize the member values to computes spherical harmonics coefficients depending of a decomposition order. The order must be at least 1.
         
            @param     order	The order.
         */
        Encoder(unsigned long order) noexcept : Ambisonic(order)
        {
            setAzimuth(0.);
            setElevation(0.);
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
         @see       setElevation()
         */
        inline void setAzimuth(const double azimuth) noexcept
        {
            m_azimuth = wrap_twopi(azimuth);
            m_cos_phi = std::cos(m_azimuth);
            m_sin_phi = std::sin(m_azimuth);
        }
        
        //! Get the azimuth angle
        /** The method returns the azimuth between 0 and 2π.
         @return     The azimuth.
         */
        inline double getAzimuth() const noexcept
        {
            return m_azimuth;
        }
        
        //! This method set the angle of elevation.
        /**	The angle of elevation in radian and you should prefer to use it between -π and π to avoid recursive wrapping of the value. The direction of rotation is from bottom to the top. The 0 radian is centered at the "front" of the soundfield, then π/2 is at the top, -π/2 is at the bottom and π is behind. Note that if the angle of elevation is between π/2 and 3*π/2, the azimuth is reversed.
         @param     elevation The elevation.
         @see       setAzimutHamonic [)
         */
        inline void setElevation(const double elevation) noexcept
        {
            m_elevation = wrap(elevation, -HOA_PI, HOA_PI);
            m_cos_theta = std::cos(HOA_PI2 + m_elevation);
            m_sqrt_rmin = std::sqrt(1 - m_cos_theta * m_cos_theta);
        }
        
        //!	Get the elevation angle
        /** The method returns the elevation between -π and π.
         @return     The elevation.
         */
        inline double getElevation()  const noexcept
        {
            return m_elevation;
        }
        
        //! Get the normalization of an harmonic.
        /** The method returns the normalization of an harmonic.
         @param index The index of the harmonic.
         @return The normalization.
         */
        inline double getNormalization(const unsigned long index) const
        {
            return 1.;
        }
        
        //! This method performs the encoding with single precision.
        /**	You should use this method for in-place or not-in-place processing and performs the encoding sample by sample. The outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics. For the elevation, the function uses three recurrence formulas :
         \f[P(l, l)(x) = (-1)^l \times (2l - 1)!! \times (1 - x^2)^{0.5l}\f]
         \f[P(l + 1, l)(x) = x \times (2l + 1) \times P(l, l)\f]
         \f[P(l + 1, m)(x) = \frac{(2l + 1) \times x \times P(m, l) - (l + m) \times P(m, l - 1)}{(l - m + 1)}\f]
         with \f$0 \leq l\f$ and \f$-l \leq m \leq +l\f$ and \f[P(0, 0] = 1\f]
         @param     input    The input sample.
         @param     outputs  The outputs array.
         */
        void process(const float input, float* outputs) const noexcept
        {
            float leg1, leg2, leg3;
            float cos_xphi = m_cos_phi, tcos_xphi = m_cos_phi, sin_xphi = m_sin_phi;
            outputs[0] = 1.f;                       // Hamonic [0,0]    = 1
            outputs[1] = -m_sqrt_rmin * sin_xphi;   // Hamonic [1,-1]   = -√(1-x²) * sin(ϕ)
            outputs[2] = m_cos_theta;               // Hamonic [1,0]    = x
            outputs[3] = -m_sqrt_rmin * cos_xphi;   // Hamonic [1,1]    = -√(1-x²) * cos(ϕ)
            
            for(unsigned long i = 4, l = 2, h = 2 * l + 1; l <= m_order_of_decomposition; l++, h = 2 * l + 1)
            {
                // P(l+1,l+1)   = -(2l+1)√(1-x²)P(l,l)
                leg1  = -m_sqrt_rmin * h * outputs[i-1];
                // P(l,l+1)     = x(2l+1)P(l,l)
                leg2  = m_cos_theta  * h * outputs[i-1];
                outputs[i++]= leg1;                 // Hamonic [l,-l]      = P(l+1,l+1)
                outputs[i++]= leg2;                 // Hamonic [l,-l+1]    = P(l,l+1)
             
                // From m = l - 2 to m = l + 2 without m = 0;
                for(unsigned long m = l - 2, index = i - 3; m; m--, index--)
                {
                    // Use of the symmetry of the harmonics
                    // P(l+1,m) = (x(2l+1)P(l,m) - (l+m)P(l-1,m))/(l-m+1)
                    leg3  = m_cos_theta * h - (l + m) * outputs[index] / (l - m + 1);
                    outputs[i]      = leg3;         // Hamonic [l,-m]      = P(l+1,m)
                    outputs[i+m*2]  = leg3;         // Hamonic [l,-m]      = P(l+1,m)
                    i++;
                }
                // Hamonic [l,0) = P(l+1,0) = (x(2l+1)P(l,0) - l * P(l-1,0))/(l+1)
                unsigned long index = i - h;
                //post("Hamonic [%lu,0) %i",l, i);
                outputs[i] = m_cos_theta * h * outputs[index] - l * outputs[index - h + 2] / (l + 1);
                i += l-1;
                // Use of the symmetry of the harmonics
                outputs[i++] = leg2;                // Hamonic [l,l-1]     = P(l,l+1)
                outputs[i++] = leg1;                // Hamonic [l,l]       = P(l+1,l+1)
    
                /*
                float cos_x = m_cos_phi;
                float sin_x = m_sin_phi;
                float tcos_x = cos_x;
                
                for(unsigned long j = 1; j < i * 2 + 1; j += 2)
                {
                    outputs[i+j] = input * sin_x;
                    outputs[i+j+1] = input * cos_x;
                    cos_x = tcos_x * m_cos_phi - sin_x * m_sin_phi; // cos(x + b) = cos(x) * cos(b) - sin(x) * sin(b)
                    sin_x = tcos_x * m_sin_phi + sin_x * m_cos_phi; // sin(x + b) = cos(x) * sin(b) + sin(x) * cos(b)
                    tcos_x = cos_x;
                }*/
            }
            // Multiplication of the factors by the input
            cblas_sscal(m_number_of_harmonics, input, outputs, 1);
        }
        
        //! This method performs the encoding with double precision.
        /**	You should use this method for in-place or not-in-place processing and performs the encoding sample by sample. The outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         
            @param     input	The input sample.
            @param     outputs  The outputs array.
         */
        void process(const double input, double* outputs) const noexcept
        {
            ;
        }
    };
}

#endif



