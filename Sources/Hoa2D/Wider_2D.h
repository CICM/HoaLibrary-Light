/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_2D_WIDER
#define DEF_HOA_2D_WIDER

#include "Ambisonic_2D.h"

namespace Hoa2D
{
    //! The ambisonic wider.
    /** The wider should be used to widen the sound propagation with fractional order simulution. The sound field precision depends to the decomposition order. The zero decomposition order has 1 omnidirectionnal harmonic and all the sounds seem to come from all the directions. While the order increases, the number of harmonics increases, the lobes of an encoded sounds narrow and the origin of the sounds is more accurate. Then fractional order can be used to decrease the sound field precision and to wide the sound field propagation. 
     */
    class Wider : public Ambisonic
    {
    private:
        double  m_factor;
        double  m_gain;
        
    public:
        
        //! The wider constructor.
        /**	The wider constructor allocates and initialize the member values to computes circular harmonics weighted coefficients depending of a decomposition order. The order must be at least 1.
         
            @param     order	The order.
         */
        Wider(unsigned long order) noexcept : Ambisonic(order)
        {
            setWideningValue(0.);
        }
        
        //! The wider destructor.
        /**	The wider destructor free the memory.
         */
        ~Wider()
        {
            ;
        }
        
        //! This method set the widening value.
        /**	The widening value is clipped between 0 and 1. At 1, the sound field has no changes. At 0, all the sound field is omnidirectionnal, only the harmonic [0 0] remains. From 0 to 1, the circular hamronics appears in logarithmic way to linearly increase the sound field precision.
         
            @param     value The widening value.
         */
        inline void setWideningValue(const double value) noexcept
        {
            m_factor = (1. - clip_minmax(value, 0., 1.)) * HOA_PI;
            m_gain   = (std::sin(clip_minmax((m_factor - 0.5f) * HOA_PI, -HOA_PI2, HOA_PI2)) + 1.) * 0.5;
        }
        
        //! This method retreive the widening value.
        /**	The method returns the widening value.
         
            @return     The widening value.
         */
        inline double getWideningValue() const noexcept
        {
            return m_factor / HOA_PI;
        }
        
        //! This method performs the widening with single precision.
        /**	You should use this method for in-place or not-in-place processing and performs the widening sample by sample. The inputs array and outputs array contains the circular harmonics samples and the minimum size must be the number of harmonics.
         
            @param     inputs   The inputs array.
            @param     outputs  The outputs array.
         */
        inline void process(const float* inputs, float* outputs) const noexcept
        {
            outputs[0] = inputs[0] * (m_gain * m_order_of_decomposition + 1.f);
            for(unsigned long i = 2, j = 1; i < m_number_of_harmonics; i += 2, j++)
            {
                const float factor  = (std::cos(clip_max(m_factor * j, HOA_PI)) + 1.f) * 0.5f * (m_gain * (m_order_of_decomposition - j) + 1.f);
                outputs[i-1]        = inputs[i-1] * factor;
                outputs[i]          = inputs[i] * factor;
            }
        }
        
        //! This method performs the widening with double precision.
        /**	You should use this method for in-place or not-in-place processing and performs the widening sample by sample. The inputs array and outputs array contains the circular harmonics samples and the minimum size must be the number of harmonics.
         
            @param     inputs   The inputs array.
            @param     outputs  The outputs array.
         */
        inline void process(const double* inputs, double* outputs) const noexcept
        {
            outputs[0] = inputs[0] * (m_gain * m_order_of_decomposition + 1.f);
            for(unsigned long i = 2, j = 1; i < m_number_of_harmonics; i += 2, j++)
            {
                const double factor  = (std::cos(clip_max(m_factor * j, HOA_PI)) + 1.) * 0.5 * (m_gain * (m_order_of_decomposition - j) + 1.);
                outputs[i-1]        = inputs[i-1] * factor;
                outputs[i]          = inputs[i] * factor;
            }
        }
    };
}

#endif



