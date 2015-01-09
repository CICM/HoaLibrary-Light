/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_WIDER_LIGHT
#define DEF_HOA_WIDER_LIGHT

#include "Harmonics.hpp"

namespace hoa
{
    template <Dimension D, typename T> class Wider;
    
    //! The ambisonic encoder with distance compensation.
    /** The encoder with distance compensation should be used to encode a source in the spherical harmonics domain depending of an order of decomposition. It allows to control the azimuth and the radius of the source.
     */
    template <typename T> class Wider<Hoa2d, T> : public Harmonic<Hoa2d, T>::Processor
    {
    private:
        T   m_widening;
        T   m_gain;
        T   m_factor;
    public:
        
        //! The encoder constructor.
        /**	The encoder constructor allocates and initialize the member values to computes circular harmonics coefficients for the encoding. The order must be at least 1.
         @param     order	The order.
         */
        Wider(const ulong order) :
        Harmonic<Hoa2d, T>::Processor(order)
        {
            setWidening(1.);
        }
        
        //! The encoder destructor.
        /**	The encoder destructor free the memory.
         */
        ~Wider()
        {
            ;
        }
        
        //! This method set the radius.
        /**	The radius is between 0 and infinity. At 0, the source is in the center of the ambisonic circle and at 1, the source is at the limit of the ambisonic circle. Over 1, the source get away the ambisonic circle.
         @param     radius   The radius.
         @see       setAzimuth()
         */
        inline void setWidening(const T radius) noexcept
        {
            m_widening  = Math<T>::clip(radius, (T)0., (T)1.);
            m_factor    = (1. - m_widening) * HOA_PI;
            m_gain      = (sin(m_factor - HOA_PI2) + 1.) * 0.5;
        }
        
        //! Get the azimuth.
        /** The method returns the azimuth between 0 and 2π.
         @return     The azimuth.
         */
        inline T getWidening() const noexcept
        {
            return m_widening;
        }
        
        //! This method performs the encoding.
        /**	You should use this method for in-place or not-in-place processing and performs the encoding sample by sample. The outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         // cos(x + b) = cos(x) * cos(b) - sin(x) * sin(b)
         // sin(x + b) = cos(x) * sin(b) + sin(x) * cos(b)
         @param     input	The input sample.
         @param     outputs The output array.
         */
        inline void process(const T* inputs, T* outputs) const noexcept
        {
            const T gain1   = (m_gain * Harmonic<Hoa2d, T>::Processor::getDecompositionOrder());
            const T factor1 = (cos(Math<T>::clip(m_factor, 0., HOA_PI)) + 1.) * 0.5 * ((gain1 - m_gain)) + 1;
            
            (*outputs++) = (*inputs) * (gain1) + 1.;            // Hamonic [0,0]
            (*outputs++) = (*inputs) * factor1;                 // Hamonic [1,-1]
            (*outputs++) = (*inputs) * factor1;                 // Hamonic [1,1]
            for(ulong i = 2; i <= Harmonic<Hoa2d, T>::Processor::getDecompositionOrder(); i++)
            {
                const T gain    = (m_gain * (Harmonic<Hoa2d, T>::Processor::getDecompositionOrder() - i)) + 1.;
                const T factor  = (cos(Math<T>::clip(m_factor * i, 0., HOA_PI)) + 1.) * 0.5 ;
                
                (*outputs++)    = (*inputs) * factor * gain;    // Hamonic [i,-i]
                (*outputs++)    = (*inputs) * factor * gain;    // Hamonic [i,i]
            }
        }

    };
}

#endif



