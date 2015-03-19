/*
// Copyright (c) 2012-2015 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_WIDER_LIGHT
#define DEF_HOA_WIDER_LIGHT

#include "Harmonics.hpp"

namespace hoa
{
    //! The wider class wides the diffusion of the sound field in the harmonics domain.
    /** The wider simulates fractional orders to diffuse a sound field from omni directional when the widening factor is \f$0\f$ to the most directional when the widening factor is \f$1\f$.
     */
    template <Dimension D, typename T> class Wider : public Harmonic<D, T>::Processor
    {
    public:
        
        //! The wider constructor.
        /**	The wider constructor allocates and initialize the member values. The order must be at least 1.
         @param     order	The order.
         */
        Wider(const ulong order) noexcept = 0;
        
        //! The wider destructor.
        /**	The wider destructor free the memory.
         */
        virtual ~Wider() noexcept = 0;
        
        //! This method set the widening value.
        /**	The the widening value is between \f$0\f$ and \f$1\f$. At \f$0\f$, the sound field is omni directional and at \f$1\f$ the sound field is instact.
         @param     radius   The radius.
         @see       setAzimuth()
         */
        virtual void setWidening(const T radius) noexcept = 0;
        
        //! Get the the widening value.
        /** The method returns the the widening value.
         @return     The widening value.
         */
        virtual T getWidening() const noexcept = 0;
        
        //! This method perform the widening.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs	The input array.
         @param     outputs The output array.
         */
        virtual void process(const T* inputs, T* outputs) const noexcept = 0;

    };
    
    template <typename T> class Wider<Hoa2d, T> : public Processor< Harmonic<Hoa2d, T> >
    {
    private:
        T   m_widening;
        T   m_gain;
        T   m_factor;
    public:
        
        //! The wider constructor.
        /**	The wider constructor allocates and initialize the member values. The order must be at least 1.
         @param     order	The order.
         */
        Wider(const ulong order) noexcept : Processor< Harmonic<Hoa2d, T> >(order)
        {
            setWidening(1.);
        }
        
        //! The wider destructor.
        /**	The wider destructor free the memory.
         */
        ~Wider() noexcept
        {
            ;
        }
        
        //! This method set the widening value.
        /**	The the widening value is between \f$0\f$ and \f$1\f$. At \f$0\f$, the sound field is omni directional and at \f$1\f$ the sound field is instact.
         @param     radius   The radius.
         @see       setAzimuth()
         */
        inline void setWidening(const T widening) noexcept
        {
            m_widening  = Math<T>::clip(widening, (T)0., (T)1.);
            m_factor    = (1. - m_widening) * HOA_PI;
            m_gain      = (sin(m_factor - HOA_PI2) + 1.) * 0.5;
        }
        
        //! Get the the widening value.
        /** The method returns the the widening value.
         @return     The widening value.
         */
        inline T getWidening() const noexcept
        {
            return m_widening;
        }
        
        //! This method perform the widening.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs	The input array.
         @param     outputs The output array.
         */
        inline void process(const T* inputs, T* outputs) const noexcept
        {
            T gain   = (m_gain * Processor< Harmonic<Hoa2d, T> >::getDecompositionOrder());
            T factor = (cos(Math<T>::clip(m_factor, 0., HOA_PI)) + 1.) * 0.5 * ((gain - m_gain) + 1.);
            
            (*outputs++) = (*inputs++) * (gain + 1.);            // Hamonic [0,0]
            (*outputs++) = (*inputs++) * factor;                 // Hamonic [1,-1]
            (*outputs++) = (*inputs++) * factor;                 // Hamonic [1,1]
            for(ulong i = 2; i <= Processor< Harmonic<Hoa2d, T> >::getDecompositionOrder(); i++)
            {
                gain    = (m_gain * (Processor< Harmonic<Hoa2d, T> >::getDecompositionOrder() - i) + 1.);
                factor  = (cos(Math<T>::clip(m_factor * i, 0., HOA_PI)) + 1.) * 0.5 ;
                
                (*outputs++)    = (*inputs++) * factor * gain;    // Hamonic [i,-i]
                (*outputs++)    = (*inputs++) * factor * gain;    // Hamonic [i,i]
            }
        }
    };
    
    template <typename T> class Wider<Hoa3d, T> : public Processor< Harmonic<Hoa3d, T> >
    {
    private:
        T   m_widening;
        T   m_gain;
        T   m_factor;
    public:
        
        //! The wider constructor.
        /**	The wider constructor allocates and initialize the member values. The order must be at least 1.
         @param     order	The order.
         */
        Wider(const ulong order) noexcept : Processor< Harmonic<Hoa3d, T> >(order)
        {
            setWidening(1.);
        }
        
        //! The wider destructor.
        /**	The wider destructor free the memory.
         */
        ~Wider() noexcept
        {
            ;
        }
        
        //! This method set the widening value.
        /**	The the widening value is between \f$0\f$ and \f$1\f$. At \f$0\f$, the sound field is omni directional and at \f$1\f$ the sound field is instact.
         @param     radius   The radius.
         @see       setAzimuth()
         */
        inline void setWidening(const T radius) noexcept
        {
            m_widening  = Math<T>::clip(radius, (T)0., (T)1.);
            m_factor    = (1. - m_widening) * HOA_PI;
            m_gain      = (sin(m_factor - HOA_PI2) + 1.) * 0.5;
        }
        
        //! Get the the widening value.
        /** The method returns the the widening value.
         @return     The widening value.
         */
        inline T getWidening() const noexcept
        {
            return m_widening;
        }
        
        //! This method perform the widening.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs	The input array.
         @param     outputs The output array.
         */
        inline void process(const T* inputs, T* outputs) const noexcept
        {
            T gain   = (m_gain * Processor< Harmonic<Hoa3d, T> >::getDecompositionOrder());
            T factor = (cos(Math<T>::clip(m_factor, 0., HOA_PI)) + 1.) * 0.5 * ((gain - m_gain) + 1.);
            
            (*outputs++) = (*inputs++) * (gain + 1.);            // Hamonic [0,0]
            (*outputs++) = (*inputs++) * factor;                 // Hamonic [1,-1]
            (*outputs++) = (*inputs++) * factor;                 // Hamonic [1,0]
            (*outputs++) = (*inputs++) * factor;                 // Hamonic [1,1]
            for(ulong i = 2; i <= Processor< Harmonic<Hoa3d, T> >::getDecompositionOrder(); i++)
            {
                gain    = (m_gain * (Processor< Harmonic<Hoa3d, T> >::getDecompositionOrder() - i) + 1.);
                factor  = (cos(Math<T>::clip(m_factor * i, 0., HOA_PI)) + 1.) * 0.5 ;
                
                for(ulong j = 0; j < 2 * i + 1; j++)
                {
                    (*outputs++)    = (*inputs++) * factor * gain;    // Hamonic [i, ~j]
                }
            }
        }
    };
}

#endif



