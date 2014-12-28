/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_ENCODER
#define DEF_HOA_ENCODER

#include "Harmonics.hpp"

namespace hoa
{
    //! The ambisonic encoder.
    /** The encoder should be used to encode a source in the spherical harmonics domain depending of an order of decomposition. It allows to control the azimuth of the source.
     */
    template <Dimension D, typename T> class Encoder : public Harmonic<D>::Processor
    {
    private:
        T    m_azimuth;
        T    m_cosx;
        T    m_sinx;
        bool m_muted;
    public:
        
        //! The encoder constructor.
        /**	The encoder constructor allocates and initialize the member values to computes circular harmonics coefficients for the encoding. The order must be at least 1.
            @param     order	The order.
         */
        Encoder(const ulong order) :
        Harmonic<D>::Processor(order)
        {
            setMute(false);
            setAzimuth(0.);
        }
        
        //! The encoder destructor.
        /**	The encoder destructor free the memory.
         */
        ~Encoder()
        {
            ;
        }
        
        //! This method set the azimuth.
        /**	The azimuth in radian and you should prefer to use it between 0 and 2π to avoid recursive wrapping of the value. The direction of rotation is counterclockwise. The 0 radian is π/2 phase shifted relative to a mathematical representation of a circle, then the 0 radian is at the "front" of the soundfield.
            @param     azimuth	The azimuth.
         */
        inline void setAzimuth(const T azimuth) noexcept
        {
            m_azimuth = azimuth;
            m_cosx    = std::cos(m_azimuth);
            m_sinx    = std::sin(m_azimuth);
        }
        
        //! Get the azimuth.
        /** The method returns the azimuth between 0 and 2π.
            @return     The azimuth.
         */
        inline T getAzimuth() const noexcept
        {
            return wrap_twopi(m_azimuth);
        }
        
        //! This method mute or unmute.
        /**	Mute or unmute.
         @param     muted	The mute state.
         */
        inline void setMute(const bool muted) noexcept
        {
            m_muted = muted;
        }
        
        //! This method retrieve the mute or unmute state of a source.
        /**	Get the Mute state of a source.
         @return    The mute state of the source.
         */
        inline bool getMute() const noexcept
        {
            return m_muted;
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
            if(!m_muted)
            {
                T cos_x = m_cosx, tcos_x = m_cosx, sin_x = m_sinx;
                (*outputs++)    = (*inputs);                         // Hamonic [0,0]
                (*outputs++)    = (*inputs) * sin_x;                 // Hamonic [1,-1]
                (*outputs++)    = (*inputs) * cos_x;                 // Hamonic [1,1]
                for(ulong i = 2; i <= Harmonic<D>::Processor::getDecompositionOrder(); i++)
                {
                    cos_x   = tcos_x * m_cosx - sin_x * m_sinx;
                    sin_x   = tcos_x * m_sinx + sin_x * m_cosx;
                    tcos_x  = cos_x;
                    (*outputs++)    = (*inputs) * sin_x;            // Hamonic [l,-l]
                    (*outputs++)    = (*inputs) * cos_x;            // Hamonic [l,l]
                }
            }
            else
            {
                for(ulong i = 0; i < Harmonic<D>::Processor::getNumberOfHarmonics(); i++)
                {
                    (*outputs++) = 0.;
                }
            }
        }
        
        //! This method performs the encoding.
        /**	You should use this method for in-place or not-in-place processing and performs the encoding sample by sample. The outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         // cos(x + b) = cos(x) * cos(b) - sin(x) * sin(b)
         // sin(x + b) = cos(x) * sin(b) + sin(x) * cos(b)
         @param     input	The input sample.
         @param     outputs The output array.
         */
        inline void processAdd(const T* inputs, T* outputs) const noexcept
        {
            if(!m_muted)
            {
                T cos_x = m_cosx, tcos_x = m_cosx, sin_x = m_sinx;
                (*outputs++)    += (*inputs);                         // Hamonic [0,0]
                (*outputs++)    += (*inputs) * sin_x;                 // Hamonic [1,-1]
                (*outputs++)    += (*inputs) * cos_x;                 // Hamonic [1,1]
                for(ulong i = 2; i <= Harmonic<D>::Processor::getDecompositionOrder(); i++)
                {
                    cos_x   += tcos_x * m_cosx - sin_x * m_sinx;
                    sin_x   += tcos_x * m_sinx + sin_x * m_cosx;
                    tcos_x  += cos_x;
                    (*outputs++)    += (*inputs) * sin_x;            // Hamonic [i,-i]
                    (*outputs++)    += (*inputs) * cos_x;            // Hamonic [i,i]
                }
            }
        }
    };
    
    //! The ambisonic encoder with distance compensation.
    /** The encoder with distance compensation should be used to encode a source in the spherical harmonics domain depending of an order of decomposition. It allows to control the azimuth and the radius of the source.
     */
    template <Dimension D, typename T> class EncoderDC : public Harmonic<D>::Processor
    {
    private:
        T   m_azimuth;
        T   m_cosx;
        T   m_sinx;
        T   m_factor;
        T   m_gain;
        T   m_radius;
        T   m_distance;
        bool m_muted;
    public:
        
        //! The encoder constructor.
        /**	The encoder constructor allocates and initialize the member values to computes circular harmonics coefficients for the encoding. The order must be at least 1.
         @param     order	The order.
         */
        EncoderDC(const ulong order) :
        Harmonic<D>::Processor(order)
        {
            setAzimuth(0.);
            setRadius(1.);
            setMute(false);
        }
        
        //! The encoder destructor.
        /**	The encoder destructor free the memory.
         */
        ~EncoderDC()
        {
            ;
        }
        
        //! This method set the azimuth.
        /**	The azimuth in radian and you should prefer to use it between 0 and 2π to avoid recursive wrapping of the value. The direction of rotation is counterclockwise. The 0 radian is π/2 phase shifted relative to a mathematical representation of a circle, then the 0 radian is at the "front" of the soundfield.
         @param     azimuth	The azimuth.
         */
        inline void setAzimuth(const T azimuth) noexcept
        {
            m_azimuth = azimuth;
            m_cosx    = std::cos(m_azimuth);
            m_sinx    = std::sin(m_azimuth);
        }
        
        //! Get the azimuth.
        /** The method returns the azimuth between 0 and 2π.
         @return     The azimuth.
         */
        inline T getAzimuth() const noexcept
        {
            return wrap_twopi(m_azimuth);
        }
        
        //! This method set the radius.
        /**	The radius is between 0 and infinity. At 0, the source is in the center of the ambisonic circle and at 1, the source is at the limit of the ambisonic circle. Over 1, the source get away the ambisonic circle.
         @param     radius   The radius.
         @see       setAzimuth()
         */
        inline void setRadius(const T radius) noexcept
        {
            m_radius = max(radius, (T)0.);
            if(m_radius < 1.)
            {
                m_factor    = (1. - clip(radius, 0., 1.)) * HOA_PI;
                m_gain      = (sin(m_factor - HOA_PI2) + 1.) * 0.5;
                m_distance  = 1.;
            }
            else
            {
                m_factor    = 0;
                m_gain      = 0;
                m_distance  = 1. / std::max(radius, (T)1.);
            }
        }
        
        //! This method mute or unmute.
        /**	Mute or unmute.
         @param     muted	The mute state.
         */
        inline void setMute(const bool muted) noexcept
        {
            m_muted = muted;
        }
        
        //! This method retrieve the mute or unmute state of a source.
        /**	Get the Mute state of a source.
         @return    The mute state of the source.
         */
        inline bool getMute() const noexcept
        {
            return m_muted;
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
            if(!m_muted)
            {
                T cos_x = m_cosx;
                T sin_x = m_sinx;
                T tcos_x = cos_x;
                const T gain1   = (m_gain * Harmonic<D>::Processor::getDecompositionOrder());
                const T factor1 = (cos(clip(m_factor, 0., HOA_PI)) + 1.) * 0.5 * ((gain1 - m_gain) + m_distance);
                
                (*outputs++) = (*inputs) * (gain1 + m_distance);            // Hamonic [0,0]
                (*outputs++) = (*inputs) * sin_x * factor1;                 // Hamonic [1,-1]
                (*outputs++) = (*inputs) * cos_x * factor1;                 // Hamonic [1,1]
                for(ulong i = 2; i <= Harmonic<D>::Processor::getDecompositionOrder(); i++)
                {
                    const T gain    = (m_gain * (Harmonic<D>::Processor::getDecompositionOrder() - i) + m_distance);
                    const T factor  = (cos(clip(m_factor * i, 0., HOA_PI)) + 1.) * 0.5 ;
                    
                    cos_x   = tcos_x * m_cosx - sin_x * m_sinx;
                    sin_x   = tcos_x * m_sinx + sin_x * m_cosx;
                    tcos_x  = cos_x;
                    
                    (*outputs++)    = (*inputs) * sin_x * factor * gain;    // Hamonic [i,-i]
                    (*outputs++)    = (*inputs) * cos_x * factor * gain;    // Hamonic [i,i]
                }
            }
            else
            {
                for(ulong i = 0; i < Harmonic<D>::Processor::getNumberOfHarmonics(); i++)
                {
                    (*outputs++) = 0.;
                }
            }
        }
        
        //! This method performs the encoding.
        /**	You should use this method for in-place or not-in-place processing and performs the encoding sample by sample. The outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         // cos(x + b) = cos(x) * cos(b) - sin(x) * sin(b)
         // sin(x + b) = cos(x) * sin(b) + sin(x) * cos(b)
         @param     input	The input sample.
         @param     outputs The output array.
         */
        inline void processAdd(const T* inputs, T* outputs) const noexcept
        {
            if(!m_muted)
            {
                T cos_x = m_cosx;
                T sin_x = m_sinx;
                T tcos_x = cos_x;
                const T gain1   = (m_gain * Harmonic<D>::Processor::getDecompositionOrder());
                const T factor1 = (cos(clip(m_factor, 0., HOA_PI)) + 1.) * 0.5 * ((gain1 - m_gain) + m_distance);
                
                (*outputs++) += (*inputs) * (gain1 + m_distance);            // Hamonic [0,0]
                (*outputs++) += (*inputs) * sin_x * factor1;                 // Hamonic [1,-1]
                (*outputs++) += (*inputs) * cos_x * factor1;                 // Hamonic [1,1]
                for(ulong i = 2; i <= Harmonic<D>::Processor::getDecompositionOrder(); i++)
                {
                    const T gain    = (m_gain * (Harmonic<D>::Processor::getDecompositionOrder() - i) + m_distance);
                    const T factor  = (cos(clip(m_factor * i, 0., HOA_PI)) + 1.) * 0.5 ;
                    
                    cos_x   = tcos_x * m_cosx - sin_x * m_sinx;
                    sin_x   = tcos_x * m_sinx + sin_x * m_cosx;
                    tcos_x  = cos_x;
                    
                    (*outputs++)    += (*inputs) * sin_x * factor * gain;    // Hamonic [i,-i]
                    (*outputs++)    += (*inputs) * cos_x * factor * gain;    // Hamonic [i,i]
                }
            }
        }
    };
    
    //! The ambisonic multi-encoder with distance compensation.
    /** The map is a multi encoder with distance compensation.
     */
    template <Dimension D, typename T> class EncoderMulti : public Harmonic<D>::Processor
    {
    private:
        const ulong                         m_number_of_sources;
        vector<unique_ptr<EncoderDC<D, T>>> m_encoders;
    public:
        
        //! The map constructor.
        /**	The map constructor allocates and initialize the member values and classes depending of a decomposition order and the number of sources. The order and the number of sources must be at least 1.
         
         @param     order            The order.
         @param     numberOfSources	The number of sources.
         */
        EncoderMulti(ulong order, ulong numberOfSources) noexcept :
        Harmonic<D>::Processor(order),
        m_number_of_sources(numberOfSources)
        {
            for(ulong i = 0; i < m_number_of_sources; i++)
            {
                m_encoders.push_back(unique_ptr<EncoderDC<D, T>>(new EncoderDC<D, T>(order)));
            }
        }
        
        //! The map destructor.
        /**	The map destructor free the memory and deallocate the member classes.
         */
        ~EncoderMulti()
        {
            m_encoders.clear();
        }
        
        //! This method retrieve the number of sources.
        /** Retrieve the number of sources.
         
         @return The number of sources.
         */
        inline ulong getNumberOfSources() const noexcept
        {
            return m_number_of_sources;
        }
        
        //! This method set the angle of azimuth of a source.
        /**	The angle of azimuth in radian and you should prefer to use it between 0 and 2 π to avoid recursive wrapping of the value. The direction of rotation is counterclockwise. The 0 radian is π/2 phase shifted relative to a mathematical representation of a circle, then the 0 radian is at the "front" of the soundfield. The index must be between 0 and the number of sources - 1.
         
         @param     index	The index of the source.
         @param     azimuth	The azimuth.
         @see       setRadius()
         */
        inline void setAzimuth(const ulong index, const T azimuth) noexcept
        {
            m_encoders[index]->setAzimuth(azimuth);
        }
        
        //! This method set the radius of a source.
        /**	The radius is between 0 and infinity. At 0, the source is in the center of the ambisonic circle and at 1, the source is at the limit of the ambisonic circle. Over 1, the source get away the ambisonic circle. The index must be between 0 and the number of sources - 1.
         
         @param     index	The index of the source.
         @param     radius   The radius.
         @see       setAzimuth()
         */
        inline void setRadius(const ulong index, const T radius) noexcept
        {
            m_encoders[index]->setRadius(radius);
        }
        
        //! This method mute or unmute a source.
        /**	Mute or unmute a source with a boolean value. The index must be between 0 and the number of sources - 1.
         
         @param     index	The index of the source.
         @param     muted	The mute state.
         */
        inline void setMute(const ulong index, const bool muted) noexcept
        {
            m_encoders[index]->setMute(muted);
        }
        
        //! This method retrieve the azimuth of a source.
        /** Retrieve the azimuth of a source.
         
         @param     index	The index of the source.
         @return The azimuth of the source if the source exists, otherwise the function generates an error.
         */
        inline T getAzimuth(const ulong index) const noexcept
        {
            return m_encoders[index]->getAzimuth();
        }
        
        //! This method retrieve the radius of a source.
        /** Retrieve the radius of a source.
         
         @param     index	The index of the source.
         @return The radius of the source if the source exists, otherwise the function generates an error.
         */
        inline T getRadius(const ulong index) const noexcept
        {
            return m_encoders[index]->getRadius();
        }
        
        //! This method retrieve the mute or unmute state of a source.
        /**	Get the Mute state of a source.
         
         @param     index	The index of the source.
         @return    The mute state of the source if the source exists, otherwise the function generates an error.
         @see       setMute()
         */
        inline bool getMute(const ulong index) const noexcept
        {
            return m_encoders[index]->getMute();
        }
        
        
        //! This method performs the encoding with distance compensation.
        /**	You should use this method for in-place or not-in-place processing and performs the encoding with distance compensation sample by sample. The inputs array contains the samples of the sources and the minimum size sould be the number of sources. The outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs  The inputs array.
         @param     outputs The outputs array.
         */
        inline void process(const T* inputs, T* outputs) const noexcept
        {
            m_encoders[0]->process(inputs, outputs);
            for(ulong i = 1; i < m_number_of_sources; i++)
            {
                m_encoders[i]->processAdd(inputs++, outputs);
            }
        }
    };
}

#endif



