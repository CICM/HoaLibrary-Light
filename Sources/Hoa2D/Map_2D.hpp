/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_2D_MAP
#define DEF_HOA_2D_MAP

#include "Wider_2D.h"

namespace hoa
{
    //! The ambisonic encoder with distance compensation.
    /** The map is a Encoder with distance compensation. It uses intances of the Wider class to decrease the directionnality of sources by simulating fractionnal orders when the sources are inside the ambisonic circle and a simple diminution of the gain when the sources get away from the ambisonic circle.
     @see Encoder
     @see Wider
     */
    template <typename T> class MapUnique : public Ambisonic2D<T>
    {
    private:
        T    m_azimuth;
        T    m_cosx;
        T    m_sinx;
        T    m_factor;
        T    m_gain;
        T    m_radius;
        bool m_muted;
        
    public:
        
        //! The unique map constructor.
        /**	The map constructor allocates and initialize the member values and classes depending of a decomposition order. The order must be at least 1.
         @param     order            The order.
         */
        MapUnique(unsigned long order) noexcept : Ambisonic2D<T>(order)
        {
            setMute(0);
            setAzimuth(0.);
            setRadius(1.);
        }
        
        //! The map destructor.
        /**	The map destructor free the memory.
         */
        ~MapUnique()
        {
            ;
        }
        
        //! This method set the angle of azimuth.
        /**	The angle of azimuth in radian and you should prefer to use it between 0 and 2 π to avoid recursive wrapping of the value. The direction of rotation is counterclockwise. The 0 radian is π/2 phase shifted relative to a mathematical representation of a circle, then the 0 radian is at the "front" of the soundfield.
         @param     azimuth	The azimuth.
         */
        inline void setAzimuth(const T azimuth) noexcept
        {
            m_azimuth = wrap_twopi(azimuth);
            m_cosx    = std::cos(m_azimuth);
            m_sinx    = std::sin(m_azimuth);
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
                m_factor = (1. - clip(radius, 0., 1.)) * HOA_PI;
                m_gain   = (sin(m_factor - HOA_PI2) + 1.) * 0.5;
            }
            else
            {
                m_factor = 0;
                m_gain   = 0;
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
        
        //! Get the azimuth angle
        /** The method returns the last angle of encoding between 0 and 2π.
         @return     The azimuth.
         */
        inline T getAzimuth() const noexcept
        {
            return m_azimuth;
        }
        
        //! This method retrieve the radius of a source.
        /** Retrieve the radius of a source.
         @param     index	The index of the source.
         @return The radius of the source if the source exists, otherwise the function generates an error.
         */
        inline T getRadius() const noexcept
        {
            return m_radius;
        }
        
		//! This method retrieve the mute or unmute state of a source.
        /**	Get the Mute state of a source.
         @return    The mute state of the source.
         */
        inline bool getMute() const noexcept
        {
            return m_muted;
        }
        
        //! This method performs the encoding with distance compensation with single precision.
        /**	You should use this method for in-place or not-in-place processing and performs the encoding with distance compensation sample by sample. The input contains the samples of the source. The outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     input  The input.
         @param     outputs The outputs array.
         */
        inline void process(const T* inputs, T* outputs) const noexcept
        {
            if(!m_muted)
            {
                T cos_x = m_cosx;
                T sin_x = m_sinx;
                T tcos_x = cos_x;
                const T gain1   = (m_gain * Ambisonic<T>::m_order_of_decomposition + 1.f);
                const T factor1 = (cos(clip(m_factor, 0., HOA_PI)) + 1.f) *
                0.5f * (m_gain * (Ambisonic<T>::m_order_of_decomposition - 1) + 1.f);
                (*outputs++) = (*inputs) * gain1;                           // Hamonic [0,0]
                (*outputs++) = (*inputs) * sin_x * factor1;                 // Hamonic [1,-1]
                (*outputs++) = (*inputs) * cos_x * factor1;                 // Hamonic [1,1]
                for(unsigned long i = 2; i <= Ambisonic<T>::m_order_of_decomposition; i++)
                {
                    const T factor  = (cos(clip(m_factor * i, 0., HOA_PI)) + 1.f) *
                    0.5f * (m_gain * (Ambisonic<T>::m_order_of_decomposition - i) + 1.f);
                    cos_x   = tcos_x * m_cosx - sin_x * m_sinx;
                    sin_x   = tcos_x * m_sinx + sin_x * m_cosx;
                    tcos_x  = cos_x;
                    (*outputs++)    = (*inputs) * sin_x * factor;
                    (*outputs++)    = (*inputs) * cos_x * factor;
                }
            }
            else
            {
                for(unsigned long i = 0; i < Ambisonic<T>::m_number_of_harmonics; i++)
                {
                    outputs[i] = 0.f;
                }
            }
        }
        
        //! This method performs the encoding with distance compensation with single precision.
        /**	You should use this method for in-place or not-in-place processing and performs the encoding with distance compensation sample by sample. The input contains the samples of the source. The outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     input  The input.
         @param     outputs The outputs array.
         */
        inline void processAdd(const T* inputs, T* outputs) const noexcept
        {
            if(!m_muted)
            {
                T cos_x = m_cosx;
                T sin_x = m_sinx;
                T tcos_x = cos_x;
                const T gain1   = (m_gain * Ambisonic<T>::m_order_of_decomposition + 1.f);
                const T factor1 = (cos(clip(m_factor, 0., HOA_PI)) + 1.f) *
                0.5f * (m_gain * (Ambisonic<T>::m_order_of_decomposition - 1) + 1.f);
                (*outputs++) += (*inputs) * gain1;                           // Hamonic [0,0]
                (*outputs++) += (*inputs) * sin_x * factor1;                 // Hamonic [1,-1]
                (*outputs++) += (*inputs) * cos_x * factor1;                 // Hamonic [1,1]
                for(unsigned long i = 2; i <= Ambisonic<T>::m_order_of_decomposition; i++)
                {
                    const T factor  = (cos(clip(m_factor * i, 0., HOA_PI)) + 1.f) *
                    0.5f * (m_gain * (Ambisonic<T>::m_order_of_decomposition - i) + 1.f);
                    cos_x   = tcos_x * m_cosx - sin_x * m_sinx;
                    sin_x   = tcos_x * m_sinx + sin_x * m_cosx;
                    tcos_x  = cos_x;
                    (*outputs++)    += (*inputs) * sin_x * factor;
                    (*outputs++)    += (*inputs) * cos_x * factor;
                }
            }
        }
    };
    
    //! The ambisonic multi-encoder with distance compensation.
    /** The map is a multi Encoder with distance compensation. It uses intances of the Wider class to decrease the directionnality of sources by simulating fractionnal orders when the sources are inside the ambisonic circle and a simple diminution of the gain when the sources get away from the ambisonic circle.
     
        @see Encoder
     */
    template <typename T> class Map : public Ambisonic2D<T>
    {
    private:
        const unsigned long m_number_of_sources;
        vector<unique_ptr<MapUnique<T>>>  m_maps;
    public:
        
        //! The map constructor.
        /**	The map constructor allocates and initialize the member values and classes depending of a decomposition order and the number of sources. The order and the number of sources must be at least 1.
         
            @param     order            The order.
            @param     numberOfSources	The number of sources.
         */
        Map(unsigned long order, unsigned long numberOfSources) noexcept : Ambisonic2D<T>(order),
        m_number_of_sources(numberOfSources)
        {
            for(unsigned long i = 0; i < m_number_of_sources; i++)
            {
                m_maps.push_back(unique_ptr<MapUnique<T>>(new MapUnique<T>(order)));
            }
        }
        
        //! The map destructor.
        /**	The map destructor free the memory and deallocate the member classes.
         */
        ~Map()
        {
            m_maps.clear();
        }
        
        //! This method retrieve the number of sources.
        /** Retrieve the number of sources.
         
            @return The number of sources.
         */
        inline unsigned long getNumberOfSources() const noexcept
        {
            return m_number_of_sources;
        }
        
        //! This method set the angle of azimuth of a source.
        /**	The angle of azimuth in radian and you should prefer to use it between 0 and 2 π to avoid recursive wrapping of the value. The direction of rotation is counterclockwise. The 0 radian is π/2 phase shifted relative to a mathematical representation of a circle, then the 0 radian is at the "front" of the soundfield. The index must be between 0 and the number of sources - 1.
         
            @param     index	The index of the source.
            @param     azimuth	The azimuth.
            @see       setRadius()
         */
        inline void setAzimuth(const unsigned long index, const T azimuth) noexcept
        {
            m_maps[index]->setAzimuth(azimuth);
        }
        
        //! This method set the radius of a source.
        /**	The radius is between 0 and infinity. At 0, the source is in the center of the ambisonic circle and at 1, the source is at the limit of the ambisonic circle. Over 1, the source get away the ambisonic circle. The index must be between 0 and the number of sources - 1.
         
            @param     index	The index of the source.
            @param     radius   The radius.
            @see       setAzimuth()
         */
        inline void setRadius(const unsigned long index, const T radius) noexcept
        {
            m_maps[index]->setRadius(radius);
        }
		
		//! This method mute or unmute a source.
        /**	Mute or unmute a source with a boolean value. The index must be between 0 and the number of sources - 1.
         
            @param     index	The index of the source.
            @param     muted	The mute state.
         */
        inline void setMute(const unsigned long index, const bool muted) noexcept
        {
            m_maps[index]->setMute(muted);
        }
        
        //! This method retrieve the azimuth of a source.
        /** Retrieve the azimuth of a source.
         
            @param     index	The index of the source.
            @return The azimuth of the source if the source exists, otherwise the function generates an error.
         */
        inline T getAzimuth(const unsigned long index) const noexcept
        {
            return m_maps[index]->getAzimuth();
        }
		
        //! This method retrieve the radius of a source.
        /** Retrieve the radius of a source.
         
            @param     index	The index of the source.
            @return The radius of the source if the source exists, otherwise the function generates an error.
         */
        inline T getRadius(const unsigned long index) const noexcept
        {
            return m_maps[index]->getRadius();
        }
        
		//! This method retrieve the mute or unmute state of a source.
        /**	Get the Mute state of a source.
         
            @param     index	The index of the source.
            @return    The mute state of the source if the source exists, otherwise the function generates an error.
            @see       setMute()
         */
        inline bool getMute(const unsigned long index) const noexcept
        {
            return m_maps[index]->getMute();
        }
		
        
        //! This method performs the encoding with distance compensation.
        /**	You should use this method for in-place or not-in-place processing and performs the encoding with distance compensation sample by sample. The inputs array contains the samples of the sources and the minimum size sould be the number of sources. The outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs  The inputs array.
         @param     outputs The outputs array.
         */
        inline void process(const T* inputs, T* outputs) const noexcept
        {
            m_maps[0]->process(inputs, outputs);
            for(unsigned long i = 1; i < m_number_of_sources; i++)
            {
                m_maps[i]->processAdd(inputs++, outputs);
            }
        }
    };
}

#endif
