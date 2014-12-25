/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_2D_MAP
#define DEF_HOA_2D_MAP

#include "Ambisonic_2D.h"
#include "Encoder_2D.h"
#include "Wider_2D.h"

namespace Hoa2D
{
    //! The ambisonic encoder with distance compensation.
    /** The map is a Encoder with distance compensation. It uses intances of the Wider class to decrease the directionnality of sources by simulating fractionnal orders when the sources are inside the ambisonic circle and a simple diminution of the gain when the sources get away from the ambisonic circle.
     @see Encoder
     @see Wider
     */
    template <typename Type> class MapUnique : public Ambisonic
    {
    private:
        Type    m_azimuth;
        Type    m_cosx;
        Type    m_sinx;
        Type    m_factor;
        Type    m_gain;
        Type    m_radius;
        bool    m_muted;
        
    public:
        
        //! The unique map constructor.
        /**	The map constructor allocates and initialize the member values and classes depending of a decomposition order. The order must be at least 1.
         @param     order            The order.
         */
        MapUnique(unsigned long order) noexcept : Ambisonic(order)
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
        inline void setAzimuth(const Type azimuth) noexcept
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
        inline void setRadius(const Type radius) noexcept
        {
            m_radius = max(radius, (Type)0.);
            if(m_radius < 1.)
            {
                m_factor = (1. - clip(radius, 0., 1.)) * HOA_PI;;
                m_gain   = (std::sin(clip((m_factor - 0.5f) * HOA_PI, -HOA_PI2, HOA_PI2)) + 1.) * 0.5;
            }
            else
            {
                m_factor = HOA_PI;
                m_gain   = (std::sin(clip((m_factor - 0.5f) * HOA_PI, -HOA_PI2, HOA_PI2)) + 1.) * 0.5 / radius;
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
        inline Type getAzimuth() const noexcept
        {
            return m_azimuth;
        }
        
        //! This method retrieve the radius of a source.
        /** Retrieve the radius of a source.
         @param     index	The index of the source.
         @return The radius of the source if the source exists, otherwise the function generates an error.
         */
        inline Type getRadius() const noexcept
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
        inline void process(const Type input, Type* outputs) const noexcept
        {
            if(!m_muted)
            {
                Type cos_x = m_cosx;
                Type sin_x = m_sinx;
                Type tcos_x = cos_x;
                outputs[0] = input * (m_gain * m_order_of_decomposition + 1.f);
                for(unsigned long i = 2, j = 1; i < m_number_of_harmonics; i += 2, j++)
                {
                    const Type factor  = (std::cos(min(m_factor * j, (Type)HOA_PI)) + 1.f) * 0.5f * (m_gain * (m_order_of_decomposition - j) + 1.f);
                    outputs[i-1]    = input * sin_x * factor;
                    outputs[i]      = input * cos_x * factor;
                    cos_x = tcos_x * m_cosx - sin_x * m_sinx; // cos(x + b) = cos(x) * cos(b) - sin(x) * sin(b)
                    sin_x = tcos_x * m_sinx + sin_x * m_cosx; // sin(x + b) = cos(x) * sin(b) + sin(x) * cos(b)
                    tcos_x = cos_x;
                }
            }
            else
            {
                for(unsigned long i = 0; i < m_number_of_harmonics; i++)
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
        inline void processAdd(const Type input, Type* outputs) noexcept
        {
            if(!m_muted)
            {
                Type cos_x = m_cosx;
                Type sin_x = m_sinx;
                Type tcos_x = cos_x;
                outputs[0] = input * (m_gain * m_order_of_decomposition + 1.f);
                for(unsigned long i = 2, j = 1; i < m_number_of_harmonics; i += 2, j++)
                {
                    const Type factor  = (std::cos(min(m_factor * j, (Type)HOA_PI)) + 1.f) * 0.5f * (m_gain * (m_order_of_decomposition - j) + 1.f);
                    outputs[i-1]    += input * sin_x * factor;
                    outputs[i]      += input * cos_x * factor;
                    cos_x = tcos_x * m_cosx - sin_x * m_sinx; // cos(x + b) = cos(x) * cos(b) - sin(x) * sin(b)
                    sin_x = tcos_x * m_sinx + sin_x * m_cosx; // sin(x + b) = cos(x) * sin(b) + sin(x) * cos(b)
                    tcos_x = cos_x;
                }
            }
        }
    };
    
    //! The ambisonic multi-encoder with distance compensation.
    /** The map is a multi Encoder with distance compensation. It uses intances of the Wider class to decrease the directionnality of sources by simulating fractionnal orders when the sources are inside the ambisonic circle and a simple diminution of the gain when the sources get away from the ambisonic circle.
     
        @see Encoder
     */
    template <typename Type> class Map : public Ambisonic
    {
    private:
        const unsigned long m_number_of_sources;
        vector<MapUnique<Type>*>  m_maps;
    public:
        
        //! The map constructor.
        /**	The map constructor allocates and initialize the member values and classes depending of a decomposition order and the number of sources. The order and the number of sources must be at least 1.
         
            @param     order            The order.
            @param     numberOfSources	The number of sources.
         */
        Map(unsigned long order, unsigned long numberOfSources) noexcept : Ambisonic(order),
        m_number_of_sources(numberOfSources)
        {
            for(unsigned long i = 0; i < m_number_of_sources; i++)
            {
                m_maps.push_back(new MapUnique<Type>(order));
            }
        }
        
        //! The map destructor.
        /**	The map destructor free the memory and deallocate the member classes.
         */
        ~Map()
        {
            for(unsigned long i = 0; i < m_number_of_sources; i++)
            {
                delete m_maps[i];
            }
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
        inline void setAzimuth(const unsigned long index, const Type azimuth) noexcept
        {
            m_maps[index]->setAzimuth(azimuth);
        }
        
        //! This method set the radius of a source.
        /**	The radius is between 0 and infinity. At 0, the source is in the center of the ambisonic circle and at 1, the source is at the limit of the ambisonic circle. Over 1, the source get away the ambisonic circle. The index must be between 0 and the number of sources - 1.
         
            @param     index	The index of the source.
            @param     radius   The radius.
            @see       setAzimuth()
         */
        inline void setRadius(const unsigned long index, const Type radius) noexcept
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
        inline Type getAzimuth(const unsigned long index) const noexcept
        {
            return m_maps[index]->getAzimuth();
        }
		
        //! This method retrieve the radius of a source.
        /** Retrieve the radius of a source.
         
            @param     index	The index of the source.
            @return The radius of the source if the source exists, otherwise the function generates an error.
         */
        inline Type getRadius(const unsigned long index) const noexcept
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
        inline void process(const Type* inputs, Type* outputs) noexcept
        {
            m_maps[0]->process(inputs[0], outputs);
            for(unsigned long i = 1; i < m_number_of_sources; i++)
            {
                m_maps[i]->processAdd(inputs[i], outputs);
            }
        }
    };
}

#endif
