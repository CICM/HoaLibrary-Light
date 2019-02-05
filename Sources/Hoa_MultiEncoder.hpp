/*
// Copyright (c) 2012-2017 CICM - Universite Paris 8 - Labex Arts H2H.
// Authors :
// 2012: Pierre Guillot, Eliott Paris & Julien Colafrancesco.
// 2012-2015: Pierre Guillot & Eliott Paris.
// 2015: Pierre Guillot & Eliott Paris & Thomas Le Meur (Light version)
// 2016-2017: Pierre Guillot.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#pragma once

#include "Hoa_Encoder.hpp"

namespace hoa
{
    // ================================================================================ //
    // MULTI ENCODER //
    // ================================================================================ //
    
    //! @brief The class generates manages a set of encoders.
    //! @details The class is just a wrapper to manage a set encoders that can be muted.
    template <Dimension D, typename T>
    class MultiEncoder
    : public ProcessorHarmonics<D, T>
    {
    public:

        //! @brief Constructor.
        //! @param order The order.
        //! @param nsources The number of sources.
        MultiEncoder(size_t order, size_t nsources)
        : ProcessorHarmonics<D, T>(order)
        , m_encoders(nsources)
        {
            m_temp = new T[ProcessorHarmonics<D, T>::getNumberOfHarmonics()];
            for(size_t i = 0; i < nsources; ++i)
            {
                m_encoders[i].encoder = new Encoder<D, T>(order);
                m_encoders[i].encoder->setAzimuth(T(i) * (T(HOA_2PI) / T(nsources)));
            }
        }

        //! @brief Destructor.
        ~MultiEncoder()
        {
            for(size_t i = 0; i < getNumberOfSources(); ++i)
            {
                delete m_encoders[i].encoder;
            }
            
            delete [] m_temp;
        }

        //! @brief Returns the number of sources.
        inline size_t getNumberOfSources() const
        {
            return m_encoders.size();
        }
        
        //! @brief Sets the radius of a source.
        //! @param     index	The index of the source.
        //! @param     radius	The radius.
        inline void setRadius(size_t index, T radius)
        {
            m_encoders[index].setRadius(radius);
        }
        
        //! @brief Sets the widening factor of a source (wrapper for radius).
        //! @param     index	The index of the source.
        //! @param     widening	The widening value.
        inline void setWidening(size_t index, T widening)
        {
            m_encoders[index].setRadius((widening > T(1.)) ? T(1.) : widening);
        }

        //! @brief Sets the azimuth of a source.
        //! @param     index	The index of the source.
        //! @param     azimuth	The azimuth.
        inline void setAzimuth(size_t index, T azimuth)
        {
            m_encoders[index].setAzimuth(azimuth);
        }

        //! @brief Sets the elevation of a source.
        //! @param     index	The index of the source.
        //! @param     elevation	The elevation.
        inline void setElevation(size_t index, T elevation)
        {
            m_encoders[index].setElevation(elevation);
        }

        //! @brief Mutes or unmutes a source.
        //! @param     index	The index of the source.
        //! @param     muted	The mute status.
        inline void setMute(size_t index, const bool muted)
        {
            m_encoders[index].setMute(muted);
        }
        
        //! @brief Applies a fisheye effect on the sources positions.
        //! @details The fishEye value is between \f$0\f$ and \f$1\f$. At \f$0\f$, the sound
        //! field all the sources are equally dispatched around the equator and at \f$1\f$,
        //! the the sources are are concentrer in front of the audience.
        inline void setFisheye(T fisheye)
        {
            const size_t nsources = getNumberOfSources();
            const T factor = 1. - ((fisheye > T(1.)) ? T(1.) : ((fisheye < T(0.)) ? T(0.) : fisheye));
            for(size_t i = 0; i < nsources; ++i)
            {
                const T azimuth = static_cast<T>(i) / static_cast<T>(nsources) * static_cast<T>(HOA_2PI);
                if(azimuth < T(HOA_PI))
                {
                    m_encoders[i].setAzimuth(azimuth * factor);
                }
                else
                {
                    m_encoders[i].setAzimuth(T(HOA_2PI) - ((T(HOA_2PI) - azimuth) * factor));
                }
            }
        }

        //! @brief Returns the radius of a source.
        inline T getRadius(size_t index) const
        {
            return m_encoders[index].getRadius();
        }
        
        //! @brief Returns the widening value of a source.
        inline T getWidening(size_t index) const
        {
            return m_encoders[index].getRadius() > T(1.) ? T(1.) : m_encoders[index].getRadius();
        }
        
        //! @brief Returns the azimuth of a source.
        inline T getAzimuth(size_t index) const
        {
            return m_encoders[index].getAzimuth();
        }
        
        //! @brief Returns the elevation of a source.
        inline T getElevation(size_t index) const
        {
            return m_encoders[index].getElevation();
        }
        
        //! @brief Returns the mute status of a source.
        inline bool getMute(size_t index) const
        {
            return m_encoders[index].getMute();
        }

        //! @brief The method performs the encoding of the harmonics signal.
        //! @details The inputs array contains the samples of the sources to encode and the
        //! outputs array contains the spherical harmonics samples thus the minimum size of
        //! the array must be the number of sources and the number of harmonics.
        //! @param input   The inputs array.
        //! @param outputs The outputs array.
        void process(const T* input, T* outputs) noexcept override
        {
            const size_t nencoders  = m_encoders.size();
            const size_t nharmos    = ProcessorHarmonics<D, T>::getNumberOfHarmonics();
            if(m_encoders[0].getMute())
            {
                memset(outputs, 0, nharmos * sizeof(T));
            }
            else
            {
                m_encoders[0].process(input, outputs);
            }
            
            for(size_t i = 1; i < nencoders; i++)
            {
                m_encoders[i].process(++input, m_temp);
                sigadd(nharmos, m_temp, outputs);
            }
        }
        
    private:
        
        static inline void sigadd(const size_t size, const T* in, T* out) noexcept
        {
            for(size_t i = size>>3; i; --i, in += 8, out += 8)
            {
                out[0] += in[0]; out[1] += in[1]; out[2] += in[2]; out[3] += in[3];
                out[4] += in[4]; out[5] += in[5]; out[6] += in[6]; out[7] += in[7];
            }
            for(size_t i = size&7; i; --i, in++, out++)
            {
                out[0] += in[0];
            }
        }
        
        struct EncoderWrap
        {
            EncoderWrap()
            : encoder(nullptr)
            , muted(false)
            {}
            
            ~EncoderWrap() = default;
            
            inline void setRadius(T radius) noexcept { encoder->setRadius(radius); }
            inline void setAzimuth(T azimuth) noexcept { encoder->setAzimuth(azimuth); }
            inline void setElevation(T elevation) noexcept { encoder->setElevation(elevation); }
            inline void setMute(bool m) noexcept { muted = m; }
            inline T getRadius() const { return encoder->getRadius(); }
            inline T getAzimuth() const { return encoder->getAzimuth(); }
            inline T getElevation() const { return encoder->getElevation(); }
            inline bool getMute() const { return muted; }
            
            inline void process(const T* input, T* outputs) noexcept
            {
                if(!muted)
                {
                    encoder->process(input, outputs);
                }
            }
            
            Encoder<D, T>* encoder = nullptr;
            bool muted = false;
        };
        
    private:
        
        T* m_temp = nullptr;
        std::vector<EncoderWrap> m_encoders {};
    };
}
