/*
// Copyright (c) 2012-2016 CICM - Universite Paris 8 - Labex Arts H2H.
// Authors :
// 2012: Pierre Guillot, Eliott Paris & Julien Colafrancesco.
// 2012-2015: Pierre Guillot & Eliott Paris.
// 2015: Pierre Guillot & Eliott Paris & Thomas Le Meur (Light version)
// 2016: Pierre Guillot & Eliott Paris.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_MULTIENCODER_LIGHT
#define DEF_HOA_MULTIENCODER_LIGHT

#include "Hoa_Encoder.hpp"
#include <iostream>

namespace hoa
{
    //! The multi encoder class generates the harmonics for several signals according to an azimuth, an elevation and a radius for each one.
    /** The multi encoder should be used to encode several signals in the harmonics domain depending on an order of decomposition. It allows to control the azimuth, the elevation and the radius of each signal. The class uses a set of dc encoders.
     */
    template <Dimension D, typename T> class MultiEncoder : public ProcessorHarmonics<D, T>
    {
    public:

        //! @brief The map constructor.
        //! @param order            The order.
        //! @param numberOfSources  The number of sources.
        MultiEncoder(size_t order, size_t numberOfSources) hoa_noexcept : ProcessorHarmonics<D, T>(order),
        m_encoders(numberOfSources, order) {
            m_temp = new T[ProcessorHarmonics<D, T>::getNumberOfHarmonics()];
        }

        //! @brief The map destructor.
        ~MultiEncoder() {
            delete [] m_temp;
        }

        //! @brief Returns the number of sources.
        inline size_t getNumberOfSources() const { return m_encoders.size(); }
        
        //! @brief Sets the radius of a source.
        //! @param     index	The index of the source.
        //! @param     radius	The radius.
        inline void setRadius(size_t index, T radius) hoa_noexcept {
            m_encoders[index].setRadius(radius);
        }

        //! @brief Sets the azimuth of a source.
        //! @param     index	The index of the source.
        //! @param     azimuth	The azimuth.
        inline void setAzimuth(size_t index, T azimuth) hoa_noexcept {
            m_encoders[index].setAzimuth(azimuth);
        }

        //! @brief Sets the elevation of a source.
        //! @param     index	The index of the source.
        //! @param     elevation	The elevation.
        inline void setElevation(size_t index, T elevation) {
            m_encoders[index].setElevation(elevation);
        }


        //! @brief Mutes or unmutes a source.
        //! @param     index	The index of the source.
        //! @param     muted	The mute status.
        inline void setMute(size_t index, const bool muted) {
            m_encoders[index].setMute(muted);
        }

        //! @brief Returns the radius of a source.
        inline T getRadius(size_t index) const {
            return m_encoders[index].getRadius();
        }
        
        //! @brief Returns the azimuth of a source.
        inline T getAzimuth(size_t index) const {
            return m_encoders[index].getAzimuth();
        }
        
        //! @brief Returns the elevation of a source.
        inline T getElevation(size_t index) const {
            return m_encoders[index].getElevation();
        }
        
        //! @brief Returns the mute status of a source.
        inline bool getMute(size_t index) const {
            return m_encoders[index].getMute();
        }

        //! @brief The method performs the encoding of the harmonics signal.
        //! @details The inputs array contains the samples of the sources to encode and the
        //! outputs array contains the spherical harmonics samples thus the minimum size of
        //! the array must be the number of sources and the number of harmonics.
        //! @param input   The inputs array.
        //! @param outputs The outputs array.
        inline void process(const T* input, T* outputs) hoa_noexcept hoa_override
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
        
        static inline void sigadd(const size_t size, const T* in, T* out) hoa_noexcept
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
        
        struct EncoderWrap {
            Encoder<D, T>*  encoder;
            bool            muted;
            
            EncoderWrap(size_t order) : encoder(new Encoder<D, T>(order)), muted(false) {}
            ~EncoderWrap() { delete encoder; }
            
            inline void setRadius(T radius) hoa_noexcept { encoder->setRadius(radius); }
            inline void setAzimuth(T azimuth) hoa_noexcept { encoder->setAzimuth(azimuth); }
            inline void setElevation(T elevation) hoa_noexcept { encoder->setElevation(elevation); }
            inline void setMute(bool m) hoa_noexcept { muted = m; }
            inline T getRadius() const { return encoder->getRadius(); }
            inline T getAzimuth() const { return encoder->getAzimuth(); }
            inline T getElevation() const { return encoder->getElevation(); }
            inline T getMute() const { return muted; }
            inline void process(const T* input, T* outputs) hoa_noexcept {
                if(!muted) { encoder->process(input, outputs); }
            }
        };
        T*                       m_temp;
        std::vector<EncoderWrap> m_encoders;
    };
}

#endif
