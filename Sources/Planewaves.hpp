/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef __DEF_HOA_PLANEWAVES__
#define __DEF_HOA_PLANEWAVES__

#include "HoaMath.hpp"

namespace hoa
{
    //! The planewaves class.
    /**
     The planewaves classes, that process on a set of channels (or planewaves), inherit from this class. It store basic informations like the number of channels, the coordinates and the names of channels.
     */
    template <typename T> class Planewaves
    {
    protected:
        struct Channel
        {
            unsigned long index;
            T             azimuth;
            static bool compare(Channel const& i, Channel const& j) noexcept
            {
                return i.index < j.index;
            }
        };
        
        const unsigned long m_number_of_channels;
        vector<Channel>     m_channels;
        T*                  m_channels_abscissa;
        T*                  m_channels_ordinate;
        T                   m_offset;
    public:
        
        //! The planewaves constructor.
        /** The lanewaves constructor allocates and initializes the general member values depending on a number of channels. The number of loudspkeakers must a least 1.
         @param     numberOfChannels	The number of channels.
         */
        Planewaves(unsigned long numberOfChannels) noexcept :
        m_number_of_channels(numberOfChannels)
        {
            m_channels_abscissa = new T[m_number_of_channels];
            m_channels_ordinate = new T[m_number_of_channels];
            for(unsigned long i = 0; i < m_number_of_channels; i++)
            {
                m_channels.push_back({i, (T)i / (T)m_number_of_channels * (T)HOA_2PI});
                m_channels_abscissa[i] = abscissa(1., m_channels[i].azimuth + m_offset);
                m_channels_ordinate[i] = ordinate(1., m_channels[i].azimuth + m_offset);
            }
            m_offset = 0.;
        }
        
        //! The planewaves destructor.
        /** The Planewaves destructor free the memorie allocated.
         */
        virtual ~Planewaves()
        {
            m_channels.clear();
            delete [] m_channels_abscissa;
            delete [] m_channels_ordinate;
        }
        
        //! Retrieve the number of channels.
        /** Retrieve the number of channels.
         @return The number of channels.
         */
        inline unsigned long getNumberOfChannels() const noexcept
        {
            return m_number_of_channels;
        }
        
        //! Set the offset of the channels.
        /**	Set the azimuth offset of the channels in radian.
         @param     offset		An azimuth value.
         */
        virtual inline void setChannelsOffset(T offset) noexcept
        {
            m_offset = wrap_twopi(offset);
            for(unsigned long i = 0; i < m_number_of_channels; i++)
            {
                m_channels_abscissa[i] = abscissa(1., m_channels[i].azimuth + m_offset);
                m_channels_ordinate[i] = ordinate(1., m_channels[i].azimuth + m_offset);
            }
        }
        
        //! Set the azimuth of a channel.
        /** Set the azimuth of a channel. The azimuth is in radian between 0 and 2 π, O is the front of the soundfield and π is the back of the sound field. The maximum index must be the number of channel - 1.
         @param     index		The index of the channel.
         @param     azimuth		The azimuth.
         */
        virtual inline void setChannelAzimuth(unsigned long index, T azimuth) noexcept
        {
            m_channels[index].azimuth = wrap_twopi(azimuth);
            m_channels_abscissa[index] = abscissa(1., m_channels[index].azimuth + m_offset);
            m_channels_ordinate[index] = ordinate(1., m_channels[index].azimuth + m_offset);
        }
        
        //! Get the offset of the channels.
        /**	Retreive the azimuth offset of the channels in radian.
         @return    The offset of the channels.
         */
        inline T getChannelsOffset() const noexcept
        {
            return m_offset;
        }
        
        //! Retrieve the azimuth of a channel.
        /** Retrieve the azimuth of a channel. The azimuth of the channel is in radian, 0 radian is at the front of the soundfield and π is at the back of the sound field. The maximum index must be the number of channels - 1.
         
         @param      index   The index of the channel.
         @return     The azimuth of the channel.
         */
        inline T getChannelAzimuth(unsigned long index) const noexcept
        {
            return m_channels[index].azimuth;
        }
        
        //! Retrieve the abscissa of a channel.
        /** Retrieve the abscissa of a channel. The abscissa is between -1 and 1, -1 is the left of the soundfield, 0 is the center of the soundfield and 1 is the right of the soundfield. The maximum index must be the number of channels - 1.
         @param     index    The index of the channel.
         @return    The abscissa of the channel.
         */
        inline T getChannelAbscissa(unsigned long index) const noexcept
        {
            return m_channels_abscissa[index];
        }
        
        //! Retrieve the ordinate of a channel.
        /** Retrieve the ordinate of a channel. The ordinate is between -1 and 1, -1 is the back of the soundfield, 0 is the center of the soundfield and 1 is the front of the soundfield. The maximum index must be the number of channels - 1.
         @param     index	The index of the channel.
         @return    The ordinate of the channel.
         */
        inline T getChannelOrdinate(unsigned long index) const noexcept
        {
            return m_channels_ordinate[index];
        }
        
        //! Retrieve a name for a channel.
        /** Retrieve a name for a channel in a std::string format that will be "Channel index azimuth (in degrees)".
         @param     index	The index of a channel.
         @return    The method returns a name for the channel.
         */
        inline std::string getChannelName(unsigned long index) const noexcept
        {
            return "Channel " + to_string(index + 1) + " : " + to_string((long)(m_channels[index].azimuth / HOA_2PI * 360.)) + "°";
        };
    };
}

#endif


