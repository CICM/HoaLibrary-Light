/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_2D_PLANEWAVES
#define DEF_HOA_2D_PLANEWAVES

#include "../Hoa.h"

namespace Hoa2D
{
    //! The planewaves class.
    /**
     The planewaves classes, that process on a set of channels (or planewaves), inherit from this class. It store basic informations like the number of channels, the coordinates and the names of channels.
     */
    class Planewaves
    {
	protected:
        
		unsigned int    m_number_of_channels;
        double*         m_channels_azimuth;
        
        //! Set the azimuth of a channel.
        /** Set the azimuth of a channel. The azimuth is in radian between 0 and 2 Pi, O is the front of the soundfield and Pi is the back of the sound field. The maximum index must be the number of channel - 1.
         
            @param     index		The index of the channel.
            @param     azimuth		The azimuth.
         */
		void setChannelAzimuth(unsigned int index, double azimuth);
		
        //! Set the azimtuh of all the channels.
        /** Set the azimtuh of all the channels. It is more efficient to set all the channels azimuths at the same time because even if only one channel has changed, all the decoding matrix have to be recomputed. The azimuths are in radian between 0 and 2 Pi, O is the front of the soundfield and Pi is the back of the sound field. The azimtuhs array must have a minimum size of the number of channels.
         
            @param     azimuths		The azimuths array.
         
            @see    setChannelAzimuth
         */
        void setChannelsAzimuth(double* azimuths);
    public:
        
		//! The planewaves constructor.
        /** The lanewaves constructor allocates and initializes the general member values depending on a number of channels. The number of loudspkeakers must a least 1.
         
            @param     numberOfChannels	The number of channels.
         */
		Planewaves(unsigned int numberOfChannels);
        
        //! The planewaves destructor.
        /** The Planewaves destructor free the memorie allocated.
         */
        ~Planewaves();
        
        //! Retrieve the number of channels.
		/** Retrieve the number of channels of the planewave class.
            
            @return The number of channels.
         */
		inline unsigned int getNumberOfChannels() const
        {
            return m_number_of_channels;
        }
        
        //! Retrieve the azimuth of a channel.
        /** Retrieve the azimuth of a channel. The azimuth of the channel is in radian, 0 radian is at the front of the soundfield and Pi is at the back of the sound field. The maximum index must be the number of channels - 1.
         
            @param      index   The index of the channel.
            @return     The azimuth of the channel if the channel exists, otherwise the function generates an error.
         
            @see getChannelAbscissa
            @see getChannelOrdinate
            @see getChannelName
         */
		inline double getChannelAzimuth(unsigned int index) const
        {
            assert(index < m_number_of_channels);
            return m_channels_azimuth[index];
        }

		
        //! Retrieve the abscissa of a channel.
		/** Retrieve the abscissa of a channel. The abscissa is between -1 and 1, -1 is the left of the soundfield, 0 is the center of the soundfield and 1 is the right of the soundfield. The maximum index must be the number of channels - 1.
         
            @param     index    The index of the channel.
            @return    The abscissa of the channel if the channel exists, otherwise the function generates an error.
         
            @see getChannelAzimuth
            @see getChannelOrdinate
            @see getChannelName
         */
		inline double getChannelAbscissa(unsigned int index) const
        {
            assert(index < m_number_of_channels);
            return abscissa(1., m_channels_azimuth[index]);
        }
		
        //! Retrieve the ordinate of a channel.
		/** Retrieve the ordinate of a channel. The ordinate is between -1 and 1, -1 is the back of the soundfield, 0 is the center of the soundfield and 1 is the front of the soundfield. The maximum index must be the number of channels - 1.
         
            @param     index	The index of the channel.
            @return    The ordinate of the channel if the channel exists, otherwise the function generates an error.
         
            @see getChannelAzimuth
            @see getChannelAbscissa
            @see getChannelName
         */
		inline double getChannelOrdinate(unsigned int index) const
        {
            assert(index < m_number_of_channels);
            return ordinate(1., m_channels_azimuth[index]);
        }
        
        //! Retrieve a name for a channel.
        /** Retrieve a name for a channel in a std::string format that will be "Channel index azimuth (in degrees)".
         
            @param     index	The index of a channel.
            @return    The method returns a name for the channel that contains its index and its azimuth if the channel exists, otherwise the function generates an error.
         
            @see getChannelAzimuth
            @see getChannelAbscissa
            @see getChannelOrdinate
         */
		inline std::string getChannelName(unsigned int index)
        {
            assert(index < m_number_of_channels);
            return "Channel " + int_to_string(index + 1) + " : " + int_to_string((int)(getChannelAzimuth(index) / HOA_2PI * 360.)) + "°";
        };
    };
	
}

#endif