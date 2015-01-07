/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include "Planewaves.h"

namespace Hoa2D
{
    Planewaves::Planewaves(unsigned int numberOfChannels)
    {
		assert(numberOfChannels > 0);
		m_number_of_channels    = numberOfChannels;
        m_channels_azimuth      = new double[m_number_of_channels];
        for(unsigned int i = 0; i < m_number_of_channels; i++)
        {
            m_channels_azimuth[i] = (double)i / (double)m_number_of_channels * HOA_2PI;
        }
    }
	
	void Planewaves::setChannelAzimuth(unsigned int index, double azimuth)
	{
		assert(index < m_number_of_channels);
		m_channels_azimuth[index] = wrap_twopi(azimuth);
        vector_sort(m_number_of_channels, m_channels_azimuth);
	}
    
    void Planewaves::setChannelsAzimuth(double* azimuths)
	{
		for(unsigned int i = 0; i < m_number_of_channels; i++)
            m_channels_azimuth[i] = wrap_twopi(azimuths[i]);
        vector_sort(m_number_of_channels, m_channels_azimuth);
	}
    
    Planewaves::~Planewaves()
    {
		delete [] m_channels_azimuth;
    }
}

