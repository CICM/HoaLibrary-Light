/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include "Planewaves_3D.h"

namespace Hoa3D
{
    Planewaves::Planewaves(unsigned int numberOfChannels)
    {
		assert(numberOfChannels > 0);
		m_number_of_channels    = numberOfChannels;
        m_channels_azimuth      = new double[m_number_of_channels];
        m_channels_elevation    = new double[m_number_of_channels];
        m_channels_rotated_azimuth      = new double[m_number_of_channels];
        m_channels_rotated_elevation    = new double[m_number_of_channels];
        sphere_discretize(numberOfChannels, m_channels_azimuth, m_channels_elevation);
        setChannelsRotation(0, 0, 0);
    }

	void Planewaves::setChannelsPosition(double* azimuths, double* elevations)
	{
		for(unsigned int i = 0; i < m_number_of_channels; i++)
		{
			m_channels_azimuth[i] = wrap_twopi(azimuths[i]);
			m_channels_elevation[i] = wrap(elevations[i], -HOA_PI, HOA_PI);
			if(m_channels_elevation[i] > HOA_PI2)
			{
				m_channels_azimuth[i] = wrap_twopi(m_channels_azimuth[i] + HOA_PI);
				m_channels_elevation[i] = HOA_PI2 - (m_channels_elevation[i] - HOA_PI2);
			}
			else if (m_channels_elevation[i] < -HOA_PI2)
			{
				m_channels_azimuth[i] = wrap_twopi(m_channels_azimuth[i] + HOA_PI);
				m_channels_elevation[i] = -HOA_PI2 + (-m_channels_elevation[i] + HOA_PI2);
			}
            rotateChannel(i);
		}
	}
	
	void Planewaves::setChannelPosition(unsigned int index, double azimuth, double elevation)
	{
		assert(index < m_number_of_channels);
		m_channels_azimuth[index] = wrap_twopi(azimuth);
        m_channels_elevation[index] = wrap(elevation, -HOA_PI, HOA_PI);

		if(m_channels_elevation[index] > HOA_PI2)
		{
			m_channels_azimuth[index] = wrap_twopi(m_channels_azimuth[index] + HOA_PI);
			m_channels_elevation[index] -= HOA_PI2;
		}
		else if (m_channels_elevation[index] < -HOA_PI2)
		{
			m_channels_azimuth[index] = wrap_twopi(m_channels_azimuth[index] + HOA_PI);
			m_channels_elevation[index] += HOA_PI2;
		}
        rotateChannel(index);
	}
    
    void Planewaves::setChannelsRotation(double axis_x, double axis_y, double axis_z)
    {
        m_channels_rotation_x = wrap_twopi(axis_x);
        m_channels_rotation_y = wrap_twopi(axis_y);
        m_channels_rotation_z = wrap_twopi(axis_z);
        
        for(unsigned int i = 0; i < m_number_of_channels; i++)
            rotateChannel(i);
    }
    
    void Planewaves::rotateChannel(unsigned int index)
    {
        double x = getChannelAbscissa(index);
        double y = getChannelOrdinate(index);
        double z = getChannelHeight(index);
        
        // Rotation around x
        double cosAngle = cos(m_channels_rotation_x);
        double sinAngle = sin(m_channels_rotation_x);
        double ry = y * cosAngle - z * sinAngle;
        double rz = y * sinAngle + z * cosAngle;
        y = ry;
        z = rz;
        
        // Rotation around z
        cosAngle = cos(m_channels_rotation_z);
        sinAngle = sin(m_channels_rotation_z);
        double rx = x * cosAngle - y * sinAngle;
        ry = x * sinAngle + y * cosAngle;
        x = rx;
        y = ry;
        
        // Rotation around y
        cosAngle = cos(m_channels_rotation_y);
        sinAngle = sin(m_channels_rotation_y);
        rx = x * cosAngle - z * sinAngle;
        rz = x * sinAngle + z * cosAngle;
        x = rx;
        z = rz;
        
        m_channels_rotated_azimuth[index] = azimuth(x, y, z);
        m_channels_rotated_elevation[index] = elevation(x, y, z);
    }
    
    Planewaves::~Planewaves()
    {
		delete [] m_channels_azimuth;
        delete [] m_channels_elevation;
        delete [] m_channels_rotated_azimuth;
        delete [] m_channels_rotated_elevation;
    }
}

