/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include "Vector_3D.h"

namespace Hoa3D
{
    Vector::Vector(unsigned int numberOfChannels) : Planewaves(numberOfChannels)
    {
        m_channels_float = new float[m_number_of_channels];
        m_channels_double = new double[m_number_of_channels];
        m_channels_abscissa_float = new float[m_number_of_channels];
        m_channels_abscissa_double = new double[m_number_of_channels];
        m_channels_ordinate_float = new float[m_number_of_channels];
        m_channels_ordinate_double = new double[m_number_of_channels];
        m_channels_height_float = new float[m_number_of_channels];
        m_channels_height_double = new double[m_number_of_channels];
        setChannelPosition(0, m_channels_azimuth[0], m_channels_elevation[0]);
    }

    void Vector::setChannelPosition(unsigned int index, double azimuth, double elevation)
    {
        Planewaves::setChannelPosition(index, azimuth, elevation);
        for(unsigned int i = 0; i < m_number_of_channels; i++)
        {
            m_channels_abscissa_float[i] = m_channels_abscissa_double[i] = getChannelRotatedAbscissa(i);
            m_channels_ordinate_float[i] = m_channels_ordinate_double[i] = getChannelRotatedOrdinate(i);
            m_channels_height_float[i] = m_channels_height_double[i] = getChannelRotatedHeight(i);
        }
    }
    
    void Vector::setChannelsPosition(double* azimuths, double* elevations)
	{
		Planewaves::setChannelsPosition(azimuths, elevations);
		setChannelPosition(0, m_channels_azimuth[0], m_channels_elevation[0]);
	}
    
    void Vector::setChannelsRotation(double axis_x, double axis_y, double axis_z)
    {
        Planewaves::setChannelsRotation(axis_x, axis_y, axis_z);
        setChannelPosition(0, m_channels_azimuth[0], m_channels_elevation[0]);
    }

    void Vector::processVelocity(const float* inputs, float* outputs)
    {
        float veclocitySum = 0.f, velocityAbscissa = 0.f, velocityOrdinate = 0.f, velocityElevation = 0.f;

        veclocitySum = cblas_sasum(m_number_of_channels, inputs, 1);
        velocityAbscissa = cblas_sdot(m_number_of_channels, inputs, 1, m_channels_abscissa_float, 1);
        velocityOrdinate = cblas_sdot(m_number_of_channels, inputs, 1, m_channels_ordinate_float, 1);
        velocityElevation= cblas_sdot(m_number_of_channels, inputs, 1, m_channels_height_float, 1);
        if(veclocitySum)
        {
            outputs[0] = velocityAbscissa / veclocitySum;
            outputs[1] = velocityOrdinate / veclocitySum;
            outputs[2] = velocityElevation / veclocitySum;
        }
        else
        {
            outputs[0] = 0.;
            outputs[1] = 0.;
            outputs[2] = 0.;

        }
    }

    void Vector::processVelocity(const double* inputs, double* outputs)
    {
        double veclocitySum = 0., velocityAbscissa = 0., velocityOrdinate = 0., velocityElevation = 0.;
        veclocitySum = cblas_dasum(m_number_of_channels, inputs, 1);
        velocityAbscissa = cblas_ddot(m_number_of_channels, inputs, 1, m_channels_abscissa_double, 1);
        velocityOrdinate = cblas_ddot(m_number_of_channels, inputs, 1, m_channels_ordinate_double, 1);
        velocityElevation= cblas_ddot(m_number_of_channels, inputs, 1, m_channels_height_double, 1);

        if(veclocitySum)
        {
            outputs[0] = velocityAbscissa / veclocitySum;
            outputs[1] = velocityOrdinate / veclocitySum;
            outputs[2] = velocityElevation / veclocitySum;
        }
        else
        {
            outputs[0] = 0.;
            outputs[1] = 0.;
            outputs[2] = 0.;
        }
    }

    void Vector::processEnergy(const float* inputs, float* outputs)
    {
        float energySum = 0.f, energyAbscissa = 0.f, energyOrdinate = 0.f, energyElevation = 0.f;
        cblas_scopy(m_number_of_channels, inputs, 1, m_channels_float, 1);
        for(unsigned int i = 0; i < m_number_of_channels; i++)
            m_channels_float[i] *= m_channels_float[i];

        energySum = cblas_sasum(m_number_of_channels, m_channels_float, 1);
        energyAbscissa = cblas_sdot(m_number_of_channels, m_channels_float, 1, m_channels_abscissa_float, 1);
        energyOrdinate = cblas_sdot(m_number_of_channels, m_channels_float, 1, m_channels_ordinate_float, 1);
        energyElevation = cblas_sdot(m_number_of_channels, m_channels_float, 1, m_channels_height_float, 1);

        if(energySum)
        {
            outputs[0] = energyAbscissa / energySum;
            outputs[1] = energyOrdinate / energySum;
            outputs[2] = energyElevation / energySum;
        }
        else
        {
            outputs[0] = 0.;
            outputs[1] = 0.;
            outputs[2] = 0.;
        }
    }

    void Vector::processEnergy(const double* inputs, double* outputs)
    {
        double energySum = 0., energyAbscissa = 0., energyOrdinate = 0., energyElevation = 0.;

        cblas_dcopy(m_number_of_channels, inputs, 1, m_channels_double, 1);
        for(unsigned int i = 0; i < m_number_of_channels; i++)
            m_channels_double[i] *= m_channels_double[i];

        energySum = cblas_dasum(m_number_of_channels, m_channels_double, 1);
        energyAbscissa = cblas_ddot(m_number_of_channels, m_channels_double, 1, m_channels_abscissa_double, 1);
        energyOrdinate = cblas_ddot(m_number_of_channels, m_channels_double, 1, m_channels_ordinate_double, 1);
        energyElevation = cblas_ddot(m_number_of_channels, m_channels_double, 1, m_channels_height_double, 1);

        if(energySum)
        {
            outputs[0] = energyAbscissa / energySum;
            outputs[1] = energyOrdinate / energySum;
            outputs[2] = energyElevation / energySum;
        }
        else
        {
            outputs[0] = 0.;
            outputs[1] = 0.;
            outputs[2] = 0.;
        }
    }

    void Vector::process(const float* inputs, float* outputs)
    {
        processVelocity(inputs, outputs);
        processEnergy(inputs, outputs+3);
    }

    void Vector::process(const double* inputs, double* outputs)
    {
        processVelocity(inputs, outputs);
        processEnergy(inputs, outputs+3);
    }

    Vector::~Vector()
    {
        delete [] m_channels_float;
        delete [] m_channels_double;
        delete [] m_channels_ordinate_float;
        delete [] m_channels_ordinate_double;
        delete [] m_channels_abscissa_float;
        delete [] m_channels_abscissa_double;
        delete [] m_channels_height_float;
        delete [] m_channels_height_double;
    }
}

