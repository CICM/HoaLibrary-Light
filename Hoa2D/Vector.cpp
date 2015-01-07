/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include "Vector.h"

namespace Hoa2D
{
    Vector::Vector(unsigned int numberOfChannels) : Planewaves(numberOfChannels)
    {
        m_channels_float = new float[m_number_of_channels];
        m_channels_double = new double[m_number_of_channels];
        m_channels_abscissa_float = new float[m_number_of_channels];
        m_channels_abscissa_double = new double[m_number_of_channels];
        m_channels_ordinate_float = new float[m_number_of_channels];
        m_channels_ordinate_double = new double[m_number_of_channels];
        m_offset    = 0.;
        for(unsigned int i = 0; i < m_number_of_channels; i++)
        {
            m_channels_abscissa_float[i] = m_channels_abscissa_double[i] = getChannelAbscissa(i);
            m_channels_ordinate_float[i] = m_channels_ordinate_double[i] = getChannelOrdinate(i);
        }
    }

    void Vector::setChannelsOffset(double offset)
    {
        m_offset = wrap_twopi(offset);
        for (int i = 0; i < m_number_of_channels; i++)
		{
			m_channels_abscissa_float[i] = m_channels_abscissa_double[i] = abscissa(1., m_channels_azimuth[i] + m_offset);
			m_channels_ordinate_float[i] = m_channels_ordinate_double[i] = ordinate(1., m_channels_azimuth[i] + m_offset);
		}
    }

    void Vector::setChannelAzimuth(unsigned int index, double azimuth)
    {
        Planewaves::setChannelAzimuth(index, azimuth);
        m_channels_abscissa_float[index] = m_channels_abscissa_double[index] = abscissa(1., m_channels_azimuth[index] + m_offset);
        m_channels_ordinate_float[index] = m_channels_ordinate_double[index] = ordinate(1., m_channels_azimuth[index] + m_offset);
    }

	void Vector::setChannelsAzimuth(double* azimuths)
    {
        Planewaves::setChannelsAzimuth(azimuths);
		for (int i = 0; i < m_number_of_channels; i++)
		{
			m_channels_abscissa_float[i] = m_channels_abscissa_double[i] = abscissa(1., m_channels_azimuth[i] + m_offset);
			m_channels_ordinate_float[i] = m_channels_ordinate_double[i] = ordinate(1., m_channels_azimuth[i] + m_offset);
		}
    }

    void Vector::processVelocity(const float* inputs, float* outputs)
    {
        float veclocitySum, velocityAbscissa, velocityOrdinate;
        veclocitySum = velocityAbscissa = velocityOrdinate = 0.f;
		
        for(int i = 0; i < m_number_of_channels; i++)
            veclocitySum += inputs[i];
		
        velocityAbscissa = cblas_sdot(m_number_of_channels, inputs, 1, m_channels_abscissa_float, 1);
        velocityOrdinate = cblas_sdot(m_number_of_channels, inputs, 1, m_channels_ordinate_float, 1);
        if(veclocitySum)
        {
            outputs[0] = velocityAbscissa / veclocitySum;
            outputs[1] = velocityOrdinate / veclocitySum;
        }
        else
        {
            outputs[0] = 0.;
            outputs[1] = 0.;
        }
    }

    void Vector::processVelocity(const double* inputs, double* outputs)
    {
        double veclocitySum = 0, velocityAbscissa = 0, velocityOrdinate = 0;
		
        for(int i = 0; i < m_number_of_channels; i++)
            veclocitySum += inputs[i];
		
        velocityAbscissa = cblas_ddot(m_number_of_channels, inputs, 1, m_channels_abscissa_double, 1);
        velocityOrdinate = cblas_ddot(m_number_of_channels, inputs, 1, m_channels_ordinate_double, 1);

        if(veclocitySum)
        {
            outputs[0] = velocityAbscissa / veclocitySum;
            outputs[1] = velocityOrdinate / veclocitySum;
        }
        else
        {
            outputs[0] = 0.;
            outputs[1] = 0.;
        }
    }

    void Vector::processEnergy(const float* inputs, float* outputs)
    {
        float energySum = 0.f, energyAbscissa = 0.f, energyOrdinate = 0.f;
        cblas_scopy(m_number_of_channels, inputs, 1, m_channels_float, 1);
		
        for(int i = 0; i < m_number_of_channels; i++)
            m_channels_float[i] *= m_channels_float[i];

        energySum = cblas_sasum(m_number_of_channels, m_channels_float, 1);
        energyAbscissa = cblas_sdot(m_number_of_channels, m_channels_float, 1, m_channels_abscissa_float, 1);
        energyOrdinate = cblas_sdot(m_number_of_channels, m_channels_float, 1, m_channels_ordinate_float, 1);

        if(energySum)
        {
            outputs[0] = energyAbscissa / energySum;
            outputs[1] = energyOrdinate / energySum;
        }
        else
        {
            outputs[0] = 0.;
            outputs[1] = 0.;
        }
    }

    void Vector::processEnergy(const double* inputs, double* outputs)
    {
        double energySum = 0., energyAbscissa = 0., energyOrdinate = 0.;

        cblas_dcopy(m_number_of_channels, inputs, 1, m_channels_double, 1);
		
        for(int i = 0; i < m_number_of_channels; i++)
            m_channels_double[i] *= m_channels_double[i];

        energySum = cblas_dasum(m_number_of_channels, m_channels_double, 1);
        energyAbscissa = cblas_ddot(m_number_of_channels, m_channels_double, 1, m_channels_abscissa_double, 1);
        energyOrdinate = cblas_ddot(m_number_of_channels, m_channels_double, 1, m_channels_ordinate_double, 1);

        if(energySum)
        {
            outputs[0] = energyAbscissa / energySum;
            outputs[1] = energyOrdinate / energySum;
        }
        else
        {
            outputs[0] = 0.;
            outputs[1] = 0.;
        }
    }

    void Vector::process(const float* inputs, float* outputs)
    {
        processVelocity(inputs, outputs);
        processEnergy(inputs, outputs+2);
    }

    void Vector::process(const double* inputs, double* outputs)
    {
        processVelocity(inputs, outputs);
        processEnergy(inputs, outputs+2);
    }

    Vector::~Vector()
    {
        delete [] m_channels_float;
        delete [] m_channels_double;
        delete [] m_channels_ordinate_float;
        delete [] m_channels_ordinate_double;
        delete [] m_channels_abscissa_float;
        delete [] m_channels_abscissa_double;
    }
}

