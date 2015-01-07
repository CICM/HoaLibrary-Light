/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include "Recomposer.h"

namespace Hoa2D
{
    Recomposer::Recomposer(unsigned int order, unsigned int numberOfChannels) : Ambisonic(order), Planewaves(numberOfChannels)
    {
        m_harmonics_float           = new float[m_number_of_harmonics];
        m_harmonics_double          = new double[m_number_of_harmonics];
        m_recomposer_matrix_float   = new float[m_number_of_harmonics * m_number_of_channels];
        m_recomposer_matrix_double  = new double[m_number_of_harmonics * m_number_of_channels];
        m_encoders                  = new Encoder*[m_number_of_channels];
        m_widers                    = new Wider*[m_number_of_channels];
        for(unsigned int i = 0; i < m_number_of_channels; i++)
        {
            m_encoders[i]   = new Encoder(m_order);
            m_widers[i]     = new Wider(m_order);
            m_encoders[i]->setAzimuth((double)i / (double)m_number_of_channels * HOA_2PI);
            m_widers[i]->setWideningValue(1.);
        }

        for(unsigned int i = 0; i < m_number_of_channels; i++)
        {
            m_encoders[0]->setAzimuth(m_channels_azimuth[i]);
            m_encoders[0]->process(1., m_harmonics_double);
            for(unsigned int j = 0; j < m_number_of_harmonics; j++)
            {
                m_recomposer_matrix_float[j * m_number_of_channels + i] = m_recomposer_matrix_double[j * m_number_of_channels + i] = m_harmonics_double[j];
            }
        }
    }

    void Recomposer::setAzimuth(unsigned int index, const double azimuth)
    {
        assert(index < m_number_of_channels);
        m_encoders[index]->setAzimuth(azimuth);
    }
	
    void Recomposer::setWideningValue(unsigned int index, const double value)
    {
        assert(index < m_number_of_channels);
        m_widers[index]->setWideningValue(value);
    }
    
    void Recomposer::setFisheye(const double fisheye)
    {
        double factor = 1. - clip_minmax(fisheye, 0., 1.);
        for(unsigned int i = 0; i < m_number_of_channels; i++)
        {
            double azimuth = (double)i / (double)m_number_of_channels * HOA_2PI;
            if(azimuth < HOA_PI)
                azimuth *= factor;
            else
                azimuth = HOA_2PI - ((HOA_2PI - azimuth) * factor);
            m_encoders[i]->setAzimuth(azimuth);
        }
    }
    
    void Recomposer::processFixe(const float* input, float* output)
	{
		cblas_sgemv(CblasRowMajor, CblasNoTrans, m_number_of_harmonics, m_number_of_channels, 1.f, m_recomposer_matrix_float, m_number_of_channels, input, 1, 0.f, output, 1);
	}
	
	void Recomposer::processFixe(const double* input, double* output)
	{
		cblas_dgemv(CblasRowMajor, CblasNoTrans, m_number_of_harmonics, m_number_of_channels, 1.f, m_recomposer_matrix_double, m_number_of_channels, input, 1, 0.f, output, 1);
	}
    
    void Recomposer::processFisheye(const float* inputs, float* outputs)
	{
		m_encoders[0]->process(inputs[0], outputs);
        for(unsigned int i = 1; i < m_number_of_channels; i++)
        {
            m_encoders[i]->process(inputs[i], m_harmonics_float);
            cblas_saxpy(m_number_of_harmonics, 1., m_harmonics_float, 1, outputs, 1);
        }
	}
	
	void Recomposer::processFisheye(const double* inputs, double* outputs)
	{
		m_encoders[0]->process(inputs[0], outputs);
        for(unsigned int i = 1; i < m_number_of_channels; i++)
        {
            m_encoders[i]->process(inputs[i], m_harmonics_double);
            cblas_daxpy(m_number_of_harmonics, 1., m_harmonics_double, 1, outputs, 1);
        }
	}
    
    void Recomposer::processFree(const float* inputs, float* outputs)
	{
		m_encoders[0]->process(inputs[0], m_harmonics_float);
        m_widers[0]->process(m_harmonics_float, outputs);
        for(unsigned int i = 1; i < m_number_of_channels; i++)
        {
            m_encoders[i]->process(inputs[i], m_harmonics_float);
            m_widers[i]->process(m_harmonics_float, m_harmonics_float);
            cblas_saxpy(m_number_of_harmonics, 1.f, m_harmonics_float, 1, outputs, 1);
        }
	}
	
	void Recomposer::processFree(const double* inputs, double* outputs)
	{
		m_encoders[0]->process(inputs[0], m_harmonics_double);
        m_widers[0]->process(m_harmonics_double, outputs);
        for(unsigned int i = 1; i < m_number_of_channels; i++)
        {
            m_encoders[i]->process(inputs[i], m_harmonics_double);
            m_widers[i]->process(m_harmonics_double, m_harmonics_double);
            cblas_daxpy(m_number_of_harmonics, 1., m_harmonics_double, 1, outputs, 1);
        }
	}
	
	Recomposer::~Recomposer()
	{
        delete [] m_harmonics_double;
        delete [] m_harmonics_float;
        delete [] m_recomposer_matrix_double;
        delete [] m_recomposer_matrix_float;
        for(unsigned int i = 0; i < m_number_of_channels; i++)
        {
            delete m_encoders[i];
            delete m_widers[i];
        }
        delete [] m_encoders;
        delete [] m_widers;
	}
}

