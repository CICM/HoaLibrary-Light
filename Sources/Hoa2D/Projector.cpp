/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include "Projector.h"

namespace Hoa2D
{
    Projector::Projector(unsigned int order, unsigned int numberOfChannels) : Ambisonic(order), Planewaves(numberOfChannels)
    {
        m_projector_matrix_double   = new double[m_number_of_channels * m_number_of_harmonics];
        m_projector_matrix_float    = new float[m_number_of_channels * m_number_of_harmonics];
        
        double*         m_harmonics_vector;
        Encoder*        m_encoder;
        m_harmonics_vector          = new double[m_number_of_harmonics];
        m_encoder                   = new Encoder(m_order);
        for(unsigned int i = 0; i < m_number_of_channels; i++)
        {
            m_channels_azimuth[i] = (double)i / (double)m_number_of_channels * HOA_2PI;
            m_encoder->setAzimuth(m_channels_azimuth[i]);
            m_encoder->process(1., m_harmonics_vector);
            
            m_projector_matrix_float[i * m_number_of_harmonics] = m_projector_matrix_double[i * m_number_of_harmonics] = 0.5 / (double)(m_order + 1.);
            for(unsigned int j = 1; j < m_number_of_harmonics; j++)
            {
                m_projector_matrix_float[i * m_number_of_harmonics + j] = m_projector_matrix_double[i * m_number_of_harmonics + j] = m_harmonics_vector[j] / (double)(m_order + 1.);
            }
        }
        delete [] m_harmonics_vector;
        delete m_encoder;
    }

    void Projector::process(const float* inputs, float* outputs)
	{
		cblas_sgemv(CblasRowMajor, CblasNoTrans, m_number_of_channels, m_number_of_harmonics, 1.f, m_projector_matrix_float, m_number_of_harmonics, inputs, 1, 0.f, outputs, 1);
	}
	
	void Projector::process(const double* inputs, double* outputs)
	{
		cblas_dgemv(CblasRowMajor, CblasNoTrans, m_number_of_channels, m_number_of_harmonics, 1., m_projector_matrix_double, m_number_of_harmonics, inputs, 1, 0., outputs, 1);
	}
	
	Projector::~Projector()
	{
		delete [] m_projector_matrix_double;
        delete [] m_projector_matrix_float;
	}
}

