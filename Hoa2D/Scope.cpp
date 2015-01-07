/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include "Scope.h"

namespace Hoa2D
{
    Scope::Scope(unsigned int order, unsigned int numberOfPoints) : Ambisonic(order)
    {
        m_number_of_points = numberOfPoints;
		m_decoder   = new DecoderRegular(order, m_number_of_points);
        m_matrix    = new double[m_number_of_points];
        m_harmonics = new double[m_number_of_harmonics];
 
        for(unsigned int i = 0; i < m_number_of_points; i++)
        {
            m_matrix[i] = 0.;
        }
    }
	
    void Scope::process(const float* inputs)
    {
        double max = 1.;
        for(unsigned int i = 0; i < m_number_of_harmonics; i++)
        {
            m_harmonics[i] = inputs[i];
        }
		m_decoder->process(m_harmonics, m_matrix);
        max = fabs(m_matrix[cblas_idamax(m_number_of_points, m_matrix, 1)]);
        if(max > 1.)
        {
            cblas_dscal(m_number_of_points, (1. / max), m_matrix, 1.);
        }
    }
    
    void Scope::process(const double* inputs)
    {
        double max = 1.;
        m_decoder->process(inputs, m_matrix);
        max = fabs(m_matrix[cblas_idamax(m_number_of_points, m_matrix, 1)]);
        if(max > 1.)
        {
            cblas_dscal(m_number_of_points, (1. / max), m_matrix, 1.);
        }
    }
    
    Scope::~Scope()
    {
		delete [] m_matrix;
        delete [] m_harmonics;
        delete m_decoder;
    }
	
}
