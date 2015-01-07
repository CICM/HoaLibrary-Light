/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include "Scope_3D.h"

namespace Hoa3D
{
    Scope::Scope(unsigned int order, unsigned int numberOfRows, unsigned int numberOfColumns) : Ambisonic(order)
    {
        m_number_of_rows = numberOfRows;
        m_number_of_columns = numberOfColumns;
		m_decoder   = new DecoderRegular(order, m_number_of_rows * m_number_of_columns);
        m_matrix    = new double[m_number_of_rows * m_number_of_columns];
        m_harmonics = new double[m_number_of_harmonics];
 
        for(unsigned int i = 0; i < m_number_of_rows; i++)
        {
            for(unsigned int j = 0; j < m_number_of_columns; j++)
            {
				m_decoder->setChannelPosition(i * m_number_of_columns + j,
											  (double)j * HOA_2PI / (double)m_number_of_columns,
											  (double)i * HOA_PI / (double)(m_number_of_rows - 1) - HOA_PI2);
				m_matrix[i * m_number_of_columns + j] = 0.;
            }
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
        max = fabs(m_matrix[cblas_idamax(m_number_of_rows * m_number_of_columns, m_matrix, 1)]);
        if(max > 1.)
        {
            cblas_dscal(m_number_of_rows * m_number_of_columns, (1. / max), m_matrix, 1.);
        }
    }
    
    void Scope::process(const double* inputs)
    {
        double max = 1.;
        m_decoder->process(inputs, m_matrix);
        max = fabs(m_matrix[cblas_idamax(m_number_of_rows * m_number_of_columns, m_matrix, 1)]);
        if(max > 1.)
        {
            cblas_dscal(m_number_of_rows * m_number_of_columns, (1. / max), m_matrix, 1.);
        }
    }
    
    Scope::~Scope()
    {
		delete [] m_matrix;
        delete [] m_harmonics;
        delete m_decoder;
    }
	
}
