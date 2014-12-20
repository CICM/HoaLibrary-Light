/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include "Ambisonic_3D.h"

namespace Hoa3D
{
    Ambisonic::Ambisonic(unsigned int order)
    {
        m_order					= order;
        m_number_of_harmonics	= (m_order + 1) * (m_order + 1);
        
        m_harmonics_orders   = new int[m_number_of_harmonics];
        m_harmonics_degrees      = new unsigned int[m_number_of_harmonics];
		
        for(unsigned int i = 0; i < m_number_of_harmonics; i++)
        {
			// [0 0] [1 -1] [1 0] [1 1] [2 -2] [2 -1] [2 0] [2 1] [2 2]
			
			int degree = m_harmonics_degrees[i] = sqrtf((float)i);
            int index = i - (degree * degree);
			m_harmonics_orders[i] = -(degree - index);
        }
    }

    Ambisonic::~Ambisonic()
    {
        delete [] m_harmonics_orders;
        delete [] m_harmonics_degrees;
    }
}

