/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include "Ambisonic.h"

namespace Hoa2D
{
    Ambisonic::Ambisonic(unsigned int order)
    {
        m_order					= order;
        m_number_of_harmonics	= m_order * 2 + 1;
        
        m_harmonics_orders   = new long[m_number_of_harmonics];
        m_harmonics_orders[0] = 0;
        for(int i = 1, j = 1; i <= m_order; i++, j += 2)
        {
            m_harmonics_orders[j] = -i;
            m_harmonics_orders[j+1] = i;
        }
    }
    
    Ambisonic::~Ambisonic()
    {
        delete [] m_harmonics_orders;
    }
}

