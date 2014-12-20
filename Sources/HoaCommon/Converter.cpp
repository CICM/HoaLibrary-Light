/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include "Converter.h"

namespace Hoa
{
    Converter::Converter(unsigned int order, Mode inputMode, Mode outputMode)
    {
        m_order					= order;
        m_input_mode            = inputMode;
        m_output_mode           = outputMode;
        if((m_input_mode == FUMA_2D || m_input_mode == FUMA_3D || m_output_mode == FUMA_2D || m_output_mode == FUMA_3D) && m_order > 3)
            m_order = 3;
        
        if(m_input_mode == HOA2D || m_input_mode == FUMA_2D || m_input_mode == ACN_N2D || m_input_mode == ACN_SN2D)
            m_number_of_inputs_harmonics = m_order * 2 + 1;
        else
            m_number_of_inputs_harmonics = (m_order + 1) * (m_order + 1);
        
        if(m_output_mode == HOA2D || m_output_mode == FUMA_2D || m_output_mode == ACN_N2D || m_output_mode == ACN_SN2D)
            m_number_of_outputs_harmonics = m_order * 2 + 1;
        else
            m_number_of_outputs_harmonics = (m_order + 1) * (m_order + 1);
        
        m_harmonics_arguments   = new int[m_number_of_harmonics];
        m_harmonics_bands      = new unsigned int[m_number_of_harmonics];
        for(unsigned int i = 0; i < m_number_of_harmonics; i++)
        {
            m_harmonics_bands[i] = sqrtf((float)i);
            int index = i - (m_harmonics_bands[i]) * (m_harmonics_bands[i]);
            
            int index2 = (index - 1) / 2. + 1.;
            if (index % 2 == 1)
                index2 = -index2;
            m_harmonics_arguments[i] = index2;
        }
    }

    Converter::~Converter()
    {
        delete [] m_harmonics_arguments;
        delete [] m_harmonics_bands;
    }
}

