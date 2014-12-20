/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include "Optim.h"

namespace Hoa2D
{
    Optim::Optim(unsigned int order, Mode mode) : Ambisonic(order)
    {
        m_harmonics = new double[m_number_of_harmonics];
        setMode(mode);
    }
    
    void Optim::setMode(Mode mode)
    {
        long double temp1;
        long double temp2;
        m_mode = mode;
        if(m_mode == Basic)
        {
            for(unsigned int i = 0; i < m_number_of_harmonics; i++)
            {
                m_harmonics[i] = 1.;
            }
        }
        else if (m_mode == MaxRe)
        {
            for(unsigned int i = 0; i < m_number_of_harmonics; i++)
            {
                m_harmonics[i] = cos(fabs((double)getHarmonicDegree(i)) * HOA_PI / (double)(2. * m_order + 2));
            }
        }
        else
        {
            for(unsigned int i = 0; i < m_number_of_harmonics; i++)
            {
                temp1 = (long double)factorial(m_order) / (long double)factorial(m_order + getHarmonicDegree(i) + 1.);
                temp2 = (long double)factorial(m_order + 1.) / (long double)factorial(m_order - fabs((double)getHarmonicDegree(i)));
                m_harmonics[i] = temp1 * temp2;
            }
        }
    }
    
    void Optim::process(const float* inputs, float* outputs)
    {
        for(unsigned int i = 0; i < m_number_of_harmonics; i++)
            outputs[i] = inputs[i] * m_harmonics[i];
    }
    
    void Optim::process(const double* inputs, double* outputs)
    {
        for(unsigned int i = 0; i < m_number_of_harmonics; i++)
            outputs[i] = inputs[i] * m_harmonics[i];
    }
    
    Optim::~Optim()
    {
        delete [] m_harmonics;
    }
}

