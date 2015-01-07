/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include "Rotate.h"

namespace Hoa2D
{
    Rotate::Rotate(unsigned int order) : Ambisonic(order)
    {
        setYaw(0.);
    }
	
	void Rotate::setYaw(const double value)
    {
		m_yaw =	wrap_twopi(value);
		m_cosx    = cos(m_yaw);
        m_sinx    = sin(m_yaw);
    }
    
    void Rotate::process(const float* inputs, float* outputs)
    {
        float cos_x = m_cosx;
        float sin_x = m_sinx;
        float tcos_x = cos_x;
        float sig;
        outputs[0] = inputs[0];
        for(int i = 2; i < m_number_of_harmonics; i += 2)
        {
            sig = inputs[i-1];
            outputs[i-1] = sin_x * inputs[i] + cos_x * sig;
            outputs[i] = cos_x * inputs[i] - sin_x * sig;
            cos_x = tcos_x * m_cosx - sin_x * m_sinx; // cos(x + b) = cos(x) * cos(b) - sin(x) * sin(b)
            sin_x = tcos_x * m_sinx + sin_x * m_cosx; // sin(x + b) = cos(x) * sin(b) + sin(x) * cos(b)
            tcos_x = cos_x;
        }
    }
    
    void Rotate::process(const double* inputs, double* outputs)
    {
        double cos_x = m_cosx;
        double sin_x = m_sinx;
        double tcos_x = cos_x;
        double sig;
        outputs[0] = inputs[0];
        for(int i = 2; i < m_number_of_harmonics; i += 2)
        {
            sig = inputs[i-1];
            outputs[i-1] = sin_x * inputs[i] + cos_x * sig;
            outputs[i] = cos_x * inputs[i] - sin_x * sig;
            cos_x = tcos_x * m_cosx - sin_x * m_sinx; // cos(x + b) = cos(x) * cos(b) - sin(x) * sin(b)
            sin_x = tcos_x * m_sinx + sin_x * m_cosx; // sin(x + b) = cos(x) * sin(b) + sin(x) * cos(b)
            tcos_x = cos_x;
        }
    }
    
    Rotate::~Rotate()
    {
        ;
    }
	
}

