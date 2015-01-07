/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include "Wider_3D.h"

namespace Hoa3D
{
    Wider::Wider(unsigned int order) : Ambisonic(order)
    {
        m_wide              = NUMBEROFLINEARPOINTS - 1;
        m_wide_matrix       = new double*[NUMBEROFLINEARPOINTS];
        
        for(int j = 0; j < NUMBEROFLINEARPOINTS; j++)
        {
            m_wide_matrix[j]    = new double[m_number_of_harmonics];
        }
        
        double weight_order = log((double)(m_order + 1) * (m_order + 1));
        
        for(int j = 0; j < NUMBEROFLINEARPOINTS; j++)
        {
            m_wide_matrix[j][0] = (1. - ((double)j / (double)(NUMBEROFLINEARPOINTS-1))) * weight_order + 1.;
        }
        for(unsigned int i = 1; i < m_number_of_harmonics; i++)
        {
            double minus =  clip_min(log((double)getHarmonicDegree(i)), 0.);
            minus = -minus;
            
            double dot	= clip_min(log((double)getHarmonicDegree(i) + 1.), 0.);
            dot += minus;
            dot  = 1. / dot;
            
            for(int j = 0; j < NUMBEROFLINEARPOINTS; j++)
            {
                double weight = (1. - ((double)j / (double)(NUMBEROFLINEARPOINTS-1))) * weight_order + 1.;
                double scale = ((double)j / (double)(NUMBEROFLINEARPOINTS-1)) * weight_order;
                double new_weight = (minus + scale) * dot;
                new_weight = clip_minmax(new_weight, 0., 1.);
                m_wide_matrix[j][i] = new_weight * weight;
            }
        }
    }
    
    void Wider::setWideningValue(const double value)
    {
         m_wide = clip_minmax(value, 0., 1.) * (double)(NUMBEROFLINEARPOINTS - 1);
    }
    
    void Wider::process(const float* inputs, float* outputs)
    {
        for(unsigned int i = 0; i < m_number_of_harmonics; i++)
            outputs[i] = inputs[i] * m_wide_matrix[m_wide][i];
    }
    
    void Wider::process(const double* inputs, double* outputs)
    {
        for(unsigned int i = 0; i < m_number_of_harmonics; i++)
            outputs[i] = inputs[i] * m_wide_matrix[m_wide][i];
    }
    
    Wider::~Wider()
    {
        delete [] m_wide_matrix;
    }
}

