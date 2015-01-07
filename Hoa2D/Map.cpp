/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include "Map.h"

namespace Hoa2D
{
    Map::Map(unsigned int order, unsigned int numberOfSources) : Ambisonic(order)
    {
        assert(numberOfSources > 0);
        
        m_number_of_sources = numberOfSources;
        m_gains             = new double[m_number_of_sources];
		m_muted				= new bool[m_number_of_sources];
        m_azimuth           = new double[m_number_of_sources];
        m_sinx              = new double[m_number_of_sources];
        m_cosx              = new double[m_number_of_sources];
        m_wide              = new long[m_number_of_sources];
        
        m_wide_matrix       = new double[NUMBEROFLINEARPOINTS * m_number_of_harmonics];
        
        double weight_order = log((double)(m_order + 1));
        
        for(unsigned int j = 0; j < NUMBEROFLINEARPOINTS; j++)
        {
            m_wide_matrix[j * m_number_of_harmonics] = (1. - ((double)j / (double)(NUMBEROFLINEARPOINTS-1))) * weight_order + 1.;
        }
        for(unsigned int i = 1; i < m_number_of_harmonics; i++)
        {
            double minus =  clip_min(log((double)getHarmonicDegree(i)), 0.);
            minus = -minus;
            
            double dot	= clip_min(log((double)getHarmonicDegree(i) + 1.), 0.);
            dot += minus;
            dot  = 1. / dot;
            
            for(unsigned int j = 0; j < NUMBEROFLINEARPOINTS; j++)
            {
                double weight = (1. - ((double)j / (double)(NUMBEROFLINEARPOINTS-1))) * weight_order + 1.;
                double scale = ((double)j / (double)(NUMBEROFLINEARPOINTS-1)) * weight_order;
                double new_weight = (minus + scale) * dot;
                new_weight = clip_minmax(new_weight, 0., 1.);
                m_wide_matrix[j * m_number_of_harmonics + i] = new_weight * weight;
            }
        }
        
		m_first_source      = 0;
        for(unsigned int i = 0; i < m_number_of_sources; i++)
        {
			setMute(i, 0);
            setAzimuth(i, 0.);
            setRadius(i, 1.);
        }
    }
    
    void Map::setAzimuth(const unsigned int index, const double azimuth)
    {
        assert(index < m_number_of_sources);
        m_azimuth[index] = wrap_twopi(azimuth);
        m_cosx[index]    = cos(m_azimuth[index]);
        m_sinx[index]    = sin(m_azimuth[index]);
    }
	
    void Map::setRadius(const unsigned int index, const double radius)
    {
        assert(index < m_number_of_sources);
        if(radius >= 1.)
        {
            m_gains[index] = 1. / (radius * radius);
            m_wide[index] = NUMBEROFLINEARPOINTS - 1;
        }
        else
        {
            m_gains[index] = 1.;
            m_wide[index] = clip_minmax(radius, 0., 1.) * (double)(NUMBEROFLINEARPOINTS - 1);
        }
    }
    
    void Map::setMute(const unsigned int index, const bool muted)
    {
        assert(index < m_number_of_sources);
        m_muted[index] = muted;
        m_first_source = -1;
        for(unsigned int i = 0; i < m_number_of_sources; i++)
        {
            if(!m_muted[i])
            {
                m_first_source = i;
                break;
            }
        }
    }
    
    void Map::process(const float* inputs, float* outputs)
    {
        int first = m_first_source;
        if(first > -1)
        {
            int index = m_wide[first] * m_number_of_harmonics;
            float cos_x = m_cosx[first];
            float sin_x = m_sinx[first];
            float tcos_x = cos_x;
            float sig = inputs[first] * m_gains[first] * m_wide_matrix[index];
            outputs[0] = sig;
            for(unsigned int i = 1; i < m_number_of_harmonics; i += 2)
            {
                outputs[i] = sig * sin_x * m_wide_matrix[index + i];
                outputs[i+1] = sig * cos_x * m_wide_matrix[index + i + 1];
                cos_x = tcos_x * m_cosx[first] - sin_x * m_sinx[first]; // cos(x + b) = cos(x) * cos(b) - sin(x) * sin(b)
                sin_x = tcos_x * m_sinx[first] + sin_x * m_cosx[first]; // sin(x + b) = cos(x) * sin(b) + sin(x) * cos(b)
                tcos_x = cos_x;
            }
            
            for(unsigned int i = first+1; i < m_number_of_sources; i++)
            {
                if (!m_muted[i])
                {
                    index = m_wide[i] * m_number_of_harmonics;
                    cos_x = m_cosx[i];
                    sin_x = m_sinx[i];
                    tcos_x = cos_x;
                    sig = inputs[i] * m_gains[i] * m_wide_matrix[index];
                    outputs[0] += sig;
                    for(unsigned int j = 1; j < m_number_of_harmonics; j += 2)
                    {
                        outputs[j] += sig * sin_x * m_wide_matrix[index + j];
                        outputs[j+1] += sig * cos_x * m_wide_matrix[index + j + 1];
                        cos_x = tcos_x * m_cosx[i] - sin_x * m_sinx[i]; // cos(x + b) = cos(x) * cos(b) - sin(x) * sin(b)
                        sin_x = tcos_x * m_sinx[i] + sin_x * m_cosx[i]; // sin(x + b) = cos(x) * sin(b) + sin(x) * cos(b)
                        tcos_x = cos_x;
                    }
                    
                }
            }
        }
        else
        {
            for(unsigned int i = 0; i < m_number_of_harmonics; i++)
                outputs[i] = 0.;
        }
    }
    
    void Map::process(const double* inputs, double* outputs)
    {
        int first = m_first_source;
		if(first > -1)
        {
            int index = m_wide[first] * m_number_of_harmonics;
            double cos_x = m_cosx[first];
            double sin_x = m_sinx[first];
            double tcos_x = cos_x;
            double sig = inputs[first] * m_gains[first] * m_wide_matrix[index];
            outputs[0] = sig;
            for(unsigned int i = 1; i < m_number_of_harmonics; i += 2)
            {
                outputs[i] = sig * sin_x * m_wide_matrix[index + i];
                outputs[i+1] = sig * cos_x * m_wide_matrix[index + i + 1];
                cos_x = tcos_x * m_cosx[first] - sin_x * m_sinx[first]; // cos(x + b) = cos(x) * cos(b) - sin(x) * sin(b)
                sin_x = tcos_x * m_sinx[first] + sin_x * m_cosx[first]; // sin(x + b) = cos(x) * sin(b) + sin(x) * cos(b)
                tcos_x = cos_x;
            }

            for(unsigned int i = first+1; i < m_number_of_sources; i++)
            {
                if (!m_muted[i])
                {
                    index = m_wide[i] * m_number_of_harmonics;
                    cos_x = m_cosx[i];
                    sin_x = m_sinx[i];
                    tcos_x = cos_x;
                    sig = inputs[i] * m_gains[i] * m_wide_matrix[index];
                    outputs[0] += sig;
                    for(unsigned int j = 1; j < m_number_of_harmonics; j += 2)
                    {
                        outputs[j] += sig * sin_x * m_wide_matrix[index + j];
                        outputs[j+1] += sig * cos_x * m_wide_matrix[index + j + 1];
                        cos_x = tcos_x * m_cosx[i] - sin_x * m_sinx[i]; // cos(x + b) = cos(x) * cos(b) - sin(x) * sin(b)
                        sin_x = tcos_x * m_sinx[i] + sin_x * m_cosx[i]; // sin(x + b) = cos(x) * sin(b) + sin(x) * cos(b)
                        tcos_x = cos_x;
                    }
                    
                }
            }
        }
        else
        {
            for(unsigned int i = 0; i < m_number_of_harmonics; i++)
                outputs[i] = 0.;
        }
    }
    
    Map::~Map()
    {
        delete [] m_gains;
        delete [] m_muted;
        delete [] m_azimuth;
        delete [] m_cosx;
        delete [] m_sinx;
        delete [] m_wide_matrix;
        delete [] m_wide;
    }
}

