/*
// Copyright (c) 2012-2014 Eliott Paris & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include "Kits.h"

namespace Hoa2D
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Kit Sources //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    KitSources::PolarLines::PolarLines(unsigned int numberOfSources)
    {
        assert(numberOfSources > 0);
        m_number_of_sources = numberOfSources;
        
        m_values_old    = new float[m_number_of_sources * 2];
        m_values_new    = new float[m_number_of_sources * 2];
        m_values_step   = new float[m_number_of_sources * 2];
    }
    
    KitSources::PolarLines::~PolarLines()
    {
        delete [] m_values_old;
        delete [] m_values_new;
        delete [] m_values_step;
    }
    
    void KitSources::PolarLines::setRamp(unsigned int ramp)
    {
        m_ramp = clip_min(ramp, (long)1);
    }
    
    void KitSources::PolarLines::setRadius(unsigned int index, double radius)
    {
        assert(index < m_number_of_sources);
        m_values_new[index]  = radius;
        m_values_step[index] = (m_values_new[index] - m_values_old[index]) / (double)m_ramp;
        m_counter = 0;
    }
    
    void KitSources::PolarLines::setAzimuth(unsigned int index, double azimuth)
    {
        assert(index < m_number_of_sources);
        double new_value = wrap_twopi(azimuth);
        double old_value = wrap_twopi(m_values_old[index + m_number_of_sources]);
        double distance;
        
        
        if(old_value > new_value)
            distance = old_value - new_value;
        else
            distance = new_value - old_value;
        
        if(distance > HOA_PI && new_value > old_value)
            new_value -= HOA_2PI;
        else if(distance > HOA_PI)
            new_value += HOA_2PI;
        
        m_values_new[index + m_number_of_sources] = new_value;
        m_values_step[index + m_number_of_sources] = (new_value - old_value) / (double)m_ramp;
        m_counter = 0;
    }
    
    void KitSources::PolarLines::setRadiusDirect(unsigned int index, double radius)
    {
        assert(index < m_number_of_sources);
        m_values_old[index] = m_values_new[index] = radius;
        m_values_step[index] = 0.;
        m_counter = 0;
    }
    
    void KitSources::PolarLines::setAzimuthDirect(unsigned int index, double azimuth)
    {
        assert(index < m_number_of_sources);
        m_values_old[index + m_number_of_sources] = m_values_new[index + m_number_of_sources] = azimuth;
        m_values_step[index + m_number_of_sources] = 0.;
        m_counter = 0;
    }
    
    void KitSources::PolarLines::process(float* vector)
    {
        cblas_saxpy(m_number_of_sources * 2, 1., m_values_step, 1, m_values_old, 1);
        if(m_counter++ >= m_ramp)
        {
            cblas_scopy(m_number_of_sources * 2, m_values_new, 1, m_values_old, 1);
            memset(m_values_step, 0, sizeof(float) * m_number_of_sources * 2);
            m_counter    = 0;
        }
        cblas_scopy(m_number_of_sources * 2, m_values_old, 1, vector, 1);
    }
    
    KitSources::KitSources()
    {
        setMaximumRadius(20);
        m_order                 = 1;
        m_number_of_sources     = 1;
        m_number_of_channels    = 2;
        
        m_decoding_mode         = DecoderMulti::Irregular;
        m_optim_mode            = Optim::InPhase;
        m_offset                = 0.;
        
        m_lines                 = new PolarLines(m_number_of_sources);
        m_map                   = new Map(m_order, m_number_of_sources);
        m_optim                 = new Optim(m_order, m_optim_mode);
        m_decoder               = new DecoderMulti(m_order);
        m_meter                 = new Meter(m_number_of_channels);
		m_decoder->setDecodingMode(m_decoding_mode);
        m_decoder->setNumberOfChannels(m_number_of_channels);
        
		sourceNewPolar(1., 0.);
        m_lines->setRadiusDirect(0, 1.);
        m_lines->setAzimuthDirect(0, 0.);
        
        m_inputs_double     = new double[8192 * 64];
        m_outputs_double    = new double[8192 * 64];
        m_harmonics_double  = new double[8192 * 64];
        m_inputs_float      = new float[8192 * 64];
        m_outputs_float     = new float[8192 * 64];
        m_harmonics_float   = new float[8192 * 64];
        m_lines_vector      = new float[64];
        
        for(int i = 0; i < m_meter->getNumberOfChannels(); i++)
        {
            m_channels_azimuth[i] = m_meter->getChannelAzimuth(i);
            m_channels_azimuth_mapped[i] = m_meter->getChannelAzimuthMapped(i);
            m_channels_width[i] = m_meter->getChannelWidth(i);
        }
    }
    
    void KitSources::setOrder(unsigned int order)
    {
        m_order = clip_min(order, 1);
        setNumberOfChannels(m_number_of_channels);
    }
    
    void KitSources::setNumberOfSources(unsigned int numberOfSources)
    {
        m_number_of_sources = clip_min(numberOfSources, 1);
    }
    
    void KitSources::setNumberOfChannels(unsigned int numberOfChannels)
    {
        if(m_decoding_mode == DecoderMulti::Irregular)
            m_number_of_channels = clip_min(numberOfChannels, 2);
        else if(m_decoding_mode == DecoderMulti::Regular)
            m_number_of_channels = clip_min(numberOfChannels, m_order * 2 + 1);
        else
            m_number_of_channels = 2;
        
    }
    
    void KitSources::setDecodingMode(DecoderMulti::Mode mode)
    {
        m_decoding_mode = mode;
        setNumberOfChannels(m_number_of_channels);
    }
    
    void KitSources::setOptimMode(Optim::Mode mode)
    {
        m_optim_mode = mode;
        m_optim->setMode(mode);
    }
    
    void KitSources::setChannelsOffset(double offset)
    {
        m_offset = wrap_twopi(offset);
    }
    
    void KitSources::setSampleRate(unsigned int sampleRate)
    {
        m_sample_rate = sampleRate;
        m_decoder->setSampleRate(m_sample_rate);
    }

    void KitSources::setVectorSize(unsigned int vectorSize)
    {
        m_vector_size = vectorSize;
        m_decoder->setVectorSize(m_vector_size);
        m_lines->setRamp(4410);
        for(int i = 0; i < m_lines->getNumberOfSources(); i++)
        {
            m_lines->setRadiusDirect(i, sourceGetRadius(i));
            m_lines->setAzimuthDirect(i, sourceGetAzimuth(i));
        }
        
        for(int  i = 0 ; i < 64; i++)
        {
            m_harmonics_float_bin[i] = m_harmonics_float+i*m_vector_size;
        }
        m_outputs_float_bin[0] = m_outputs_float;
        m_outputs_float_bin[0] = m_outputs_float+m_vector_size;
    }
    
    void KitSources::process(const float** ins, float** outs)
	{
        int numins  = m_map->getNumberOfSources();
        int numouts = m_decoder->getNumberOfChannels();
        int nharmo  = m_map->getNumberOfHarmonics();
        
        for(int i = 0; i < numins; i++)
        {
            m_lines->setRadius(i, sourceGetRadius(i));
            m_lines->setAzimuth(i, sourceGetAzimuth(i));
        }
        for(int i = 0; i < numins; i++)
        {
            cblas_scopy(m_vector_size, ins[i], 1, m_inputs_float+i, numins);
            m_map->setMute(i, sourceGetMute(i));
        }
        
        if(m_decoder->getDecodingMode() == DecoderMulti::Regular)
        {
            for(int i = 0; i < m_vector_size; i++)
            {
                m_lines->process(m_lines_vector);
                for(int j = 0; j < numins; j++)
                    m_map->setRadius(j, m_lines_vector[j]);
                for(int j = 0; j < numins; j++)
                    m_map->setAzimuth(j, m_lines_vector[j+numins]);
                
                m_map->process(m_inputs_float + numins * i, m_harmonics_float + nharmo * i);
                m_optim->process(m_harmonics_float + nharmo * i, m_harmonics_float + nharmo * i);
                m_decoder->processRegular(m_harmonics_float + nharmo * i, m_outputs_float + numouts * i);
                m_meter->process(m_outputs_float + numouts * i);
            }
        }
        else if(m_decoder->getDecodingMode() == DecoderMulti::Irregular)
        {
            for(int i = 0; i < m_vector_size; i++)
            {
                m_lines->process(m_lines_vector);
                for(int j = 0; j < numins; j++)
                    m_map->setRadius(j, m_lines_vector[j]);
                for(int j = 0; j < numins; j++)
                    m_map->setAzimuth(j, m_lines_vector[j+numins]);
                 
                m_map->process(m_inputs_float + numins * i, m_harmonics_float + nharmo * i);
                m_optim->process(m_harmonics_float + nharmo * i, m_harmonics_float + nharmo * i);
                m_decoder->processIrregular(m_harmonics_float + nharmo * i, m_outputs_float + numouts * i);
                m_meter->process(m_outputs_float + numouts * i);
            }
        }
        else
        {
            for(int i = 0; i < m_vector_size; i++)
            {
                m_lines->process(m_lines_vector);
                for(int j = 0; j < numins; j++)
                    m_map->setRadius(j, m_lines_vector[j]);
                for(int j = 0; j < numins; j++)
                    m_map->setAzimuth(j, m_lines_vector[j+numins]);
            
                m_map->process(m_inputs_float + numins * i, m_harmonics_float + nharmo * i);
                m_optim->process(m_harmonics_float + nharmo * i, m_harmonics_float + nharmo * i);
                m_decoder->processBinaural(m_harmonics_float + nharmo * i, m_outputs_float + numouts * i);
                m_meter->process(m_outputs_float + numouts * i);
            }
        }
        for(int i = 0; i < numouts; i++)
        {
            cblas_scopy(m_vector_size, m_outputs_float+i, numouts, outs[i], 1);
        }
	}
    
    void KitSources::process(const double* input, double* output)
	{
		;
	}
    
    bool KitSources::applyChanges()
    {
        bool changed = 0;
    
        if(m_order != m_map->getDecompositionOrder() || m_map->getNumberOfSources() != m_number_of_sources)
        {
            delete m_map;
            m_map       = new Map(m_order, m_number_of_sources);
            changed = 1;
        }
        if(m_order != m_optim->getDecompositionOrder())
        {
            delete m_optim;
            m_optim     = new Optim(m_order, m_optim_mode);
            changed = 1;
        }
        if(m_optim_mode != m_optim->getMode())
        {
            m_optim->setMode(m_optim_mode);
        }
        if(m_order != m_decoder->getDecompositionOrder())
        {
            delete m_decoder;
            m_decoder   = new DecoderMulti(m_order);
            m_decoder->setDecodingMode(m_decoding_mode);
            m_decoder->setVectorSize(m_vector_size);
            m_decoder->setSampleRate(m_sample_rate);
            m_decoder->setNumberOfChannels(m_number_of_channels);
            m_number_of_channels = m_decoder->getNumberOfChannels();
            
            changed = 1;
        }
        if(m_decoding_mode != m_decoder->getDecodingMode())
        {
            m_decoder->setDecodingMode(m_decoding_mode);
            m_decoder->setVectorSize(m_vector_size);
            m_decoder->setSampleRate(m_sample_rate);
            m_decoder->setNumberOfChannels(m_number_of_channels);
            m_number_of_channels = m_decoder->getNumberOfChannels();
            changed = 1;
        }
        if(m_number_of_channels != m_decoder->getNumberOfChannels())
        {
            m_decoder->setNumberOfChannels(m_number_of_channels);
            m_number_of_channels = m_decoder->getNumberOfChannels();
            changed = 1;
        }
        if(m_number_of_channels != m_meter->getNumberOfChannels())
        {
            delete m_meter;
            m_meter = new Meter(m_number_of_channels);
            for(int i = 0; i < m_meter->getNumberOfChannels(); i++)
            {
                m_meter->setChannelAzimuth(i, m_decoder->getChannelAzimuth(i));
            }
            changed = 1;
        }
        if(m_lines->getNumberOfSources() != m_number_of_sources)
        {
            delete m_lines;
            m_lines         = new PolarLines(m_number_of_sources);
            changed = 1;
        }
        
        if(changed)
        {
            for(int i = 0; i < m_meter->getNumberOfChannels(); i++)
            {
                m_meter->setChannelAzimuth(i, m_decoder->getChannelAzimuth(i));
            }
            for(int i = 0; i < m_meter->getNumberOfChannels(); i++)
            {
                m_channels_azimuth[i] = m_meter->getChannelAzimuth(i);
                m_channels_azimuth_mapped[i] = m_meter->getChannelAzimuthMapped(i);
                m_channels_width[i] = m_meter->getChannelWidth(i);
            }
        }

        return changed;
    }
	
	KitSources::~KitSources()
	{
        delete m_map;
        delete m_optim;
        delete m_decoder;
        delete m_meter;
        delete m_lines;
        
        delete [] m_inputs_double;
        delete [] m_outputs_double;
        delete [] m_harmonics_double;
        delete [] m_inputs_float;
        delete [] m_outputs_float;
        delete [] m_harmonics_float;
        delete [] m_lines_vector;
	}
}

