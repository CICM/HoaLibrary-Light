/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include "Decoder_3D.h"

namespace Hoa3D
{
	DecoderRegular::DecoderRegular(unsigned int order, unsigned int numberOfChannels) : Ambisonic(order), Planewaves(numberOfChannels)
	{
        m_harmonics_vector          = new double[m_number_of_harmonics];
        m_decoder_matrix            = new double[m_number_of_channels * m_number_of_harmonics];
        m_decoder_matrix_float      = new float[m_number_of_channels * m_number_of_harmonics];
        m_encoder                   = new Encoder(m_order);
        setChannelsPosition(m_channels_azimuth, m_channels_elevation);
	}
	
	void DecoderRegular::setChannelPosition(unsigned int index, double azimuth, double elevation)
	{
        Planewaves::setChannelPosition(index, azimuth, elevation);
        
        m_encoder->setAzimuth(m_channels_rotated_azimuth[index]);
        m_encoder->setElevation(m_channels_rotated_elevation[index]);
        m_encoder->process(12.5 / (double)((m_order+1.)*(m_order+1.)), m_harmonics_vector);
        
        for(unsigned int j = 0; j < m_number_of_harmonics; j++)
        {
            m_decoder_matrix_float[index * m_number_of_harmonics + j] = m_decoder_matrix[index * m_number_of_harmonics + j] = m_harmonics_vector[j] * m_encoder->getNormalization(j) * m_encoder->getNormalization(j);
        }
	}
    
    void DecoderRegular::setChannelsPosition(double* azimuths, double* elevations)
	{
        for(unsigned int i = 0; i < m_number_of_channels; i++)
            setChannelPosition(i, azimuths[i], elevations[i]);
	}
    
    void DecoderRegular::setChannelsRotation(double axis_x, double axis_y, double axis_z)
    {
        Planewaves::setChannelsRotation(axis_x, axis_y, axis_z);
        setChannelsPosition(m_channels_azimuth, m_channels_elevation);
    }
	
	void DecoderRegular::process(const float* input, float* output)
	{
		cblas_sgemv(CblasRowMajor, CblasNoTrans, m_number_of_channels, m_number_of_harmonics, 1.f, m_decoder_matrix_float, m_number_of_harmonics, input, 1, 0.f, output, 1);
	}
	
	void DecoderRegular::process(const double* input, double* output)
	{
		cblas_dgemv(CblasRowMajor, CblasNoTrans, m_number_of_channels, m_number_of_harmonics, 1., m_decoder_matrix, m_number_of_harmonics, input, 1, 0., output, 1);
	}
	
	DecoderRegular::~DecoderRegular()
	{
		delete [] m_decoder_matrix;
        delete [] m_decoder_matrix_float;
        delete [] m_harmonics_vector;
        delete m_encoder;
	}
	
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Decoder Binaural //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    DecoderBinaural::DecoderBinaural(unsigned int order) : Ambisonic(order), Planewaves(2)
    {
        m_channels_azimuth[0] = HOA_PI2;
        m_channels_azimuth[1] = HOA_PI + HOA_PI2;
        if(m_order == 1)
            m_decoder = new DecoderRegular(m_order, 8);
        else if(m_order == 2)
            m_decoder = new DecoderRegular(m_order, 12);
        else if(m_order == 3)
            m_decoder = new DecoderRegular(m_order, 20);
        else
            m_decoder = new DecoderRegular(m_order, (m_order + 1) * (m_order + 1));
        
        for(int i = 0; i < m_decoder->getNumberOfChannels(); i++)
        {
            double angle = m_decoder->getChannelAzimuth(i);
            m_filters_left.push_back(BinauralFilter(angle, m_decoder->getChannelElevation(i)));
            m_filters_right.push_back(BinauralFilter(-angle, m_decoder->getChannelElevation(i)));
        }
        m_pinna_size = Small;
        m_outputs_double    = new double[m_decoder->getNumberOfChannels()];
        m_outputs_float     = new float[m_decoder->getNumberOfChannels()];
    }
    
    void DecoderBinaural::setChannelsRotation(double axis_x, double axis_y, double axis_z)
    {
        Planewaves::setChannelsRotation(axis_x, axis_y, axis_z);
    }
    void DecoderBinaural::setPinnaSize(PinnaSize pinnaSize)
    {
        m_pinna_size = pinnaSize;
    }
    
    void DecoderBinaural::setSampleRate(double sampleRate)
    {
        if (m_sampleRate != sampleRate)
        {
            m_sampleRate = sampleRate;
            
            m_filters_left.clear();
            m_filters_right.clear();
            
            for(int i = 0; i < m_decoder->getNumberOfChannels(); i++)
            {
                double angle = m_decoder->getChannelAzimuth(i);
                m_filters_left.push_back(BinauralFilter(angle, m_decoder->getChannelElevation(i)));
                m_filters_right.push_back(BinauralFilter(-angle, m_decoder->getChannelElevation(i)));
            }
        }
    }
    
    void DecoderBinaural::process(const float* inputs, float* outputs)
	{
        outputs[0] = 0.f;
        outputs[1] = 0.f;
        m_decoder->process(inputs, m_outputs_float);
        for(int i = 0; i < m_decoder->getNumberOfChannels(); i++)
        {
            outputs[0] += m_filters_left[i].process(m_outputs_float[i]);
            outputs[1] += m_filters_right[i].process(m_outputs_float[i]);
        }
    }
    
    void DecoderBinaural::process(const double* inputs, double* outputs)
	{
        outputs[0] = 0.f;
        outputs[1] = 0.f;
        m_decoder->process(inputs, m_outputs_double);
        for(int i = 0; i < m_decoder->getNumberOfChannels(); i++)
        {
            outputs[0] += m_filters_left[i].process(m_outputs_double[i]);
            outputs[1] += m_filters_right[i].process(m_outputs_double[i]);
        }
    }
    
	DecoderBinaural::~DecoderBinaural()
	{
        delete m_decoder;
        m_filters_left.clear();
        m_filters_right.clear();
        delete [] m_outputs_double;
        delete [] m_outputs_float;
	}
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Decoder Multi //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    DecoderMulti::DecoderMulti(unsigned int order) : Ambisonic(order)
    {
        m_mode = Regular;
        m_decoder_regular   = new DecoderRegular(m_order, (m_order + 1) * (m_order + 1));
        m_decoder_binaural  = new DecoderBinaural(m_order);
    }
	
	DecoderMulti::DecoderMulti(unsigned int order, unsigned int numberOfChannels) : Ambisonic(order)
    {
        m_mode = Regular;
        m_decoder_regular   = new DecoderRegular(m_order, numberOfChannels);
        m_decoder_binaural  = new DecoderBinaural(m_order);
    }
    
    void DecoderMulti::setDecodingMode(Mode mode)
    {
        m_mode = mode;
    }
    
    void DecoderMulti::setNumberOfChannels(unsigned int numberOfChannels)
    {
        if(numberOfChannels != getNumberOfChannels())
        {
            if(m_mode == Regular)
            {
                delete m_decoder_regular;
                m_decoder_regular = new DecoderRegular(m_order, numberOfChannels);
            }
        }
    }
    
    void DecoderMulti::setChannelPosition(unsigned int index, double azimuth, double elevation)
    {
        if(m_mode == Regular)
        {
            m_decoder_regular->setChannelPosition(index, azimuth, elevation);
        }
    }
    
    void DecoderMulti::setChannelsPosition(double* azimuths, double* elevations)
    {
        if(m_mode == Regular)
        {
            m_decoder_regular->setChannelsPosition(azimuths, elevations);
        }
    }
    
    void DecoderMulti::setChannelsRotation(double axis_x, double axis_y, double axis_z)
    {
        m_decoder_regular->setChannelsRotation(axis_x, axis_y, axis_z);
        m_decoder_binaural->setChannelsRotation(axis_x, axis_y, axis_z);
    }
    
    void DecoderMulti::setPinnaSize(DecoderBinaural::PinnaSize pinnaSize)
    {
        m_decoder_binaural->setPinnaSize(pinnaSize);
    }
    
    void DecoderMulti::setSampleRate(double sampleRate)
    {
        m_decoder_binaural->setSampleRate(sampleRate);
        m_sample_rate = sampleRate;
    }
    
	DecoderMulti::~DecoderMulti()
	{
        delete m_decoder_regular;
        delete m_decoder_binaural;
	}

}


