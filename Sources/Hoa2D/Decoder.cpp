/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include "Decoder.h"

namespace Hoa2D
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Decoder Regular //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    DecoderRegular::DecoderRegular(unsigned int order, unsigned int numberOfChannels) : Ambisonic(order), Planewaves(numberOfChannels)
    {
        assert(numberOfChannels >= m_number_of_harmonics);

        m_harmonics_vector          = new double[m_number_of_harmonics];
        m_decoder_matrix_double     = new double[m_number_of_channels * m_number_of_harmonics];
        m_decoder_matrix_float      = new float[m_number_of_channels * m_number_of_harmonics];
        m_encoder                   = new Encoder(m_order);
        setChannelsOffset(0.);
    }

    void DecoderRegular::setChannelsOffset(double offset)
	{
        m_offset = wrap_twopi(offset);
        for(unsigned int i = 0; i < m_number_of_channels; i++)
        {
            m_encoder->setAzimuth(m_channels_azimuth[i] + m_offset);
            m_encoder->process(1., m_harmonics_vector);

            m_decoder_matrix_float[i * m_number_of_harmonics] = m_decoder_matrix_double[i * m_number_of_harmonics] = 0.5 / (double)(m_order + 1.);
            for(unsigned int j = 1; j < m_number_of_harmonics; j++)
            {
                m_decoder_matrix_float[i * m_number_of_harmonics + j] = m_decoder_matrix_double[i * m_number_of_harmonics + j] = m_harmonics_vector[j] / (double)(m_order + 1.);
            }
        }
	}

    void DecoderRegular::process(const float* input, float* output)
	{
		cblas_sgemv(CblasRowMajor, CblasNoTrans, m_number_of_channels, m_number_of_harmonics, 1.f, m_decoder_matrix_float, m_number_of_harmonics, (float *)input, 1, 0.f, output, 1);
	}

	void DecoderRegular::process(const double* input, double* output)
	{
		cblas_dgemv(CblasRowMajor, CblasNoTrans, m_number_of_channels, m_number_of_harmonics, 1.f, m_decoder_matrix_double, m_number_of_harmonics, (double *)input, 1, 0.f, output, 1);
	}

	DecoderRegular::~DecoderRegular()
	{
		delete [] m_decoder_matrix_double;
        delete [] m_decoder_matrix_float;
        delete [] m_harmonics_vector;
        delete m_encoder;
	}

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Decoder Irregular //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    DecoderIrregular::DecoderIrregular(unsigned int order, unsigned int numberOfChannels) : Ambisonic(order), Planewaves(numberOfChannels)
    {
        m_harmonics_vector          = new double[m_number_of_harmonics];
        m_decoder_matrix            = new double[m_number_of_channels * m_number_of_harmonics];
        m_decoder_matrix_float      = new float[m_number_of_channels * m_number_of_harmonics];
        m_encoder                   = new Encoder(m_order);
        m_nearest_channel[0]        = NULL;
        m_nearest_channel[1]        = NULL;
        
        m_offset = 0;
        setChannelAzimuth(0, 0.);
    }

    void DecoderIrregular::setChannelsOffset(double offset)
	{
        m_offset = wrap_twopi(offset);
        setChannelAzimuth(0, m_channels_azimuth[0]);
    }

    void DecoderIrregular::setChannelsAzimuth(double* azimuths)
    {
        Planewaves::setChannelsAzimuth(azimuths);
        setChannelAzimuth(0, m_channels_azimuth[0]);
    }

    void DecoderIrregular::setChannelAzimuth(unsigned int index, double azimuth)
    {
        double  current_distance, minimum_distance;

        Planewaves::setChannelAzimuth(index, azimuth);

        // Get the minimum distance between the channels
        minimum_distance    = HOA_2PI + 1;
        current_distance    = distance_radian(m_channels_azimuth[0], m_channels_azimuth[m_number_of_channels-1]);
        if(current_distance < minimum_distance)
            minimum_distance    = current_distance;
        for(unsigned int i = 1; i < m_number_of_channels; i++)
        {
            current_distance  = distance_radian(m_channels_azimuth[i], m_channels_azimuth[i-1]);
            if(current_distance < minimum_distance)
                minimum_distance = current_distance;
        }

        // Get the optimal number of virtual channels
        // Always prefer the number of harmonics + 1
        if(minimum_distance > 0)
            m_number_of_virtual_channels = (HOA_2PI / minimum_distance);
        else
            m_number_of_virtual_channels = m_number_of_harmonics + 1;
        if(m_number_of_virtual_channels < m_number_of_harmonics + 1)
        {
            m_number_of_virtual_channels = m_number_of_harmonics + 1;
        }
        
        if(m_nearest_channel[0] && m_nearest_channel[1])
        {
            delete [] m_nearest_channel[0];
            delete [] m_nearest_channel[1];
        }
        m_nearest_channel[0] = new unsigned int[m_number_of_virtual_channels];
        m_nearest_channel[1] = new unsigned int[m_number_of_virtual_channels];
        
        for(unsigned int i = 0; i < m_number_of_channels * m_number_of_harmonics; i++)
        {
            m_decoder_matrix[i] = 0.;
            m_decoder_matrix_float[i] = 0.;
        }

        if(m_number_of_channels == 1)
        {
            for(unsigned int i = 0; i < m_number_of_virtual_channels; i++)
            {
                double angle = (double)i / (double)m_number_of_virtual_channels * HOA_2PI;
                m_encoder->setAzimuth(angle + m_offset);
                m_encoder->process(1., m_harmonics_vector);

                m_decoder_matrix[0] += (0.5 / (double)(m_order + 1.));
                for(unsigned int j = 1; j < m_number_of_harmonics; j++)
                {
                    m_decoder_matrix[j] += (m_harmonics_vector[j] / (double)(m_order + 1.));
                }
            }
            for(unsigned int i = 0; i < m_number_of_harmonics; i++)
            {
                m_decoder_matrix_float[i] = m_decoder_matrix[i];
            }
        }
        else if(m_number_of_channels == 2)
        {
            for(unsigned int i = 0; i < m_number_of_virtual_channels; i++)
            {
                double factor_index1 = 0, factor_index2 = 0;
                double angle = (double)i / (double)m_number_of_virtual_channels * HOA_2PI;
                m_encoder->setAzimuth(angle + m_offset);
                m_encoder->process(1., m_harmonics_vector);

                m_decoder_matrix[0] += (0.5 / (double)(m_order + 1.));
                m_decoder_matrix[m_number_of_harmonics] += (0.5 / (double)(m_order + 1.));

                factor_index1 = fabs(cos(distance_radian(angle, m_channels_azimuth[0]) / HOA_PI * HOA_PI2));
                factor_index2 = fabs(cos(distance_radian(angle, m_channels_azimuth[1]) / HOA_PI * HOA_PI2));
                for(unsigned int j = 1; j < m_number_of_harmonics; j++)
                {
                    m_decoder_matrix[j] += (m_harmonics_vector[j] / (double)(m_order + 1.)) * factor_index1;

                    m_decoder_matrix[m_number_of_harmonics + j] += (m_harmonics_vector[j] / (double)(m_order + 1.)) * factor_index2;
                }
            }


            for(unsigned int i = 0; i < m_number_of_channels * m_number_of_harmonics; i++)
            {
                m_decoder_matrix_float[i] = m_decoder_matrix[i];
            }
        }
        else
        {
            // Get the nearest channels
            for(unsigned int i = 0; i < m_number_of_virtual_channels; i++)
            {
                long   channel_index1 = 0, channel_index2 = 0;
                double distance1 = HOA_2PI, distance2 = HOA_2PI;
                double angle = (double)i / (double)m_number_of_virtual_channels * HOA_2PI;
                double factor_index1 = 0, factor_index2 = 0;
                
                for(unsigned int j = 0; j < m_number_of_channels; j++)
                {
                    if(radianClosestDistance(m_channels_azimuth[j], angle) < distance1)
                    {
                        distance1 = radianClosestDistance(m_channels_azimuth[j], angle);
                        channel_index1 = j;
                    }
                }
                
                for(unsigned int j = 0; j < m_number_of_channels; j++)
                {
                    if(radianClosestDistance(m_channels_azimuth[j], angle) < distance2 && j != channel_index1)
                    {
                        distance2 = radianClosestDistance(m_channels_azimuth[j], angle);
                        channel_index2 = j;
                    }
                }
                
                if(fabs(distance1 - distance2) < HOA_PI / (double)m_number_of_virtual_channels)
                {
                    double angle1 = m_channels_azimuth[channel_index1], angle2 = m_channels_azimuth[channel_index2];
                    double distance_index1 = radianClosestDistance(angle, angle1);
                    double distance_index2 = radianClosestDistance(angle, angle2);
                    double distance_ratio = distance_index1 + distance_index2;
                    factor_index1   = cos(distance_index1 / (distance_ratio) * HOA_PI2);
                    factor_index2   = cos(distance_index2 / (distance_ratio) * HOA_PI2);
                }
                else
                {
                    factor_index1   = 1;
                    factor_index2   = 0;
                }
                
                // Get the harmonics coefficients for the virtual channel
                m_encoder->setAzimuth(angle + m_offset);
                m_encoder->process(1., m_harmonics_vector);
                
                m_decoder_matrix[channel_index1 * m_number_of_harmonics] += (0.5 / (double)(m_order + 1.)) * factor_index1;
                m_decoder_matrix[channel_index2 * m_number_of_harmonics] += (0.5 / (double)(m_order + 1.)) * factor_index2;
                for(unsigned int j = 1; j < m_number_of_harmonics; j++)
                {
                    m_decoder_matrix[channel_index1 * m_number_of_harmonics + j] += (m_harmonics_vector[j] / (double)(m_order + 1.)) * factor_index1;
                    m_decoder_matrix[channel_index2 * m_number_of_harmonics + j] += (m_harmonics_vector[j] / (double)(m_order + 1.)) * factor_index2;
                }
            }

            for(unsigned int i = 0; i < m_number_of_channels * m_number_of_harmonics; i++)
            {
                m_decoder_matrix_float[i] = m_decoder_matrix[i];
            }
        }
    }

    void DecoderIrregular::process(const float* input, float* output)
	{
		cblas_sgemv(CblasRowMajor, CblasNoTrans, m_number_of_channels, m_number_of_harmonics, 1.f, m_decoder_matrix_float, m_number_of_harmonics, input, 1, 0.f, output, 1);
	}

	void DecoderIrregular::process(const double* input, double* output)
	{
		cblas_dgemv(CblasRowMajor, CblasNoTrans, m_number_of_channels, m_number_of_harmonics, 1.f, m_decoder_matrix, m_number_of_harmonics, input, 1, 0.f, output, 1);
	}

	DecoderIrregular::~DecoderIrregular()
	{
		delete [] m_decoder_matrix;
        delete [] m_decoder_matrix_float;
        delete [] m_harmonics_vector;
        delete [] m_nearest_channel[0];
        delete [] m_nearest_channel[1];
        delete m_encoder;
	}

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Decoder Binaural //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    DecoderBinaural::DecoderBinaural(unsigned int order) : Ambisonic(order), Planewaves(2)
    {
        m_channels_azimuth[0] = HOA_PI2;
        m_channels_azimuth[1] = HOA_PI + HOA_PI2;
        m_decoder = new DecoderRegular(m_order, m_order * 2 + 2);
        for(int i = 0; i < m_decoder->getNumberOfChannels(); i++)
        {
            double angle = m_decoder->getChannelAzimuth(i);
            m_filters_left.push_back(BinauralFilter(angle, 0));
            m_filters_right.push_back(BinauralFilter(-angle, 0));
        }
        m_pinna_size = Small;
        m_outputs_double    = new double[m_decoder->getNumberOfChannels()];
        m_outputs_float     = new float[m_decoder->getNumberOfChannels()];
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
                m_filters_left.push_back(BinauralFilter(angle, 0));
                m_filters_right.push_back(BinauralFilter(-angle, 0));
            }
        }
        
        /*
        for(int i = 0; i < m_decoder->getNumberOfChannels(); i++)
        {
            m_filters_left[i].setSampleRate(sampleRate);
            m_filters_right[i].setSampleRate(sampleRate);
        }
        */
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

    DecoderMulti::DecoderMulti(unsigned int order) : Ambisonic(order), Planewaves(order * 2 + 2)
    {
        m_mode = Regular;
        m_decoder_regular   = new DecoderRegular(m_order, m_order * 2 + 2);
        m_decoder_irregular = new DecoderIrregular(m_order, m_order * 2 + 2);
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
            if(m_mode == Regular && numberOfChannels >= m_decoder_regular->getNumberOfHarmonics())
            {
                delete m_decoder_regular;
                m_decoder_regular = new DecoderRegular(m_order, numberOfChannels);
            }
            else if(m_mode == Irregular)
            {
                delete m_decoder_irregular;
                m_decoder_irregular = new DecoderIrregular(m_order, numberOfChannels);
            }
        }
    }

    void DecoderMulti::setChannelsOffset(double offset)
	{
        if(m_mode == Regular)
        {
            m_decoder_regular->setChannelsOffset(offset);
        }
        else if(m_mode == Irregular)
        {
            m_decoder_irregular->setChannelsOffset(offset);
        }
	}

    void DecoderMulti::setChannelAzimuth(unsigned int index, double azimuth)
    {
        if(m_mode == Irregular)
        {
            m_decoder_irregular->setChannelAzimuth(index, azimuth);
        }
    }

    void DecoderMulti::setChannelsAzimuth(double* azimuths)
    {
        if(m_mode == Irregular)
        {
            m_decoder_irregular->setChannelsAzimuth(azimuths);
        }
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
        delete m_decoder_irregular;
        delete m_decoder_binaural;
	}
}

