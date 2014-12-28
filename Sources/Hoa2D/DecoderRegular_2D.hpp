/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_2D_DECODER_REGULAR
#define DEF_HOA_2D_DECODER_REGULAR

#include "Ambisonic_2D.hpp"

namespace hoa
{
    //! The ambisonic regular decoder.
    /** The regular decoder should be used to decode an ambisonic sound field for a set a channels at equal distances on a circle depending on a decomposition order. The number of channels must be at least the number of harmonics. Note that you can only change the offset of the channels.
     */
    template <typename T> class DecoderRegular : public Encoder<HOA2D, T>, public hoa::Planewaves<T>
    {
    private:
        T*  m_matrix;
    public:
        
        //! The regular decoder constructor.
        /**	The regular decoder constructor allocates and initialize the member values to the decoding matrix depending of a decomposition order and a number of channels. The order must be at least 1 and the number of channels must be at least the number of harmonics.
         
            @param     order				The order
            @param     numberOfChannels     The number of channels.
         */
		DecoderRegular(unsigned long order, unsigned long numberOfChannels) noexcept :
        Encoder<HOA2D, T>(order),
        Planewaves<T>(numberOfChannels)
        {
            m_matrix = new T[Planewaves<T>::getNumberOfChannels() * Encoder<HOA2D, T>::getNumberOfHarmonics()];
            setChannelsOffset(0.);
        }
		
        //! The regular decoder destructor.
        /**	The regular decoder destructor free the memory.
         */
		~DecoderRegular()
        {
            delete [] m_matrix;
        }
        
        //! Set the offset of the channels.
		/**	Set the azimuth offset of the channels in radian.
            @param     offset		An azimuth value.
         */
		inline void setChannelsOffset(const T offset) noexcept
        {
            Planewaves<T>::setChannelsOffset(offset);
            const T factor = 1. / (Encoder<HOA2D, T>::getDecompositionOrder() + 1.);
            for(unsigned long i = 0; i < Planewaves<T>::getNumberOfChannels(); i++)
            {
                Encoder<HOA2D, T>::setAzimuth(Planewaves<T>::getChannelAzimuth(i) + Planewaves<T>::getChannelsOffset());
                Encoder<HOA2D, T>::process(&factor, m_matrix+i*Encoder<HOA2D, T>::getNumberOfHarmonics());
                m_matrix[i * Encoder<HOA2D, T>::getNumberOfHarmonics()] = factor * 0.5f;
            }
        }
        
        //! Set the azimuth of a channel.
        /** Set the azimuth of a channel. The azimuth is in radian between 0 and 2 π, O is the front of the soundfield and π is the back of the sound field. The maximum index must be the number of channel - 1.
         @param     index		The index of the channel.
         @param     azimuth		The azimuth.
         */
        inline void setChannelAzimuth(unsigned long index, T azimuth) noexcept override
        {
            Planewaves<T>::setChannelAzimuth(index, azimuth);
        }
				
        //! This method performs the regular decoding.
		/**	You should use this method for in-place or not-in-place processing and performs the regular decoding sample by sample. The inputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics and the outputs array contains the channels samples and the minimym size must be the number of channels.
            @param     inputs  The input array that contains the samples of the harmonics.
            @param     outputs The output array that contains samples destinated to channels.
         */
		inline void process(const T* inputs, T* outputs) const noexcept
        {
            matrix_vector_mul(Encoder<HOA2D, T>::m_number_of_harmonics, Planewaves<T>::m_number_of_channels, inputs, m_matrix, outputs);
        }
        
        void compute() noexcept
        {
            double  current_distance, minimum_distance;
            /*
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
                    //m_encoder->process(1., m_harmonics_vector);
                    
                    m_decoder_matrix[0] += (0.5 / (double)(m_order_of_decomposition + 1.));
                    for(unsigned int j = 1; j < m_number_of_harmonics; j++)
                    {
                        m_decoder_matrix[j] += (m_harmonics_vector[j] / (double)(m_order_of_decomposition + 1.));
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
                    //m_encoder->process(1., m_harmonics_vector);
                    
                    m_decoder_matrix[0] += (0.5 / (double)(m_order_of_decomposition + 1.));
                    m_decoder_matrix[m_number_of_harmonics] += (0.5 / (double)(m_order_of_decomposition + 1.));
                    
                    factor_index1 = fabs(cos(distance_radian(angle, m_channels_azimuth[0]) / HOA_PI * HOA_PI2));
                    factor_index2 = fabs(cos(distance_radian(angle, m_channels_azimuth[1]) / HOA_PI * HOA_PI2));
                    for(unsigned int j = 1; j < m_number_of_harmonics; j++)
                    {
                        m_decoder_matrix[j] += (m_harmonics_vector[j] / (double)(m_order_of_decomposition + 1.)) * factor_index1;
                        
                        m_decoder_matrix[m_number_of_harmonics + j] += (m_harmonics_vector[j] / (double)(m_order_of_decomposition + 1.)) * factor_index2;
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
                    //m_encoder->process(1., m_harmonics_vector);
                    
                    m_decoder_matrix[channel_index1 * m_number_of_harmonics] += (0.5 / (double)(m_order_of_decomposition + 1.)) * factor_index1;
                    m_decoder_matrix[channel_index2 * m_number_of_harmonics] += (0.5 / (double)(m_order_of_decomposition + 1.)) * factor_index2;
                    for(unsigned int j = 1; j < m_number_of_harmonics; j++)
                    {
                        m_decoder_matrix[channel_index1 * m_number_of_harmonics + j] += (m_harmonics_vector[j] / (double)(m_order_of_decomposition + 1.)) * factor_index1;
                        m_decoder_matrix[channel_index2 * m_number_of_harmonics + j] += (m_harmonics_vector[j] / (double)(m_order_of_decomposition + 1.)) * factor_index2;
                    }
                }
                
                for(unsigned int i = 0; i < m_number_of_channels * m_number_of_harmonics; i++)
                {
                    m_decoder_matrix_float[i] = m_decoder_matrix[i];
                }
            }
             */
        }
    };
}




#endif


