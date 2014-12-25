/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_2D_DECODER_REGULAR
#define DEF_HOA_2D_DECODER_REGULAR

#include "Ambisonic_2D.h"
#include "Planewaves_2D.h"
#include "Encoder_2D.h"
#include "Rotate_2D.h"

namespace Hoa2D
{
    //! The ambisonic regular decoder.
    /** The regular decoder should be used to decode an ambisonic sound field for a set a channels at equal distances on a circle depending on a decomposition order. The number of channels must be at least the number of harmonics. Note that you can only change the offset of the channels.
     */
    class DecoderRegular : public Encoder<float>, public Planewaves
    {
    private:
        double          m_offset;
        double*         m_decoder_matrix_double;
        float*          m_decoder_matrix_float;
    public:
        
        //! The regular decoder constructor.
        /**	The regular decoder constructor allocates and initialize the member values to the decoding matrix depending of a decomposition order and a number of channels. The order must be at least 1 and the number of channels must be at least the number of harmonics.
         
            @param     order				The order
            @param     numberOfChannels     The number of channels.
         */
		DecoderRegular(unsigned long order, unsigned long numberOfChannels) noexcept :
        Encoder(order),
        Planewaves(numberOfChannels)
        {
            m_decoder_matrix_double     = new double[m_number_of_channels * m_number_of_harmonics];
            m_decoder_matrix_float      = new float[m_number_of_channels * m_number_of_harmonics];
            setChannelsOffset(0.);
        }
		
        //! The regular decoder destructor.
        /**	The regular decoder destructor free the memory.
         */
		~DecoderRegular()
        {
            delete [] m_decoder_matrix_double;
            delete [] m_decoder_matrix_float;
        }
        
        //! Set the offset of the channels.
		/**	Set the azimuth offset of the channels in radian.
            @param     offset		An azimuth value.
         */
		inline void setChannelsOffset(const double offset) noexcept
        {
            const double factor = 1. / (m_order_of_decomposition + 1.);
            m_offset = wrap_twopi(offset);
            for(unsigned long i = 0; i < m_number_of_channels; i++)
            {
                Encoder::setAzimuth(m_channels_azimuth[i] + m_offset);
                Encoder::process(factor, m_decoder_matrix_float+i*m_number_of_harmonics);
                //Encoder::process(factor, m_decoder_matrix_double+i*m_number_of_harmonics);
                m_decoder_matrix_float[i * m_number_of_harmonics] = factor * 0.5f;
                m_decoder_matrix_double[i * m_number_of_harmonics]= factor * 0.5f;
            }
        }
        
        //! Get the offset of the channels.
        /**	Retreive the azimuth offset of the channels in radian.
            @return    The offset of the channels.
         */
		inline double getChannelsOffset() const noexcept
        {
            return m_offset;
        }
				
        //! This method performs the regular decoding with single precision.
		/**	You should use this method for in-place or not-in-place processing and performs the regular decoding sample by sample. The inputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics and the outputs array contains the channels samples and the minimym size must be the number of channels.
            @param     inputs  The input array that contains the samples of the harmonics.
            @param     outputs The output array that contains samples destinated to channels.
         */
		inline void process(const float* input, float* output) const noexcept
        {
            cblas_sgemv(CblasRowMajor, CblasNoTrans, m_number_of_channels, m_number_of_harmonics, 1.f, m_decoder_matrix_float, m_number_of_harmonics, (float *)input, 1, 0.f, output, 1);
        }
		
        //! This method performs the regular decoding with double precision.
		/**	You should use this method for in-place or not-in-place processing and performs the regular decoding sample by sample. The inputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics and the outputs array contains the channels samples and the minimym size must be the number of channels.
            @param     inputs  The input array that contains the samples of the harmonics.
            @param     outputs The output array that contains the samples destinated to channels.
         */
		inline void process(const double* input, double* output) const noexcept
        {
            cblas_dgemv(CblasRowMajor, CblasNoTrans, m_number_of_channels, m_number_of_harmonics, 1.f, m_decoder_matrix_double, m_number_of_harmonics, (double *)input, 1, 0.f, output, 1);
        }
    };
}




#endif


