/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef __DEF_HOA_3D_DECODER__
#define __DEF_HOA_3D_DECODER__

#include "Ambisonic_3D.h"
#include "Planewaves_3D.h"
#include "Encoder_3D.h"

namespace Hoa3D
{
	//! The ambisonic decoder.
    /** The decoder should be used to decode a signal encoded in the spherical harmonics domain depending on a decomposition order and a number of channels.
     */
	class DecoderRegular : public Ambisonic, public Planewaves
	{
		
	private:
		double*         m_decoder_matrix;
        float*          m_decoder_matrix_float;
		double*         m_harmonics_vector;
        Encoder*        m_encoder;        
	public:
        
		/**	The decoder constructor.
         @param     order					The order, must be at least 1.
		 @param     numberOfChannels	The number of channels, must be at least (order + 1)^2.
		 @param     shape					Is a sphere or a half sphere.
         */
		DecoderRegular(unsigned int order, unsigned int numberOfChannels);
		
        /**	The decoder destructor.
         */
		~DecoderRegular();
        
		/**	Set channel position.
		 @param     index		The index of the channel.
		 @param     azimuth		An azimuth value. In radian, between 0 and 2π.
		 @param     elevation	An elevation value. In radian, between 0 and 2π.
         */
		void	setChannelPosition(unsigned int index, double azimuth, double elevation);
        
        //! Set the position of the channels.
        /** Set the position of the channels with polar coordinates. The azimtuh is in radian between 0 and 2 Pi, O is the front of the soundfield and Pi is the back of the sound field. The elevation is in radian between -1/2 Pi and 1/2 Pi, -1/2 Pi the the bottom of the sound field, 0 is the center of the sound field and 1/2 Pi is the top of the sound field. The maximum index must be the number of channels - 1.
         
         @param     azimuths		The azimuths.
         @param     elevations	The elevations.
         */
		void setChannelsPosition(double* azimuths, double* elevations);
        
        //! Set the rotation of the channels.
		/**	Set the angles in radian of the rotation of the channels around the axes x, y and z.
         
         @param     axis_x	The angle of rotation around the x axe.
         @param     axis_y	The angle of rotation around the y axe.
         @param     axis_z	The angle of rotation around the z axe.
         */
		void setChannelsRotation(double axis_x, double axis_y, double axis_z);
		
        /**	This method performs the decoding with single precision.
         @param     input	The inputs array.
         @param     outputs The output array that contains samples destinated to channels.
         */
		void process(const float* input, float* output);
		
		/**	This method performs the decoding with double precision.
         @param     input	The inputs array.
         @param     outputs The output array that contains samples destinated to channels.
         */
		void process(const double* input, double* output);
	};
    
    const float* get_mit_hrtf_3D(long samplerate, double azimuth, long elevation, bool large);
    
    //! The ambisonic binaural decoder.
    /** The binaural decoder should be used to decode an ambisonic sound field for headphones.
     */
    class DecoderBinaural : public Ambisonic, public Planewaves
    {
    public:
        
        enum PinnaSize
        {
            Small       = 0,	/**< Small Pinna Size  */
            Large       = 1,	/**< Large Pinna Size */
        };
        
    private:
        PinnaSize       m_pinna_size;
        double*         m_outputs_double;
        float*          m_outputs_float;
        DecoderRegular* m_decoder;
        double          m_sampleRate;
        std::vector<BinauralFilter> m_filters_left;
        std::vector<BinauralFilter> m_filters_right;
    public:
        
        //! The binaural decoder constructor.
        /**	The binaural decoder constructor allocates and initialize the member values to the decoding matrix depending of a decomposition order and a number of channels. The order and the number of channels must be at least 1 and the maximum order is 35. It is essential to set the sample rate and the vector size to load the impulse response and to be able to use the binaural decoding. The binaural process is optimized for block processing. The HRTF are from the MIT database.
         
         @param     order				The order
         */
		DecoderBinaural(unsigned int order);
		
        //! The binaural decoder destructor.
        /**	The binaural decoder destructor free the memory.
         */
		~DecoderBinaural();
        
        //! Set the rotation of the channels.
		/**	Set the angles in radian of the rotation of the channels around the axes x, y and z.
         
         @param     axis_x	The angle of rotation around the x axe.
         @param     axis_y	The angle of rotation around the y axe.
         @param     axis_z	The angle of rotation around the z axe.
         */
		void setChannelsRotation(double axis_x, double axis_y, double axis_z);
        
        //! Set the sample rate.
        /** Set the sample rate. The sample will change the impulse responses size and their sizes increase with it. The valid sample rate are 44100, 48000, 88200 and 9600. Setting the sample rate will load the impulse responses, it is essential to define it before the digital signal processing.
         
         @param     sampleRate		The sample rate.
         
         @see    setVectorSize
         */
        void setSampleRate(double sampleRate);
        
        //! Set the pinna size.
        /** Set the pinna size used to compute the HRTF. Setting the pinna size will re-allocate the vector to compute the binaural decoding.
         
         @param     pinnaSize		The pinna size.
         
         */
        void setPinnaSize(PinnaSize pinnaSize);
        
        //! Retrieve the pinna size.
        /** Retrieve the current size of the pinna.
         
         @return    The function returns the pinna size used to compute the HRTF.
         */
		inline PinnaSize getPinnaSize() const
        {
            return m_pinna_size;
        };
        
        //! Retrieve a name for a channel.
        /** Retrieve a name for a channel in a std::string format that will be "Headphone Left" or "Headphone Right".
         
         @param     index	The index of a channel.
         @return    The method returns a name for the channel.
         
         */
		inline std::string getChannelName(unsigned int index)
        {
            assert(index < 2);
            if(index == 0)
                return "Headphone Left";
            else
                return "Headphone Right";
        };
        
        //! This method performs the binaural decoding with single precision.
		/**	You should use this method for not-in-place processing and performs the binaural decoding sample by sample. The inputs array contains the spherical harmonics samples : inputs[number of harmonics] and the outputs array contains the headphones samples : outputs[2].
         
         @param     inputs	The input samples.
         @param     outputs  The output array that contains samples destinated to channels.
         */
		void process(const float* inputs, float* outputs);
		
        //! This method performs the binaural decoding with double precision.
		/**	You should use this method for not-in-place processing and performs the binaural decoding sample by sample. The inputs array contains the spherical harmonics samples : inputs[number of harmonics] and the outputs array contains the headphones samples : outputs[2].
         
         @param     input    The input samples.
         @param     outputs  The output array that contains samples destinated to channels.
         */
		void process(const double* inputs, double* outputs);
    };
    
    //! The ambisonic multi-decoder.
    /** The multi-decoder is a class that facilitates the use of the three decoder : regular, irregular and binaural.
     */
    class DecoderMulti : public Ambisonic
    {
    public:
        
        enum Mode
        {
            Regular     = 0,	/**< Regular Decoding   */
            Binaural    = 1     /**< Binaural Decoding  */
        };
        
    private:
        DecoderRegular*     m_decoder_regular;
        DecoderBinaural*    m_decoder_binaural;
        Mode                m_mode;
        double              m_sample_rate;
        
    public:
        
        //! The multi-decoder constructor.
        /**	The multi-decoder constructor allocates and initialize the three decoder. The default decoder will be the regular one with 2 * order + 2 number of channels.
         
         @param     order				The order
         */
        DecoderMulti(unsigned int order);
		
		//! The multi-decoder constructor.
        /**	The multi-decoder constructor allocates and initialize the three decoder.
         
         @param     order				The order
		 @param     numberOfChannels	The number of channels
         */
		DecoderMulti(unsigned int order, unsigned int numberOfChannels);
        
        //! The multi-decoder destructor.
        /**	The multi-decoder destructor free the memory.
         */
        ~DecoderMulti();
        
        //! Set the decoding mode.
        /**	Set the decoding mode. It will allocate the right decoder and initialize the class.
         
         @param     mode		The decoding mode.
         */
        void setDecodingMode(Mode mode);
        
        //! Retrieve the decoding mode.
        /** Retrieve the current decoding mode of the multi-decoder.
         
         @return    The decoding mode;
         */
        inline Mode getDecodingMode() const
        {
            return m_mode;
        };
        
        //! Set the number of channels for the regular or irregular decoding.
        /**	Set the number of channels for the regular or irregular decoding.
         
         @param     numberOfChannels		The number of channels.
         */
        void setNumberOfChannels(unsigned int numberOfChannels);
        
        //! Retrieve the number of channels.
        /** Retrieve the number of channels of the planewave class.
         
         @return The number of channels.
         */
        inline unsigned int getNumberOfChannels() const
        {
            if(m_mode == Regular)
                return m_decoder_regular->getNumberOfChannels();
            else
                return m_decoder_binaural->getNumberOfChannels();
        }
        
        /**	Set channel position.
         @param     index		The index of the channel.
         @param     azimuth		An azimuth value. In radian, between 0 and 2π.
         @param     elevation	An elevation value. In radian, between 0 and 2π.
         */
        void setChannelPosition(unsigned int index, double azimuth, double elevation);
        
        //! Set the position of the channels.
        /** Set the position of the channels with polar coordinates. The azimtuh is in radian between 0 and 2 Pi, O is the front of the soundfield and Pi is the back of the sound field. The elevation is in radian between -1/2 Pi and 1/2 Pi, -1/2 Pi the the bottom of the sound field, 0 is the center of the sound field and 1/2 Pi is the top of the sound field. The maximum index must be the number of channels - 1.
         
         @param     azimuths		The azimuths.
         @param     elevations	The elevations.
         */
		void setChannelsPosition(double* azimuths, double* elevations);
        
        //! Set the rotation of the channels.
		/**	Set the angles in radian of the rotation of the channels around the axes x, y and z.
         
         @param     axis_x	The angle of rotation around the x axe.
         @param     axis_y	The angle of rotation around the y axe.
         @param     axis_z	The angle of rotation around the z axe.
         */
		void setChannelsRotation(double axis_x, double axis_y, double axis_z);
        
        //! Get the x rotation of the channels.
        /**	Retreive the x rotation of the channels in radian.
         
         @return    The x rotation of the channels.
         */
		double getChannelsRotationX() const
        {
            if(m_mode == Regular)
                return m_decoder_regular->getChannelsRotationX();
            else
                return 0;
        }
        
        //! Get the y rotation of the channels.
        /**	Retreive the x rotation of the channels in radian.
         
         @return    The y rotation of the channels.
         */
		double getChannelsRotationY() const
        {
            if(m_mode == Regular)
                return m_decoder_regular->getChannelsRotationY();
            else
                return 0;
        }
        
        //! Get the z rotation of the channels.
        /**	Retreive the z rotation of the channels in radian.
         
         @return    The z rotation of the channels.
         */
		double getChannelsRotationZ() const
        {
            if(m_mode == Regular)
                return m_decoder_regular->getChannelsRotationZ();
            else
                return 0;
        }
        
        //! Set the sample rate.
        /** Set the sample rate. The sample will change the impulse responses size and their sizes increase with it. The valid sample rate are 44100, 48000, 88200 and 9600. Setting the sample rate will load the impulse responses, it is essential to define it before the digital signal processing.
         
         @param     sampleRate		The sample rate.
         
         @see    setVectorSize
         */
        void setSampleRate(double sampleRate);
        
        //! Set the pinna size.
        /** Set the pinna size of the binaural decoding.
         
         @param     pinnaSize		The pinna size.
         
         */
        void setPinnaSize(DecoderBinaural::PinnaSize pinnaSize);

        
        //! Retrieve if the pinna size of the binaural decoder.
        /** Retrieve if the pinna size of the binaural decoder.
         
         @return    The function returns the pinna size of the binaural decoder.
         */
        inline DecoderBinaural::PinnaSize getPinnaSize() const
        {
            return m_decoder_binaural->getPinnaSize();
        };
        
        //! Retrieve the azimuth of a channel.
        /** Retrieve the azimuth of a channel. The azimuth of the channel is in radian, 0 radian is at the front of the soundfield and Pi is at the back of the sound field. The maximum index must be the number of channels - 1.
         
         @param      index   The index of the channel.
         @return     The azimuth of the channel if the channel exists, otherwise the function generates an error.
         
         @see getChannelAbscissa
         @see getChannelOrdinate
         @see getChannelName
         */
        inline double getChannelAzimuth(unsigned int index) const
        {
            if(m_mode == Regular)
                return m_decoder_regular->getChannelAzimuth(index);
            else
                return m_decoder_binaural->getChannelAzimuth(index);
        }
        
        //! Retrieve the elevation of a channel.
        /** Retrieve the elevation of a channel. The elevation is in radian between -1/2 Pi and 1/2 Pi, -1/2 Pi the the bottom of the sound field, 0 is the center of the sound field and 1/2 Pi is the top of the sound field. The maximum index must be the number of channels - 1.
         
         @param      index   The index of the channel.
         @return     The elevation of the channel.
         @see getChannelAzimuth
         */
        inline double getChannelElevation(unsigned int index) const
        {
            if(m_mode == Regular)
                return m_decoder_regular->getChannelElevation(index);
            else
                return m_decoder_binaural->getChannelElevation(index);
        }
        
        //! Retrieve the abscissa of a channel.
        /** Retrieve the abscissa of a channel. The abscissa is between -1 and 1, -1 is the left of the soundfield, 0 is the center of the soundfield and 1 is the right of the soundfield. The maximum index must be the number of channels - 1.
         
         @param     index    The index of the channel.
         @return    The abscissa of the channel if the channel exists, otherwise the function generates an error.
         
         @see getChannelAzimuth
         @see getChannelOrdinate
         @see getChannelName
         */
        inline double getChannelAbscissa(unsigned int index) const
        {
            if(m_mode == Regular)
                return m_decoder_regular->getChannelAbscissa(index);
            else
                return m_decoder_binaural->getChannelAbscissa(index);
        }
        
        //! Retrieve the ordinate of a channel.
        /** Retrieve the ordinate of a channel. The ordinate is between -1 and 1, -1 is the back of the soundfield, 0 is the center of the soundfield and 1 is the front of the soundfield. The maximum index must be the number of channels - 1.
         
         @param     index	The index of the channel.
         @return    The ordinate of the channel if the channel exists, otherwise the function generates an error.
         
         @see getChannelAzimuth
         @see getChannelAbscissa
         @see getChannelName
         */
        inline double getChannelOrdinate(unsigned int index) const
        {
            if(m_mode == Regular)
                return m_decoder_regular->getChannelOrdinate(index);
            else
                return m_decoder_binaural->getChannelOrdinate(index);
        }
        
        //! Retrieve the height of a channel.
        /** Retrieve the height of a channel. The height is between -1 and 1, -1 is the bottom of the soundfield, 0 is the center of the soundfield and 1 is the top of the soundfield. The maximum index must be the number of channels - 1.
         
         @param      index	The index of the channel.
         @return     The height of the channel.
         @see getChannelAbscissa
         @see getChannelOrdinate
         */
        inline double getChannelHeight(unsigned int index) const
        {
            if(m_mode == Regular)
                return m_decoder_regular->getChannelHeight(index);
            else
                return m_decoder_binaural->getChannelHeight(index);
        }
        
        
        //! Retrieve a name for a channel.
        /** Retrieve a name for a channel in a std::string format, look at each decoder for further informations.
         
         @param     index	The index of a channel.
         @return    The method returns a name for the channel.
         
         */
        inline std::string getChannelName(unsigned int index)
        {
            if(m_mode == Regular)
                return m_decoder_regular->getChannelName(index);
            else
                return m_decoder_binaural->getChannelName(index);
        };
        
        //! This method performs the decoding depending of the mode with single precision.
        /**	You should use this method for not-in-place processing and performs the binaural decoding on block of samples. The inputs matrix contains the spherical harmonics samples : inputs[number of harmonics][vector size] and the outputs matrix contains the headphones samples : outputs[2][vector size].
         
         @param     inputs	The input samples.
         @param     outputs  The output matrix that contains samples destinated to channels.
         */
        inline void process(const float* inputs, float* outputs)
        {
            if(m_mode == Regular)
                m_decoder_regular->process(inputs, outputs);
            else
                m_decoder_binaural->process(inputs, outputs);
        }
        
        //! This method performs the decoding depending of the mode with double precision.
        /**	You should use this method for in-place or not-in-place processing and performs the regular decoding sample by sample. The inputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics and the outputs array contains the channels samples and the minimym size must be the number of channels.
         
         @param     input	The input sample.
         @param     outputs The output array that contains samples destinated to channels.
         */
        inline void process(const double* inputs, double* outputs)
        {
            if(m_mode == Regular)
                m_decoder_regular->process(inputs, outputs);
            else
                m_decoder_binaural->process(inputs, outputs);
        }
    };

} // end of namespace Hoa3D

#endif
