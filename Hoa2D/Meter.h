/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_2D_Meter
#define DEF_HOA_2D_Meter

#include "Planewaves.h"
#include "Vector.h"

namespace Hoa2D
{
    //! The planewaves peak level meter.
    /** The meter should be used to compute and display channels peak levels.
     */
    class Meter : public Planewaves
    {
    private:
        unsigned int    m_ramp;
        unsigned int    m_vector_size;
        double*         m_channels_peaks;
		double*			m_channels_azimuth_mapped;
		double*			m_channels_azimuth_width;
        double          m_offset;
        
		void computeAngles();
    public:
        
        //! The meter constructor.
        /**	The meter constructor allocates and initialize the member values.
         
         @param     number_of_channels	The number of channels.
         */
        Meter(unsigned int  number_of_channels);
        
        //! The meter destructor.
        /**	The meter destructor free the memory.
         */
        ~Meter();

        //! Set the azimuth of a given channel.
        /** Set the azimuth of a given channel. The azimuth is in radian between 0 and 2 Pi, O is the front of the soundfield and Pi is the back of the sound field.
         *
         * @param     index		The index of the channel.
         * @param     azimuth		The azimuth.
		 * @see getChannelAzimuth
         */
		void setChannelAzimuth(unsigned int index, double azimuth);
		
		//! Set the azimtuh of all the channels.
        /** Set the azimtuh of all the channels. The azimuths are in radian between 0 and 2 Pi, O is the front of the soundfield and Pi is the back of the sound field. The azimtuhs array must have a minimum size of the number of channels.
         
		 @param     azimuths		The azimuths array.
         
		 @see    setChannelAzimuth
         */
		void setChannelsAzimuth(double* azimuth);
        
        //! Set the offset of the channels.
		/**	Set the azimuth offset of the channels in radian.
         
         @param     offset		An azimuth value.
         */
		void setChannelsOffset(double offset);
        
        //! Get the offset of the channels.
        /**	Retreive the azimuth offset of the channels in radian.
         
         @return    The offset of the channels.
         */
		double getChannelsOffset() const
        {
            return m_offset;
        }
        
		//! Get the mapped azimuth of a given channel.
        /**
         *
         * @param     index		The index of the channel.
         * @return				The mapped azimuth.
		 * @see getChannelAzimuth, getChannelAzimuthMapped
         */
        double getChannelAzimuthMapped(unsigned int index) const
		{
			assert(index < m_number_of_channels);
            return m_channels_azimuth_mapped[index];
		}
		
		//! Get the width of a given channel.
        /**
         *
         * @param     index		The index of the channel.
         * @return				The channel width.
		 * @see getChannelAzimuth, getChannelAzimuthMapped
         */
        double getChannelWidth(unsigned int index) const
		{
			assert(index < m_number_of_channels);
            return m_channels_azimuth_width[index];
		}
        
        //! Set the vector size.
        /**
         * @param     vectorSize	The vector size.
		 * @see getChannelAzimuth
         */
        void setVectorSize(unsigned int vectorSize);
		
		//! Get the vector size.
        /**
         * @param     vectorSize	The vector size.
		 * @see setVectorSize
         */
        unsigned int getVectorSize(unsigned int vectorSize) const {return m_vector_size;}
        
		//! Retrieve the peak value of a given channel.
        /**
         * @param		index		The index of the channel.
		 * @see getChannelEnergy
         */
        inline double getChannelPeak(unsigned int index) const
        {
            assert(index < m_number_of_channels);
            return m_channels_peaks[index];
        }
        
		//! Retrieve the energy of a given channel.
        /**
         * @param		index		The index of the channel.
		 * @see getChannelPeak
         */
        inline double getChannelEnergy(unsigned int index) const
        {
            assert(index < m_number_of_channels);
            return atodb(m_channels_peaks[index]);
        }
        
        //! This method performs the metering with single precision.
        /**	You should use this method for in-place or not-in-place processing and performs the metering sample by sample. The inputs array contains channels samples and the minimum size must be the number of channels.
         
            @param     inputs   The inputs array.
         */
        void process(const float* inputs);
        
        //! This method performs the metering with double precision.
        /**	You should use this method for in-place or not-in-place processing and performs the metering sample by sample. The inputs array contains channels samples and the minimum size must be the number of channels.
         
            @param     inputs   The inputs array.
         */
        void process(const double* inputs);
    };
}

#endif



