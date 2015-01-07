/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_3D_VECTOR
#define DEF_HOA_3D_VECTOR

#include "Planewaves.h"

namespace Hoa2D
{
    //! The ambisonic vector.
    /** The vector class compute the energy and the velocity vector of a soudfield for a set of channels. It is an useful tool to characterize the quality of the sound field resitution. For futher information : Michael A. Gerzon, General metatheorie of auditory localisation. Audio Engineering Society Preprint, 3306, 1992. This class retreive the cartesian coordinates of the vectors, the abscissa and the ordinate.
     */
    class Vector : public Planewaves
    {
        
    private:
        double* m_channels_abscissa_double;
        double* m_channels_ordinate_double;
        double* m_channels_double;
        
        float* m_channels_abscissa_float;
        float* m_channels_ordinate_float;
        float* m_channels_float;
        double m_offset;
    public:
        
        //! The vector constructor.
        /**	The optimization constructor allocates and initialize the member values to computes vectors. The number of channels must be at least 1.
         
            @param     numberOfChannels	The number of channels.
         */
        Vector(unsigned int numberOfChannels);
        
        //! The optimization destructor.
        /**	The optimization destructor free the memory.
         */
        ~Vector();
        
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
        
        //! Set the azimtuh of a channel.
        /** Set the azimtuh of a channel. The azimtuh is in radian between 0 and 2 Pi, O is the front of the soundfield and Pi is the back of the sound field. The maximum index must be the number of channel - 1.
         
            @param     index		The index of the channel.
            @param     azimuth		The azimuth.
         */
		void setChannelAzimuth(unsigned int index, double azimuth);
		
		//! Set the azimtuh of all the channels.
        /** Set the azimtuh of all the channels. The azimuths are in radian between 0 and 2 Pi, O is the front of the soundfield and Pi is the back of the sound field. The azimtuhs array must have a minimum size of the number of channels.
         
		 @param     azimuths		The azimuths array.
         
		 @see    setChannelAzimuth
         */
		void setChannelsAzimuth(double* azimuth);
        
        //! This method compute the energy and velocity vectors with single precision.
        /**	You should use this method for in-place or not-in-place processing and compute the vectors sample by sample. The inputs array contains the channels samples and the minimum size must be the number of channels. The outputs array contains the vectors cartesian coordinates and the minimum size must be 4. The coordinates arrangement in the outputs array is velocity abscissa, velocity ordinate, energy abscissa, energy ordinate.
         
            @param     inputs   The inputs array.
            @param     outputs  The outputs array.
         */
        void process(const float* inputs, float* outputs);
        
        //! This method compute the energy and velocity vectors with double precision.
        /**	You should use this method for in-place or not-in-place processing and compute the vectors sample by sample. The inputs array contains the channels samples and the minimum size must be the number of channels. The outputs array contains the vectors cartesian coordinates and the minimum size must be 4. The coordinates arrangement in the outputs array is velocity abscissa, velocity ordinate, energy abscissa, energy ordinate.
         
            @param     inputs   The inputs array.
            @param     outputs  The outputs array.
         */
        void process(const double* inputs, double* outputs);
        
        //! This method compute the velocity vector with single precision.
        /**	You should use this method for in-place or not-in-place processing and compute the vectors sample by sample. The inputs array and contains the channels samples and the minimum size must be the number of channels. The outputs array contains the vectors cartesian coordinates and the minimum size must be 2. The coordinates arrangement in the outputs array is velocity abscissa, velocity ordinate.
         
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        void processVelocity(const float* inputs, float* outputs);
        
        //! This method compute the velocity vector with double precision.
        /**	You should use this method for in-place or not-in-place processing and compute the vectors sample by sample. The inputs array and contains the channels samples and the minimum size must be the number of channels. The outputs array contains the vectors cartesian coordinates and the minimum size must be 2. The coordinates arrangement in the outputs array is velocity abscissa, velocity ordinate.
         
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        void processVelocity(const double* inputs, double* outputs);
        
        //! This method compute the energy vector with single precision.
        /**	You should use this method for in-place or not-in-place processing and compute the vectors sample by sample. The inputs array and contains the channels samples and the minimum size must be the number of harmonics. The outputs array contains the vectors cartesian coordinates and the minimum size must be 2. The coordinates arrangement in the outputs array is energy abscissa, energy ordinate.
         
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        void processEnergy(const float* inputs, float* outputs);
        
        //! This method compute the energy vector with double precision.
        /**	You should use this method for in-place or not-in-place processing and compute the vectors sample by sample. The inputs array and contains the channels samples and the minimum size must be the number of harmonics. The outputs array contains the vectors cartesian coordinates and the minimum size must be 2. The coordinates arrangement in the outputs array is energy abscissa, energy ordinate.
         
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        void processEnergy(const double* inputs, double* outputs);
    };
}

#endif



