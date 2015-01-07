/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef __DEF_HOA_3D_VECTOR__
#define __DEF_HOA_3D_VECTOR__

#include "Planewaves_3D.h"

namespace Hoa3D
{
    //! The planewaves vector.
    /** The vector class compute the energy and the velocity vector of a soudfield with signal of a spjerical set of channels. It is an useful tool to characterize the quality of the sound field resitution. For futher information : Michael A. Gerzon, General metatheorie of auditory localisation. Audio Engineering Society Preprint, 3306, 1992. This class retreive the cartesian coordinates of the vectors, the abscissa, the ordinate and the height.
     */
    class Vector : public Planewaves
    {
        
    private:
        double* m_channels_abscissa_double;
        double* m_channels_ordinate_double;
        double* m_channels_height_double;
        double* m_channels_double;
        
        float* m_channels_abscissa_float;
        float* m_channels_ordinate_float;
        float* m_channels_height_float;
        float* m_channels_float;
    public:
        
        //! The vector constructor.
        /**	The optimization constructor allocates and initialize the member values to computes vectors. The numberOfChannels must be at least 1.
         
            @param     numberOfChannels	The number of channels.
         */
        Vector(unsigned int numberOfChannels);
        
        //! The optimization destructor.
        /**	The optimization destructor free the memory.
         */
        ~Vector();
        
        //! Set the position of a channel.
        /** Set the position of a channel with polar coordinates. The azimtuh is in radian between 0 and 2 Pi, O is the front of the soundfield and Pi is the back of the sound field. The elevation is in radian between -1/2 Pi and 1/2 Pi, -1/2 Pi the the bottom of the sound field, 0 is the center of the sound field and 1/2 Pi is the top of the sound field. The maximum index must be the number of channels - 1.
         
         @param     index		The index of the channel.
         @param     azimuth		The azimuth.
         @param     elevation	The elevation.
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
        
        //! This method compute the energy and velocity vectors with single precision.
        /**	You should use this method for in-place or not-in-place processing and compute the vectors sample by sample. The inputs array and contains the spherical harmonics samples and the minimum size must be the number of harmonics. The outputs array contains the vectors cartesian coordinates and the minimum size must be 6. The coordinates arrengement in the outputs array is velocity abscissa, velocity ordinate, velocity height, energy abscissa, energy ordinate, energy height.
         
            @param     inputs   The inputs array.
            @param     outputs  The outputs array.
         */
        void process(const float* inputs, float* outputs);
        
        //! This method compute the energy and velocity vectors with double precision.
        /**	You should use this method for in-place or not-in-place processing and compute the vectors sample by sample. The inputs array and contains the spherical harmonics samples and the minimum size must be the number of harmonics. The outputs array contains the vectors cartesian coordinates and the minimum size must be 6. The coordinates arrengement in the outputs array is velocity abscissa, velocity ordinate, velocity height, energy abscissa, energy ordinate, energy height.
         
            @param     inputs   The inputs array.
            @param     outputs  The outputs array.
         */
        void process(const double* inputs, double* outputs);
        
        //! This method compute the velocity vector with single precision.
        /**	You should use this method for in-place or not-in-place processing and compute the vectors sample by sample. The inputs array and contains the spherical harmonics samples and the minimum size must be the number of harmonics. The outputs array contains the vectors cartesian coordinates and the minimum size must be 3. The coordinates arrengement in the outputs array is velocity abscissa, velocity ordinate, velocity height.
         
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        void processVelocity(const float* inputs, float* outputs);
        
        //! This method compute the velocity vector with double precision.
        /**	You should use this method for in-place or not-in-place processing and compute the vectors sample by sample. The inputs array and contains the spherical harmonics samples and the minimum size must be the number of harmonics. The outputs array contains the vectors cartesian coordinates and the minimum size must be 3. The coordinates arrengement in the outputs array is velocity abscissa, velocity ordinate, velocity height.
         
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        void processVelocity(const double* inputs, double* outputs);
        
        //! This method compute the energy vector with single precision.
        /**	You should use this method for in-place or not-in-place processing and compute the vectors sample by sample. The inputs array and contains the spherical harmonics samples and the minimum size must be the number of harmonics. The outputs array contains the vectors cartesian coordinates and the minimum size must be 3. The coordinates arrengement in the outputs array is energy abscissa, energy ordinate, energy height.
         
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        void processEnergy(const float* inputs, float* outputs);
        
        //! This method compute the energy vector with double precision.
        /**	You should use this method for in-place or not-in-place processing and compute the vectors sample by sample. The inputs array and contains the spherical harmonics samples and the minimum size must be the number of harmonics. The outputs array contains the vectors cartesian coordinates and the minimum size must be 3. The coordinates arrengement in the outputs array is energy abscissa, energy ordinate, energy height.
         
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        void processEnergy(const double* inputs, double* outputs);
    };
}

#endif



