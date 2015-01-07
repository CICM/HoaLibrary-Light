/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef __DEF_HOA_3D_METER__
#define __DEF_HOA_3D_METER__

#include "Planewaves_3D.h"

namespace Hoa3D
{
    class MeterPoint
    {
		private :
        
        double xyzae[5];
        double xyzae_rel[5];
        
        public :
        // Plutot faire en spherique - de calculs
        MeterPoint(double y_azimuth = 0, double z_elevation = 0)
        {
            xyzae[0] = Hoa::abscissa(1., y_azimuth, z_elevation);
            xyzae[1] = Hoa::ordinate(1., y_azimuth, z_elevation);
            xyzae[2] = Hoa::elevation(1., y_azimuth, z_elevation);
            xyzae[3] = y_azimuth;
            xyzae[4] = z_elevation;
        }
        
        void setRelativePoint(MeterPoint const& pt)
        {
            xyzae_rel[0] = x() - pt.x();
            xyzae_rel[1] = y() - pt.y();
            xyzae_rel[2] = z();
            xyzae_rel[3] = Hoa::azimuth(xyzae_rel[0], xyzae_rel[1]);
            xyzae_rel[4] = Hoa::elevation(xyzae_rel[0], xyzae_rel[1], xyzae_rel[2]);
        }
        
        double x() const
        {
            return xyzae[0];
        }
        
        double y() const
        {
            return xyzae[1];
        }
        
        double z() const
        {
            return xyzae[2];
        }
        
        double radius() const
        {
            return 1;
        }
        
        double azimuth() const
        {
            return xyzae[3];
        }
        
        double elevation() const
        {
            return xyzae[4];
        }
        
        double x_rel() const
        {
            return xyzae_rel[0];
        }
        
        double y_rel() const
        {
            return xyzae_rel[1];
        }
        
        double z_rel() const
        {
            return xyzae_rel[2];
        }
        
        double radius_rel() const
        {
            return 1;
        }
        
        double azimuth_rel() const
        {
            return xyzae_rel[3];
        }
        
        static bool compareRelativeAzimuth(MeterPoint pt1, MeterPoint pt2)
        {
            return pt1.azimuth_rel() > pt2.azimuth_rel();
        }
        
        ~MeterPoint(){};
    };
    
    //! The planewaves peak level meter.
    /** The meter should be used to widen the sound propagation.
     */
    class Meter : public Planewaves
    {
    private:
        unsigned int    m_ramp;
        unsigned int    m_vector_size;
        unsigned int    m_number_of_rows;
        unsigned int    m_number_of_columns;
        double*         m_channels_peaks;
        
        std::vector<MeterPoint> m_points_top[256];
        std::vector<MeterPoint> m_points_bottom[256];

		void find_channels_boundaries();
    public:
        
        //! The meter constructor.
        /**	The meter constructor allocates and initialize the member values to computes spherical harmonics weighted coefficients depending of a decomposition order. The order must be at least 1.
         
            @param     order	The order.
         */
        Meter(unsigned int numberOfChannels, unsigned int numberOfRows, unsigned int numberOfColumns);
        
        //! The meter destructor.
        /**	The meter destructor free the memory.
         */
        ~Meter();
        
        void setVectorSize(unsigned int vectorSize);
        
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
        
		inline unsigned int getChannelNumberOfPoints(unsigned int index, bool top = 1) const
        {
            assert(index < m_number_of_channels);
			if(top)
				return m_points_top[index].size();
			else
				return m_points_bottom[index].size();
        }

		inline double getChannelPointAzimuth(unsigned int index, unsigned int pointindex, bool top = 1) const
        {
            assert(index < m_number_of_channels);
            if(top)
				return m_points_top[index][pointindex].azimuth();
			else
				return m_points_bottom[index][pointindex].azimuth();
            
        }

		inline double getChannelPointElevation(unsigned int index, unsigned int pointindex, bool top = 1) const
        {
            assert(index < m_number_of_channels);
            if(top)
				return m_points_top[index][pointindex].elevation();
			else
				return m_points_bottom[index][pointindex].elevation();
        }

        inline double getChannelPeak(unsigned int index) const
        {
            assert(index < m_number_of_channels);
            return m_channels_peaks[index];
        }
        
        inline double getChannelEnergy(unsigned int index) const
        {
            assert(index < m_number_of_channels);
            if(m_channels_peaks[index] == 0.)
                return -91;
            else
                return clip_min(20. * log10(m_channels_peaks[index]), -90.);
        }
       
        //! This method performs the widening with single precision.
        /**	You should use this method for in-place or not-in-place processing and performs the widening sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         
            @param     inputs   The inputs array.
         */
        void process(const float* inputs);
        
        //! This method performs the widening with double precision.
        /**	You should use this method for in-place or not-in-place processing and performs the widening sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         
            @param     inputs   The inputs array.
         */
        void process(const double* inputs);
    };
}

#endif



