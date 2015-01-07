/*
 // Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
 // For information on usage and redistribution, and for a DISCLAIMER OF ALL
 // WARRANTIES, see the file, "LICENSE.txt," in this distribution.
 */

#include "Meter_3D.h"

namespace Hoa3D
{
    Meter::Meter(unsigned int numberOfChannels, unsigned int numberOfRows, unsigned int numberOfColumns) : Planewaves(numberOfChannels)
    {
        m_ramp                  = 0;
        m_vector_size           = 256;
        m_channels_peaks		= new double[m_number_of_channels];
        m_number_of_rows        = numberOfRows;
        m_number_of_columns     = numberOfColumns;
        if(m_number_of_rows % 2 != 1)
            m_number_of_rows++;
        if(m_number_of_columns % 2 != 0)
            m_number_of_columns++;
        
        for(unsigned int i = 0; i < numberOfChannels; i++)
            m_channels_peaks[i] = 0.;
        
		setChannelPosition(0, m_channels_azimuth[0], m_channels_elevation[0]);
    }
    
    void Meter::setVectorSize(unsigned int vectorSize)
    {
        m_vector_size   = vectorSize;
        m_ramp          = 0;
    }
    
	void Meter::setChannelsPosition(double* azimuths, double* elevations)
	{
		Planewaves::setChannelsPosition(azimuths, elevations);
		setChannelPosition(0, m_channels_azimuth[0], m_channels_elevation[0]);
	}
    
    void Meter::setChannelPosition(unsigned int index, double azimuth, double elevation)
	{
        Planewaves::setChannelPosition(index, azimuth, elevation);
		find_channels_boundaries();
	}
    
    void Meter::setChannelsRotation(double axis_x, double axis_y, double axis_z)
    {
        Planewaves::setChannelsRotation(axis_x, axis_y, axis_z);
        setChannelPosition(0, m_channels_azimuth[0], m_channels_elevation[0]);
    }
    
	void Meter::find_channels_boundaries()
	{
		int indices[8];
		double dist1, dist2, azi, ele;
		unsigned int numberOfRows = m_number_of_rows;
		unsigned int numberOfColumns = m_number_of_columns;
        unsigned int numberOfChannels = m_number_of_channels;
		int* sphere = new int[numberOfRows * numberOfColumns];

		double* azimuths = m_channels_rotated_azimuth;
		double* elevations = m_channels_rotated_elevation;
        
        // Fill a matrix that discretize a sphere with the indices of the closest loudspeakers
        for(unsigned int i = 0; i < numberOfRows; i++)
        {
            for(unsigned int j = 0; j < numberOfColumns; j++)
            {
                azi = (double)j / (double)numberOfColumns * HOA_2PI;
                ele = (double)(i * HOA_PI) / (double)(numberOfRows - 1) - HOA_PI2;
                dist1 = distance_spherical(azimuths[0], elevations[0], azi, ele);
                sphere[i * numberOfColumns + j] = 0;
                for(unsigned int k = 0; k < numberOfChannels; k++)
                {
                    dist2 = distance_spherical(azimuths[k], elevations[k], azi, ele);
                    if(dist2 < dist1)
                    {
                        dist1 = dist2;
                        sphere[i * numberOfColumns + j] = k;
                    }
                }
            }
        }
        
        for(unsigned int l = 0; l < numberOfChannels; l++)
        {
            MeterPoint ChannelPt = MeterPoint(azimuths[l], elevations[l]);
            m_points_top[l].clear();
            m_points_bottom[l].clear();
            
            bool test = 1;
            for(unsigned int j = 0; j < numberOfColumns && test; j++)
            {
                if(sphere[j] == l)
                {
                    for(unsigned int k = 0; k < numberOfColumns; k++)
                    {
                        if(sphere[numberOfColumns + k] != l)
                        {
                            MeterPoint Pt = MeterPoint(0, -HOA_PI2);
                            Pt.setRelativePoint(ChannelPt);
                            m_points_bottom[l].push_back(MeterPoint(Pt));
                            test = 0;
                            break;
                        }
                    }
                }
            }
            for(unsigned int i = 1; i < numberOfRows - 1; i++)
            {
                for(unsigned int j = 0; j < numberOfColumns; j++)
                {
                    if(sphere[i * numberOfColumns + j] == l)
                    {
                        int j0 = j-1, j1 = j - 1, j2 = j -1, j3 = j, j4 = j + 1, j5 = j + 1, j6 = j + 1, j7 = j;
                        int i0 = i + 1, i1 = i, i2 = i - 1, i3 = i - 1, i4 = i - 1, i5 = i, i6 = i + 1, i7 = i + 1;
                        
                        if(j == 0)
                            j0 = j1 = j2 = numberOfColumns-1;
                        else if(j == numberOfColumns - 1)
                            j4 = j5 = j6 = 0;
                        
                        indices[0] = sphere[i0 * numberOfColumns + j0];
                        indices[1] = sphere[i1 * numberOfColumns + j1];
                        indices[2] = sphere[i2 * numberOfColumns + j2];
                        indices[3] = sphere[i3 * numberOfColumns + j3];
                        indices[4] = sphere[i4 * numberOfColumns + j4];
                        indices[5] = sphere[i5 * numberOfColumns + j5];
                        indices[6] = sphere[i6 * numberOfColumns + j6];
                        indices[7] = sphere[i7 * numberOfColumns + j7];
                        
                        if(i == 90)
                        {
                            azi = (double)j / (double)numberOfColumns * HOA_2PI;
                            MeterPoint Pt = MeterPoint(azi, 0);
                            Pt.setRelativePoint(ChannelPt);
                            m_points_top[l].push_back(Pt);
                            m_points_bottom[l].push_back(Pt);
                        }
                        else
                        {
                            for(int k = 0; k < 8; k++)
                            {
                                if(indices[k] != l)
                                {
                                    azi = (double)j / (double)numberOfColumns * HOA_2PI;
                                    ele = (double)i * HOA_PI / (double)(numberOfRows - 1) - HOA_PI2;
                                    MeterPoint Pt = MeterPoint(azi, ele);
                                    Pt.setRelativePoint(ChannelPt);
                                    if(ele > 0)
                                        m_points_top[l].push_back(Pt);
                                    else
                                        m_points_bottom[l].push_back(Pt);
                                    break;
                                }
                            }
                        }
                        
                    }
                }
            }
            test = 1;
            for(unsigned int j = 0; j < numberOfColumns && test; j++)
            {
                if(sphere[(numberOfRows - 1) * numberOfColumns + j] == l)
                {
                    for(unsigned int k = 0; k < numberOfColumns; k++)
                    {
                        if(sphere[(numberOfRows - 2) * numberOfColumns + k] != l)
                        {
                            MeterPoint Pt = MeterPoint(0, HOA_PI2);
                            Pt.setRelativePoint(ChannelPt);
                            m_points_top[l].push_back(MeterPoint(Pt));
                            test = 0;
                            break;
                        }
                    }

                }
            }
            test = 1;
            int size = m_points_top[l].size();
            for(unsigned int i = 0; i < size; i++)
            {
                if(m_points_top[l][i].elevation() > 0.)
                {
                    test = 0;
                    break;
                }
            }
            if(test)
                m_points_top[l].clear();
            
            test = 1;
            size = m_points_bottom[l].size();
            for(unsigned int i = 0; i < size; i++)
            {
                if(m_points_bottom[l][i].elevation() < 0.)
                {
                    test = 0;
                    break;
                }
            }
            if(test)
                m_points_bottom[l].clear();
            
            std::sort(m_points_top[l].begin(), m_points_top[l].end(), MeterPoint::compareRelativeAzimuth);
            std::sort(m_points_bottom[l].begin(), m_points_bottom[l].end(), MeterPoint::compareRelativeAzimuth);
        }
        
		delete [] sphere;
	}
    
    void Meter::process(const float* inputs)
    {
        if(m_ramp++ == m_vector_size)
        {
            m_ramp = 0;
            for(unsigned int i = 0; i < m_number_of_channels; i++)
            {
                m_channels_peaks[i] = fabs(inputs[i]);
            }
        }
        else
        {
            for(unsigned int i = 0; i < m_number_of_channels; i++)
            {
                if(fabs(inputs[i]) > m_channels_peaks[i])
                    m_channels_peaks[i] = fabs(inputs[i]);
            }
        }
    }
    
    void Meter::process(const double* inputs)
    {
        if(m_ramp++ == m_vector_size)
        {
            m_ramp = 0;
            for(unsigned int i = 0; i < m_number_of_channels; i++)
            {
                m_channels_peaks[i] = fabs(inputs[i]);
            }
        }
        else
        {
            for(unsigned int i = 0; i < m_number_of_channels; i++)
            {
                if(fabs(inputs[i]) > m_channels_peaks[i])
                    m_channels_peaks[i] = fabs(inputs[i]);
            }
        }
    }
    
    Meter::~Meter()
    {
		for(unsigned int i = 0; i < m_number_of_channels; i++)
        {
            //m_points[i].clear();
            m_points_top[i].clear();
            m_points_bottom[i].clear();
        }
    }
}
