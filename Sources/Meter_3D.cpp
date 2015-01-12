/*
 // Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
 // For information on usage and redistribution, and for a DISCLAIMER OF ALL
 // WARRANTIES, see the file, "LICENSE.txt," in this distribution.
 */

#include "Meter_3D.h"

namespace Hoa3D
{
    Planewaves::Planewaves(unsigned int numberOfChannels)
    {
        assert(numberOfChannels > 0);
        m_number_of_channels    = numberOfChannels;
        m_channels_azimuth      = new double[m_number_of_channels];
        m_channels_elevation    = new double[m_number_of_channels];
        m_channels_rotated_azimuth      = new double[m_number_of_channels];
        m_channels_rotated_elevation    = new double[m_number_of_channels];
        sphere_discretize(numberOfChannels, m_channels_azimuth, m_channels_elevation);
        setChannelsRotation(0, 0, 0);
    }
    
    void Planewaves::setChannelsPosition(double* azimuths, double* elevations)
    {
        for(unsigned int i = 0; i < m_number_of_channels; i++)
        {
            m_channels_azimuth[i] = wrap_twopi(azimuths[i]);
            m_channels_elevation[i] = wrap(elevations[i], -HOA_PI, HOA_PI);
            if(m_channels_elevation[i] > HOA_PI2)
            {
                m_channels_azimuth[i] = wrap_twopi(m_channels_azimuth[i] + HOA_PI);
                m_channels_elevation[i] = HOA_PI2 - (m_channels_elevation[i] - HOA_PI2);
            }
            else if (m_channels_elevation[i] < -HOA_PI2)
            {
                m_channels_azimuth[i] = wrap_twopi(m_channels_azimuth[i] + HOA_PI);
                m_channels_elevation[i] = -HOA_PI2 + (-m_channels_elevation[i] + HOA_PI2);
            }
            rotateChannel(i);
        }
    }
    
    void Planewaves::setChannelPosition(unsigned int index, double azimuth, double elevation)
    {
        assert(index < m_number_of_channels);
        m_channels_azimuth[index] = wrap_twopi(azimuth);
        m_channels_elevation[index] = wrap(elevation, -HOA_PI, HOA_PI);
        
        if(m_channels_elevation[index] > HOA_PI2)
        {
            m_channels_azimuth[index] = wrap_twopi(m_channels_azimuth[index] + HOA_PI);
            m_channels_elevation[index] -= HOA_PI2;
        }
        else if (m_channels_elevation[index] < -HOA_PI2)
        {
            m_channels_azimuth[index] = wrap_twopi(m_channels_azimuth[index] + HOA_PI);
            m_channels_elevation[index] += HOA_PI2;
        }
        rotateChannel(index);
    }
    
    void Planewaves::setChannelsRotation(double axis_x, double axis_y, double axis_z)
    {
        m_channels_rotation_x = wrap_twopi(axis_x);
        m_channels_rotation_y = wrap_twopi(axis_y);
        m_channels_rotation_z = wrap_twopi(axis_z);
        
        for(unsigned int i = 0; i < m_number_of_channels; i++)
            rotateChannel(i);
    }
    
    void Planewaves::rotateChannel(unsigned int index)
    {
        double x = getChannelAbscissa(index);
        double y = getChannelOrdinate(index);
        double z = getChannelHeight(index);
        
        // Rotation around x
        double cosAngle = cos(m_channels_rotation_x);
        double sinAngle = sin(m_channels_rotation_x);
        double ry = y * cosAngle - z * sinAngle;
        double rz = y * sinAngle + z * cosAngle;
        y = ry;
        z = rz;
        
        // Rotation around z
        cosAngle = cos(m_channels_rotation_z);
        sinAngle = sin(m_channels_rotation_z);
        double rx = x * cosAngle - y * sinAngle;
        ry = x * sinAngle + y * cosAngle;
        x = rx;
        y = ry;
        
        // Rotation around y
        cosAngle = cos(m_channels_rotation_y);
        sinAngle = sin(m_channels_rotation_y);
        rx = x * cosAngle - z * sinAngle;
        rz = x * sinAngle + z * cosAngle;
        x = rx;
        z = rz;
        
        m_channels_rotated_azimuth[index] = azimuth(x, y, z);
        m_channels_rotated_elevation[index] = elevation(x, y, z);
    }
    
    Planewaves::~Planewaves()
    {
        delete [] m_channels_azimuth;
        delete [] m_channels_elevation;
        delete [] m_channels_rotated_azimuth;
        delete [] m_channels_rotated_elevation;
    }
    
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
    
    Vector::Vector(unsigned int numberOfChannels) : Planewaves(numberOfChannels)
    {
        m_channels_float = new float[m_number_of_channels];
        m_channels_double = new double[m_number_of_channels];
        m_channels_abscissa_float = new float[m_number_of_channels];
        m_channels_abscissa_double = new double[m_number_of_channels];
        m_channels_ordinate_float = new float[m_number_of_channels];
        m_channels_ordinate_double = new double[m_number_of_channels];
        m_channels_height_float = new float[m_number_of_channels];
        m_channels_height_double = new double[m_number_of_channels];
        setChannelPosition(0, m_channels_azimuth[0], m_channels_elevation[0]);
    }
    
    void Vector::setChannelPosition(unsigned int index, double azimuth, double elevation)
    {
        Planewaves::setChannelPosition(index, azimuth, elevation);
        for(unsigned int i = 0; i < m_number_of_channels; i++)
        {
            m_channels_abscissa_float[i] = m_channels_abscissa_double[i] = getChannelRotatedAbscissa(i);
            m_channels_ordinate_float[i] = m_channels_ordinate_double[i] = getChannelRotatedOrdinate(i);
            m_channels_height_float[i] = m_channels_height_double[i] = getChannelRotatedHeight(i);
        }
    }
    
    void Vector::setChannelsPosition(double* azimuths, double* elevations)
    {
        Planewaves::setChannelsPosition(azimuths, elevations);
        setChannelPosition(0, m_channels_azimuth[0], m_channels_elevation[0]);
    }
    
    void Vector::setChannelsRotation(double axis_x, double axis_y, double axis_z)
    {
        Planewaves::setChannelsRotation(axis_x, axis_y, axis_z);
        setChannelPosition(0, m_channels_azimuth[0], m_channels_elevation[0]);
    }
    
    void Vector::processVelocity(const float* inputs, float* outputs)
    {
        float veclocitySum = 0.f, velocityAbscissa = 0.f, velocityOrdinate = 0.f, velocityElevation = 0.f;
        
        veclocitySum = cblas_sasum(m_number_of_channels, inputs, 1);
        velocityAbscissa = cblas_sdot(m_number_of_channels, inputs, 1, m_channels_abscissa_float, 1);
        velocityOrdinate = cblas_sdot(m_number_of_channels, inputs, 1, m_channels_ordinate_float, 1);
        velocityElevation= cblas_sdot(m_number_of_channels, inputs, 1, m_channels_height_float, 1);
        if(veclocitySum)
        {
            outputs[0] = velocityAbscissa / veclocitySum;
            outputs[1] = velocityOrdinate / veclocitySum;
            outputs[2] = velocityElevation / veclocitySum;
        }
        else
        {
            outputs[0] = 0.;
            outputs[1] = 0.;
            outputs[2] = 0.;
            
        }
    }
    
    void Vector::processVelocity(const double* inputs, double* outputs)
    {
        double veclocitySum = 0., velocityAbscissa = 0., velocityOrdinate = 0., velocityElevation = 0.;
        veclocitySum = cblas_dasum(m_number_of_channels, inputs, 1);
        velocityAbscissa = cblas_ddot(m_number_of_channels, inputs, 1, m_channels_abscissa_double, 1);
        velocityOrdinate = cblas_ddot(m_number_of_channels, inputs, 1, m_channels_ordinate_double, 1);
        velocityElevation= cblas_ddot(m_number_of_channels, inputs, 1, m_channels_height_double, 1);
        
        if(veclocitySum)
        {
            outputs[0] = velocityAbscissa / veclocitySum;
            outputs[1] = velocityOrdinate / veclocitySum;
            outputs[2] = velocityElevation / veclocitySum;
        }
        else
        {
            outputs[0] = 0.;
            outputs[1] = 0.;
            outputs[2] = 0.;
        }
    }
    
    void Vector::processEnergy(const float* inputs, float* outputs)
    {
        float energySum = 0.f, energyAbscissa = 0.f, energyOrdinate = 0.f, energyElevation = 0.f;
        cblas_scopy(m_number_of_channels, inputs, 1, m_channels_float, 1);
        for(unsigned int i = 0; i < m_number_of_channels; i++)
            m_channels_float[i] *= m_channels_float[i];
        
        energySum = cblas_sasum(m_number_of_channels, m_channels_float, 1);
        energyAbscissa = cblas_sdot(m_number_of_channels, m_channels_float, 1, m_channels_abscissa_float, 1);
        energyOrdinate = cblas_sdot(m_number_of_channels, m_channels_float, 1, m_channels_ordinate_float, 1);
        energyElevation = cblas_sdot(m_number_of_channels, m_channels_float, 1, m_channels_height_float, 1);
        
        if(energySum)
        {
            outputs[0] = energyAbscissa / energySum;
            outputs[1] = energyOrdinate / energySum;
            outputs[2] = energyElevation / energySum;
        }
        else
        {
            outputs[0] = 0.;
            outputs[1] = 0.;
            outputs[2] = 0.;
        }
    }
    
    void Vector::processEnergy(const double* inputs, double* outputs)
    {
        double energySum = 0., energyAbscissa = 0., energyOrdinate = 0., energyElevation = 0.;
        
        cblas_dcopy(m_number_of_channels, inputs, 1, m_channels_double, 1);
        for(unsigned int i = 0; i < m_number_of_channels; i++)
            m_channels_double[i] *= m_channels_double[i];
        
        energySum = cblas_dasum(m_number_of_channels, m_channels_double, 1);
        energyAbscissa = cblas_ddot(m_number_of_channels, m_channels_double, 1, m_channels_abscissa_double, 1);
        energyOrdinate = cblas_ddot(m_number_of_channels, m_channels_double, 1, m_channels_ordinate_double, 1);
        energyElevation = cblas_ddot(m_number_of_channels, m_channels_double, 1, m_channels_height_double, 1);
        
        if(energySum)
        {
            outputs[0] = energyAbscissa / energySum;
            outputs[1] = energyOrdinate / energySum;
            outputs[2] = energyElevation / energySum;
        }
        else
        {
            outputs[0] = 0.;
            outputs[1] = 0.;
            outputs[2] = 0.;
        }
    }
    
    void Vector::process(const float* inputs, float* outputs)
    {
        processVelocity(inputs, outputs);
        processEnergy(inputs, outputs+3);
    }
    
    void Vector::process(const double* inputs, double* outputs)
    {
        processVelocity(inputs, outputs);
        processEnergy(inputs, outputs+3);
    }
    
    Vector::~Vector()
    {
        delete [] m_channels_float;
        delete [] m_channels_double;
        delete [] m_channels_ordinate_float;
        delete [] m_channels_ordinate_double;
        delete [] m_channels_abscissa_float;
        delete [] m_channels_abscissa_double;
        delete [] m_channels_height_float;
        delete [] m_channels_height_double;
    }
}
