/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include "SourcesGroup.h"

namespace HoaCommon
{
	SourcesGroup::SourcesGroup(SourcesManager* sourcesManager, bool state)
	{
		m_source_manager = sourcesManager;
		setExistence(state);
		m_color = new double[4];
		setColor(0.2, 0.2, 0.2, 1.);
		setDescription("");
		computeCentroid();
		m_maximum_radius = -1;
		m_mute = 0;
	}
	
	void SourcesGroup::setExistence(bool state)
	{
		m_exist = state;
	}
	
	void SourcesGroup::setDescription(std::string description)
	{
		m_description = description;
	}
	
	void SourcesGroup::setColor(double red, double green, double blue, double alpha)
	{
		m_color[0]	=  clip_minmax(red, 0., 1.);
		m_color[1]	=  clip_minmax(green, 0., 1.);
		m_color[2]	=  clip_minmax(blue, 0., 1.);
		m_color[3]	=  clip_minmax(alpha, 0., 1.);
	}
	
	void SourcesGroup::setMaximumRadius(double limitValue)
	{
		m_maximum_radius = clip_min(limitValue, 0.0000001);
	}
	
	void SourcesGroup::computeCentroid()
	{
		m_centroid_x = 0.;
		m_centroid_y = 0.;
        m_centroid_z = 0.;
		if(m_sources.size())
		{
			for(int i = 0; i < m_sources.size(); i++)
			{
				if(m_source_manager->sourceGetExistence(m_sources[i]))
				{
					m_centroid_x += m_source_manager->sourceGetAbscissa(m_sources[i]);
					m_centroid_y += m_source_manager->sourceGetOrdinate(m_sources[i]);
                    m_centroid_z += m_source_manager->sourceGetHeight(m_sources[i]);
				}
			}
			m_centroid_x /= m_sources.size();
			m_centroid_y /= m_sources.size();
            m_centroid_z /= m_sources.size();
		}
	}
	
	void SourcesGroup::addSource(long sourceIndex)
	{
		for(int i = 0; i < m_sources.size(); i++)
		{
			if(m_sources[i] == sourceIndex)
				return;
		}
		m_sources.push_back(sourceIndex);
		
		computeCentroid();
	}
	
	void SourcesGroup::sourceHasMoved()
	{
		computeCentroid();
	}
	
	void SourcesGroup::removeSource(long sourceIndex)
	{
		if(m_sources.size() > 0)
		{
			int size = m_sources.size();
			
			for(int i = 0; i < size; i++)
			{
				if(m_sources[i] == sourceIndex)
				{
					for(int j = i; j < size - 1; j++)
					{
						m_sources[j] = m_sources[j+1];
					}
					m_sources.pop_back();
				}
			}
		}
		computeCentroid();
	}
	
	void SourcesGroup::shiftPolar(double radius, double azimuth)
	{
		shiftRadius(radius);
		shiftAzimuth(azimuth);
	}
    
    void SourcesGroup::shiftPolar(double radius, double azimuth, double elevation)
	{
		shiftRadius(radius);
		shiftAzimuth(azimuth);
        shiftElevation(elevation);
	}
	
	void SourcesGroup::shiftRadius(double radius)
	{
		if(m_maximum_radius >= 0)
		{
			if(radius < 0.)
			{
				double refRadius = m_maximum_radius;
				for(int i = 0; i < m_sources.size(); i++)
				{
					if(m_source_manager->sourceGetRadius(m_sources[i]) < refRadius)
					{
						refRadius = m_source_manager->sourceGetRadius(m_sources[i]);
					}
				}
				if(radius + refRadius < 0.)
				{
					radius = -refRadius;
				}
			}
			else if(radius >= 0.)
			{
				double refRadius = -m_maximum_radius;
				for(int i = 0; i < m_sources.size(); i++)
				{
					if(m_source_manager->sourceGetRadius(m_sources[i]) > refRadius)
					{
						refRadius = m_source_manager->sourceGetRadius(m_sources[i]);
					}
				}
				if(radius + refRadius > m_maximum_radius)
				{
					radius = m_maximum_radius - refRadius;
				}
			}
		}
		for(int i = 0; i < m_sources.size(); i++)
		{
			m_source_manager->sourceSetRadius(m_sources[i], radius + m_source_manager->sourceGetRadius(m_sources[i]));
		}
	}
	
	void SourcesGroup::shiftAzimuth(double azimuth)
	{
		for(int i = 0; i < m_sources.size(); i++)
		{
			m_source_manager->sourceSetAzimuth(m_sources[i], azimuth + m_source_manager->sourceGetAzimuth(m_sources[i]));
		}
	}
    
    void SourcesGroup::shiftElevation(double elevation)
	{
		for(int i = 0; i < m_sources.size(); i++)
		{
			m_source_manager->sourceSetElevation(m_sources[i], elevation + m_source_manager->sourceGetElevation(m_sources[i]));
		}
	}
	
    void SourcesGroup::shiftCartesian(double abscissa, double ordinate)
	{
		if(m_maximum_radius >= 0)
		{
			if(abscissa < 0. &&  ordinate < 0.)
			{
				double refAbcsissa = -m_maximum_radius * 2.;
				double refOrdinate = -m_maximum_radius * 2.;
				for(int i = 0; i < m_sources.size(); i++)
				{
					double circleAbscissa = -sqrt(m_maximum_radius * m_maximum_radius - m_source_manager->sourceGetOrdinate(m_sources[i]) * m_source_manager->sourceGetOrdinate(m_sources[i]));
					double circleOrdinate = -sqrt(m_maximum_radius * m_maximum_radius - m_source_manager->sourceGetAbscissa(m_sources[i]) * m_source_manager->sourceGetAbscissa(m_sources[i]));
					if(circleAbscissa - m_source_manager->sourceGetAbscissa(m_sources[i]) > refAbcsissa)
					{
						refAbcsissa = circleAbscissa - m_source_manager->sourceGetAbscissa(m_sources[i]);
					}
					if(circleOrdinate - m_source_manager->sourceGetOrdinate(m_sources[i]) > refOrdinate)
					{
						refOrdinate = circleOrdinate - m_source_manager->sourceGetOrdinate(m_sources[i]);
					}
				}
				if(abscissa < refAbcsissa)
				{
					abscissa = refAbcsissa;
				}
				if(ordinate < refOrdinate)
				{
					ordinate = refOrdinate;
				}
			}
			else if(abscissa >= 0.)
			{
				double refAbcsissa = m_maximum_radius * 2.;
				for(int i = 0; i < m_sources.size(); i++)
				{
					double circleAbscissa = sqrt(m_maximum_radius * m_maximum_radius - m_source_manager->sourceGetOrdinate(m_sources[i]) * m_source_manager->sourceGetOrdinate(m_sources[i]));
					if(circleAbscissa - m_source_manager->sourceGetAbscissa(m_sources[i]) < refAbcsissa)
					{
						refAbcsissa = circleAbscissa - m_source_manager->sourceGetAbscissa(m_sources[i]);
					}
				}
				if(abscissa > refAbcsissa)
				{
					abscissa = refAbcsissa;
				}
			}
		}
		
		for(int i = 0; i < m_sources.size(); i++)
		{
			m_source_manager->sourceSetAbscissa(m_sources[i], abscissa + m_source_manager->sourceGetAbscissa(m_sources[i]));
			m_source_manager->sourceSetOrdinate(m_sources[i], ordinate + m_source_manager->sourceGetOrdinate(m_sources[i]));
		}
	}
    
    void SourcesGroup::shiftCartesian(double abscissa, double ordinate, double height)
	{
		if(m_maximum_radius >= 0)
		{
			if(abscissa < 0. &&  ordinate < 0.)
			{
				double refAbcsissa = -m_maximum_radius * 2.;
				double refOrdinate = -m_maximum_radius * 2.;
				for(int i = 0; i < m_sources.size(); i++)
				{
					double circleAbscissa = -sqrt(m_maximum_radius * m_maximum_radius - m_source_manager->sourceGetOrdinate(m_sources[i]) * m_source_manager->sourceGetOrdinate(m_sources[i]));
					double circleOrdinate = -sqrt(m_maximum_radius * m_maximum_radius - m_source_manager->sourceGetAbscissa(m_sources[i]) * m_source_manager->sourceGetAbscissa(m_sources[i]));
					if(circleAbscissa - m_source_manager->sourceGetAbscissa(m_sources[i]) > refAbcsissa)
					{
						refAbcsissa = circleAbscissa - m_source_manager->sourceGetAbscissa(m_sources[i]);
					}
					if(circleOrdinate - m_source_manager->sourceGetOrdinate(m_sources[i]) > refOrdinate)
					{
						refOrdinate = circleOrdinate - m_source_manager->sourceGetOrdinate(m_sources[i]);
					}
				}
				if(abscissa < refAbcsissa)
				{
					abscissa = refAbcsissa;
				}
				if(ordinate < refOrdinate)
				{
					ordinate = refOrdinate;
				}
			}
			else if(abscissa >= 0.)
			{
				double refAbcsissa = m_maximum_radius * 2.;
				for(int i = 0; i < m_sources.size(); i++)
				{
					double circleAbscissa = sqrt(m_maximum_radius * m_maximum_radius - m_source_manager->sourceGetOrdinate(m_sources[i]) * m_source_manager->sourceGetOrdinate(m_sources[i]));
					if(circleAbscissa - m_source_manager->sourceGetAbscissa(m_sources[i]) < refAbcsissa)
					{
						refAbcsissa = circleAbscissa - m_source_manager->sourceGetAbscissa(m_sources[i]);
					}
				}
				if(abscissa > refAbcsissa)
				{
					abscissa = refAbcsissa;
				}
			}
		}
		
		for(int i = 0; i < m_sources.size(); i++)
		{
			m_source_manager->sourceSetAbscissa(m_sources[i], abscissa + m_source_manager->sourceGetAbscissa(m_sources[i]));
			m_source_manager->sourceSetOrdinate(m_sources[i], ordinate + m_source_manager->sourceGetOrdinate(m_sources[i]));
		}
	}
    
	void SourcesGroup::shiftAbscissa(double abscissa)
	{
		if(m_maximum_radius >= 0)
		{
			if(abscissa < 0.)
			{
				double refAbcsissa = -m_maximum_radius * 2.;
				for(int i = 0; i < m_sources.size(); i++)
				{
					double circleAbscissa = -sqrt(m_maximum_radius * m_maximum_radius - m_source_manager->sourceGetOrdinate(m_sources[i]) * m_source_manager->sourceGetOrdinate(m_sources[i]));
					if(circleAbscissa - m_source_manager->sourceGetAbscissa(m_sources[i]) > refAbcsissa)
					{
						refAbcsissa = circleAbscissa - m_source_manager->sourceGetAbscissa(m_sources[i]);
					}
				}
				if(abscissa < refAbcsissa)
				{
					abscissa = refAbcsissa;
				}
			}
			else if(abscissa >= 0.)
			{
				double refAbcsissa = m_maximum_radius * 2.;
				for(int i = 0; i < m_sources.size(); i++)
				{
					double circleAbscissa = sqrt(m_maximum_radius * m_maximum_radius - m_source_manager->sourceGetOrdinate(m_sources[i]) * m_source_manager->sourceGetOrdinate(m_sources[i]));
					if(circleAbscissa - m_source_manager->sourceGetAbscissa(m_sources[i]) < refAbcsissa)
					{
						refAbcsissa = circleAbscissa - m_source_manager->sourceGetAbscissa(m_sources[i]);
					}
				}
				if(abscissa > refAbcsissa)
				{
					abscissa = refAbcsissa;
				}
			}
		}
		
		for(int i = 0; i < m_sources.size(); i++)
		{
			m_source_manager->sourceSetAbscissa(m_sources[i], abscissa + m_source_manager->sourceGetAbscissa(m_sources[i]));
		}
	}
	
	void SourcesGroup::shiftOrdinate(double ordinate)
	{
		if(m_maximum_radius >= 0)
		{
			if(ordinate < 0.)
			{
				double refOrdinate = -m_maximum_radius * 2.;
				for(int i = 0; i < m_sources.size(); i++)
				{
					double circleOrdinate = -sqrt(m_maximum_radius * m_maximum_radius - m_source_manager->sourceGetAbscissa(m_sources[i]) * m_source_manager->sourceGetAbscissa(m_sources[i]));
					if(circleOrdinate - m_source_manager->sourceGetOrdinate(m_sources[i]) > refOrdinate)
					{
						refOrdinate = circleOrdinate - m_source_manager->sourceGetOrdinate(m_sources[i]);
					}
				}
				if(ordinate < refOrdinate)
				{
					ordinate = refOrdinate;
				}
			}
			else if(ordinate >= 0.)
			{
				double refOrdinate = m_maximum_radius * 2.;
				for(int i = 0; i < m_sources.size(); i++)
				{
					double circleOrdinate = sqrt(m_maximum_radius * m_maximum_radius - m_source_manager->sourceGetAbscissa(m_sources[i]) * m_source_manager->sourceGetAbscissa(m_sources[i]));
					if(circleOrdinate - m_source_manager->sourceGetOrdinate(m_sources[i]) < refOrdinate)
					{
						refOrdinate = circleOrdinate - m_source_manager->sourceGetOrdinate(m_sources[i]);
					}
				}
				if(ordinate > refOrdinate)
				{
					ordinate = refOrdinate;
				}
			}
		}
		
		for(int i = 0; i < m_sources.size(); i++)
		{
			m_source_manager->sourceSetOrdinate(m_sources[i], ordinate + m_source_manager->sourceGetOrdinate(m_sources[i]));
		}
	}
    
    void SourcesGroup::shiftHeight(double height)
	{
		if(m_maximum_radius >= 0)
		{
			if(height < 0.)
			{
				double refHeight = -m_maximum_radius * 2.;
				for(int i = 0; i < m_sources.size(); i++)
				{
					double circleHeight = -sqrt(m_maximum_radius * m_maximum_radius - m_source_manager->sourceGetAbscissa(m_sources[i]) * m_source_manager->sourceGetAbscissa(m_sources[i]));
					if(circleHeight - m_source_manager->sourceGetHeight(m_sources[i]) > circleHeight)
					{
						refHeight = circleHeight - m_source_manager->sourceGetHeight(m_sources[i]);
					}
				}
				if(height < refHeight)
				{
					height = refHeight;
				}
			}
			else if(height >= 0.)
			{
				double refHeight = m_maximum_radius * 2.;
				for(int i = 0; i < m_sources.size(); i++)
				{
					double circleHeight = sqrt(m_maximum_radius * m_maximum_radius - m_source_manager->sourceGetAbscissa(m_sources[i]) * m_source_manager->sourceGetAbscissa(m_sources[i]));
					if(circleHeight - m_source_manager->sourceGetHeight(m_sources[i]) < refHeight)
					{
						refHeight = circleHeight - m_source_manager->sourceGetHeight(m_sources[i]);
					}
				}
				if(height > refHeight)
				{
					height = refHeight;
				}
			}
		}
		
		for(int i = 0; i < m_sources.size(); i++)
		{
			m_source_manager->sourceSetHeight(m_sources[i], height + m_source_manager->sourceGetHeight(m_sources[i]));
		}
	}
	
	void SourcesGroup::setCoordinatesPolar(double radius, double azimuth)
	{
		setCoordinatesCartesian(abscissa(radius, azimuth), ordinate(radius, azimuth));
	}
    
    void SourcesGroup::setCoordinatesPolar(double radius, double azimuth, double elevation)
	{
		setCoordinatesCartesian(abscissa(radius, azimuth, elevation), ordinate(radius, azimuth, elevation));
	}
	
	void SourcesGroup::setRadius(double radius)
	{
        setCoordinatesCartesian(abscissa(radius, getAzimuth(), getElevation()), ordinate(radius, getAzimuth(), getElevation()), height(radius, getAzimuth(), getElevation()));
	}
	
	void SourcesGroup::setAzimuth(double azimuth)
	{
        setCoordinatesCartesian(abscissa(getRadius(), azimuth, getElevation()), ordinate(getRadius(), azimuth, getElevation()), height(getRadius(), azimuth, getElevation()));
	}
    
    void SourcesGroup::setElevation(double elevation)
	{
		setCoordinatesCartesian(abscissa(getRadius(), getAzimuth(), elevation), ordinate(getRadius(), getAzimuth(), elevation), height(getRadius(), getAzimuth(), elevation));
	}
	
	void SourcesGroup::setCoordinatesCartesian(double abscissa, double ordinate)
	{
		abscissa = abscissa - getAbscissa();
		ordinate = ordinate - getOrdinate();
		shiftAbscissa(abscissa);
		shiftOrdinate(ordinate);
		computeCentroid();
	}
    
    void SourcesGroup::setCoordinatesCartesian(double abscissa, double ordinate, double height)
	{
		abscissa = abscissa - getAbscissa();
		ordinate = ordinate - getOrdinate();
        height = height - getHeight();
		shiftAbscissa(abscissa);
		shiftOrdinate(ordinate);
        shiftHeight(height);
		computeCentroid();
	}
	
	void SourcesGroup::setAbscissa(double abscissa)
	{
		double aAbscissaOffset = abscissa - getAbscissa();
		shiftAbscissa(aAbscissaOffset);
		computeCentroid();
	}
	
	void SourcesGroup::setOrdinate(double ordinate)
	{
		double aOrdinateOffset = ordinate - getOrdinate();
		shiftOrdinate(aOrdinateOffset);
		computeCentroid();
	}
    
    void SourcesGroup::setHeight(double height)
	{
		double aHeightOffset = height - getHeight();
		shiftHeight(aHeightOffset);
		computeCentroid();
	}
	
	void SourcesGroup::setRelativeCoordinatesPolar(double radius, double azimuth)
	{
		setRelativeRadius(radius);
		setRelativeAzimuth(azimuth);
	}
    
    void SourcesGroup::setRelativeCoordinatesPolar(double radius, double azimuth, double elevation)
	{
		setRelativeRadius(radius);
		setRelativeAzimuth(azimuth);
        setRelativeElevation(elevation);
	}
	
	void SourcesGroup::setRelativeRadius(double radius)
	{
		double aRadiusOffset = clip_min(radius, 0.) - getRadius();
		shiftRadius(aRadiusOffset);
		computeCentroid();
	}
	
	void SourcesGroup::setRelativeAzimuth(double azimuth)
	{
		azimuth +=  HOA_PI2;
		while (azimuth > HOA_2PI)
			azimuth -= HOA_2PI;
		while (azimuth < 0.)
			azimuth += HOA_2PI;
		
		double aAngleOffset = azimuth  - getAzimuth();
		shiftAzimuth(aAngleOffset);
		computeCentroid();
	}
    
    void SourcesGroup::setRelativeElevation(double elevation)
	{
		elevation +=  HOA_PI2;
		while (elevation > HOA_2PI)
			elevation -= HOA_2PI;
		while (elevation < 0.)
			elevation += HOA_2PI;
		
		double aAngleOffset = elevation  - getElevation();
		shiftElevation(aAngleOffset);
		computeCentroid();
	}
	
	void SourcesGroup::setMute(long aValue)
	{
		m_mute = clip_minmax(aValue, (long)0, (long)1);
	}
	
	SourcesGroup::~SourcesGroup()
	{
		delete m_color;
		m_sources.clear();
	}
}