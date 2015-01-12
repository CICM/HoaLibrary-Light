/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/


#include "Source.hpp"

namespace hoa
{
    
    Source::Source(bool existence, double radius, double azimuth, double elevation)
    {
        m_exist = existence;
        setRadius(radius);
        setAzimuth(azimuth);
        setElevation(elevation);
        m_color = new double[4];
        setColor(0.2, 0.2, 0.2, 1.);
        setDescription("");
        m_maximum_radius = -1;
        m_mute = 0;
    }
    
    void Source::setMaximumRadius(double limitValue)
    {
        m_maximum_radius = limitValue;
    }
    
    void Source::setExistence(bool state)
    {
        m_exist = state;
    }
    
    void Source::setCoordinatesPolar(double radius, double azimuth)
    {
        setRadius(radius);
        setAzimuth(azimuth);
    }
    
    void Source::setCoordinatesPolar(double radius, double azimuth, double elevation)
    {
        setRadius(radius);
        setAzimuth(azimuth);
        setElevation(elevation);
    }
    
    void Source::setRadius(double radius)
    {
        if(m_maximum_radius >= 0)
        {
            if(radius < -m_maximum_radius || radius > m_maximum_radius)
                return;
        }
        m_radius = max(radius, 0.);
    }
    
    void Source::setAzimuth(double azimuth)
    {
        m_azimuth = Math<double>::wrap_twopi(azimuth);
    }
    
    void Source::setElevation(double elevation)
    {
        m_elevation = Math<double>::wrap(elevation, -HOA_PI, HOA_PI);
        if(m_elevation > HOA_PI2)
        {
            m_azimuth = Math<double>::wrap_twopi(m_azimuth + HOA_PI);
            m_elevation = HOA_PI2 - (elevation - HOA_PI2);
        }
        else if(m_elevation < -HOA_PI2)
        {
            m_azimuth = Math<double>::wrap_twopi(m_azimuth + HOA_PI);
            m_elevation = -HOA_PI2 + (-elevation + HOA_PI2);
        }
    }
    
    void Source::setCoordinatesCartesian(double abscissa, double ordinate)
    {
        setRadius(Math<double>::radius(abscissa, ordinate));
        setAzimuth(Math<double>::azimuth(abscissa, ordinate));
    }
    
    void Source::setCoordinatesCartesian(double abscissa, double ordinate, double height)
    {
        setRadius(Math<double>::radius(abscissa, ordinate, height));
        setAzimuth(Math<double>::azimuth(abscissa, ordinate, height));
        setHeight(Math<double>::elevation(abscissa, ordinate, height));
    }
    
    void Source::setAbscissa(double abscissa)
    {
        double ordinate = getOrdinate();
        double height = getHeight();
        setRadius(Math<double>::radius(abscissa, ordinate, height));
        setAzimuth(Math<double>::azimuth(abscissa, ordinate, height));
        setElevation(Math<double>::elevation(abscissa, ordinate, height));
    }
    
    void Source::setOrdinate(double ordinate)
    {
        double abscissa = getAbscissa();
        double height = getHeight();
        setRadius(Math<double>::radius(abscissa, ordinate, height));
        setAzimuth(Math<double>::azimuth(abscissa, ordinate, height));
        setElevation(Math<double>::elevation(abscissa, ordinate, height));
    }
    
    void Source::setHeight(double height)
    {
        double abscissa = getAbscissa();
        double ordinate = getOrdinate();
        setRadius(Math<double>::radius(abscissa, ordinate, height));
        setAzimuth(Math<double>::azimuth(abscissa, ordinate, height));
        setElevation(Math<double>::elevation(abscissa, ordinate, height));
    }
    
    void Source::setColor(double red, double green, double blue, double alpha)
    {
        m_color[0]	=  Math<double>::clip(red, 0., 1.);
        m_color[1]	=  Math<double>::clip(green, 0., 1.);
        m_color[2]	=  Math<double>::clip(blue, 0., 1.);
        m_color[3]	=  Math<double>::clip(alpha, 0., 1.);
    }
    
    void Source::setDescription(std::string description)
    {
        m_description = description;
    }
    
    void Source::setGroup(long groupIndex)
    {
        for(int i = 0; i < m_groups.size(); i++)
        {
            if(m_groups[i] == groupIndex)
                return;
        }
        m_groups.push_back(groupIndex);
    }
    
    void Source::setMute(bool state)
    {
        m_mute = state;
    }
    
    void Source::removeGroup(long groupIndex)
    {
        for(int i = 0; i < m_groups.size(); i++)
        {
            if(m_groups[i] == groupIndex)
            {
                for(int j = i; j < m_groups.size() - 1; j++)
                {
                    m_groups[j] = m_groups[j+1];
                }
                m_groups.pop_back();
            }
        }
    }
    
    long Source::getGroupIndex(long index)
    {
        if(index < m_groups.size() && index >= 0)
            return m_groups[index];
        else
            return -1;
    }
    
    bool Source::isOwnedByGroup(long groupIndex)
    {
        for (int i = 0; i < m_groups.size(); i++)
            if (m_groups[i] == groupIndex) return true;
        return false;
    }
    
    Source::~Source()
    {
        delete m_color;
        m_groups.clear();
    }
    
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
        m_color[0]	=  Math<double>::clip(red, 0., 1.);
        m_color[1]	=  Math<double>::clip(green, 0., 1.);
        m_color[2]	=  Math<double>::clip(blue, 0., 1.);
        m_color[3]	=  Math<double>::clip(alpha, 0., 1.);
    }
    
    void SourcesGroup::setMaximumRadius(double limitValue)
    {
        m_maximum_radius = max(limitValue, 0.0000001);
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
        setCoordinatesCartesian(Math<double>::abscissa(radius, azimuth), Math<double>::ordinate(radius, azimuth));
    }
    
    void SourcesGroup::setCoordinatesPolar(double radius, double azimuth, double elevation)
    {
        setCoordinatesCartesian(Math<double>::abscissa(radius, azimuth, elevation), Math<double>::ordinate(radius, azimuth, elevation));
    }
    
    void SourcesGroup::setRadius(double radius)
    {
        setCoordinatesCartesian(Math<double>::abscissa(radius, getAzimuth(), getElevation()), Math<double>::ordinate(radius, getAzimuth(), getElevation()), Math<double>::height(radius, getAzimuth(), getElevation()));
    }
    
    void SourcesGroup::setAzimuth(double azimuth)
    {
        setCoordinatesCartesian(Math<double>::abscissa(getRadius(), azimuth, getElevation()), Math<double>::ordinate(getRadius(), azimuth, getElevation()), Math<double>::height(getRadius(), azimuth, getElevation()));
    }
    
    void SourcesGroup::setElevation(double elevation)
    {
        setCoordinatesCartesian(Math<double>::abscissa(getRadius(), getAzimuth(), elevation), Math<double>::ordinate(getRadius(), getAzimuth(), elevation), Math<double>::height(getRadius(), getAzimuth(), elevation));
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
        double aRadiusOffset = max(radius, 0.) - getRadius();
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
        m_mute = Math<double>::clip(aValue, (long)0, (long)1);
    }
    
    SourcesGroup::~SourcesGroup()
    {
        delete m_color;
        m_sources.clear();
    }
    
    
    SourcesManager::SourcesManager(double maximumLimitValue, bool existence)
    {
        setExistence(existence);
        setMaximumRadius(maximumLimitValue);
        m_zoom = 1.;
    }
    
    void SourcesManager::setExistence(bool state)
    {
        m_exist = state;
        if(m_exist == 0)
        {
            m_sources.clear();
            m_groups.clear();
        }
    }
    
    void SourcesManager::clearAll()
    {
        
        for(int i = 0; i <= getMaximumIndexOfSource(); i++)
            sourceRemove(i);
        for(int i = 0; i <= getMaximumIndexOfGroup(); i++)
            groupRemove(i);
    }
    
    bool SourcesManager::getExistence()
    {
        return m_exist;
    }
    
    void SourcesManager::setMaximumRadius(double limitValue)
    {
        m_maximum_radius = limitValue;
        for(int i = 0; i < m_groups.size(); i++)
            m_groups[i]->setMaximumRadius(m_maximum_radius);
    }
    
    void SourcesManager::setZoom(double aZoom)
    {
        m_zoom = Math<double>::clip(aZoom, 1. / m_maximum_radius, 1.);
    }
    
    double SourcesManager::getZoom()
    {
        return m_zoom;
    }
    
    double SourcesManager::getLimitMaximum()
    {
        return m_maximum_radius;
    }
    
    long SourcesManager::getMaximumIndexOfSource()
    {
        long index = 0;
        for (int i = 0; i < m_sources.size(); i++)
        {
            if(sourceGetExistence(i))
                index = i;
        }
        return  index;
    }
    
    long SourcesManager::getNumberOfSources()
    {
        long numberOfSources = 0;
        for(int i = 0; i < m_sources.size(); i++)
        {
            if(m_sources[i]->getExistence())
                numberOfSources++;
        }
        return  numberOfSources;
    }
    
    long SourcesManager::getMaximumIndexOfGroup()
    {
        return  m_groups.size();
    }
    
    long SourcesManager::getNumberOfGroups()
    {
        long numberOfGroups = 0;
        for(int i = 0; i < m_groups.size(); i++)
        {
            if (m_groups[i]->getExistence())
                numberOfGroups++;
        }
        return  numberOfGroups;
    }
    
    /*******************************************************************************/
    /*********************************  SOURCES ************************************/
    /*******************************************************************************/
    
    void SourcesManager::sourceRemove(long index)
    {
        if(index < m_sources.size() && index >= 0)
        {
            int numberOfGroups = m_sources[index]->getNumberOfGroups();
            int indexOfGroup = -1;
            for(int i = 0; i < numberOfGroups; i++)
            {
                indexOfGroup = m_sources[index]->getGroupIndex(i);
                
                if(indexOfGroup >= 0 && indexOfGroup <= getMaximumIndexOfGroup())
                {
                    if(groupGetExistence(indexOfGroup))
                        m_groups[indexOfGroup]->removeSource(index);
                }
            }
            m_sources[index]->setExistence(0);
            m_sources[index]->setDescription("");
            m_sources[index]->setColor(0.2, 0.2, 0.2, 1.);
            m_sources[index]->setCoordinatesCartesian(0., 1.);
            m_sources[index]->setMute(0);
        }
    }
    
    void SourcesManager::sourceNewPolar(double radius, double azimuth)
    {
        for (int i = 0; i <= getMaximumIndexOfSource()+1; i++)
        {
            if(!sourceGetExistence(i))
            {
                sourceSetPolar(i, radius, azimuth);
                return;
            }
        }
    }
    
    void SourcesManager::sourceNewPolar(double radius, double azimuth, double elevation)
    {
        for (int i = 0; i <= getMaximumIndexOfSource()+1; i++)
        {
            if(!sourceGetExistence(i))
            {
                sourceSetPolar(i, radius, azimuth, elevation);
                return;
            }
        }
    }
    
    void SourcesManager::sourceNewCartesian(double abscissa, double ordinate)
    {
        for(int i = 0; i <= getMaximumIndexOfSource()+1; i++)
        {
            if(!sourceGetExistence(i))
            {
                sourceSetCartesian(i, abscissa, ordinate);
                return;
            }
        }
    }
    
    void SourcesManager::sourceNewCartesian(double abscissa, double ordinate, double height)
    {
        for(int i = 0; i <= getMaximumIndexOfSource()+1; i++)
        {
            if(!sourceGetExistence(i))
            {
                sourceSetCartesian(i, abscissa, ordinate, height);
                return;
            }
        }
    }
    
    void SourcesManager::sourceSetPolar(long index, double radius, double azimuth)
    {
        sourceSetRadius(index, radius);
        sourceSetAzimuth(index, azimuth);
    }
    
    void SourcesManager::sourceSetPolar(long index, double radius, double azimuth, double elevation)
    {
        sourceSetRadius(index, radius);
        sourceSetAzimuth(index, azimuth);
        sourceSetElevation(index, elevation);
    }
    
    void SourcesManager::sourceSetRadius(long index, double radius)
    {
        if (index < 0)
            return;
        
        if(index >= m_sources.size())
        {
            for(int i = m_sources.size(); i < index; i++)
            {
                m_sources.push_back(new Source(0));
                m_sources[i]->setMaximumRadius(m_maximum_radius);
                
            }
            m_sources.push_back(new Source(1));
            m_sources[index]->setMaximumRadius(m_maximum_radius);
            m_sources[index]->setRadius(radius);
        }
        else if(index >= 0)
        {
            m_sources[index]->setRadius(radius);
            if(!m_sources[index]->getExistence())
                m_sources[index]->setExistence(1);
        }
        for(int i = 0; i < m_sources[index]->getNumberOfGroups(); i++)
        {
            int groupIndex = m_sources[index]->getGroupIndex(i);
            if(groupIndex >= 0 && groupIndex < m_groups.size())
                m_groups[groupIndex]->sourceHasMoved();
        }
    }
    
    void SourcesManager::sourceSetAzimuth(long index, double azimuth)
    {
        if (index < 0)
            return;
        
        if(index >= m_sources.size())
        {
            for(int i = m_sources.size(); i < index; i++)
            {
                m_sources.push_back(new Source(0));
                m_sources[i]->setMaximumRadius(m_maximum_radius);
            }
            m_sources.push_back(new Source(1));
            m_sources[index]->setMaximumRadius(m_maximum_radius);
            m_sources[index]->setAzimuth(azimuth);
        }
        else if(index >= 0)
        {
            m_sources[index]->setAzimuth(azimuth);
            if(!m_sources[index]->getExistence())
                m_sources[index]->setExistence(1);
        }
        for(int i = 0; i < m_sources[index]->getNumberOfGroups(); i++)
        {
            int groupIndex = m_sources[index]->getGroupIndex(i);
            if(groupIndex >= 0 && groupIndex < m_groups.size())
                m_groups[groupIndex]->sourceHasMoved();
        }
    }
    
    void SourcesManager::sourceSetElevation(long index, double elevation)
    {
        if (index < 0)
            return;
        
        if(index >= m_sources.size())
        {
            for(int i = m_sources.size(); i < index; i++)
            {
                m_sources.push_back(new Source(0));
                m_sources[i]->setMaximumRadius(m_maximum_radius);
            }
            m_sources.push_back(new Source(1));
            m_sources[index]->setMaximumRadius(m_maximum_radius);
            m_sources[index]->setElevation(elevation);
        }
        else if(index >= 0)
        {
            m_sources[index]->setElevation(elevation);
            if(!m_sources[index]->getExistence())
                m_sources[index]->setExistence(1);
        }
        for(int i = 0; i < m_sources[index]->getNumberOfGroups(); i++)
        {
            int groupIndex = m_sources[index]->getGroupIndex(i);
            if(groupIndex >= 0 && groupIndex < m_groups.size())
                m_groups[groupIndex]->sourceHasMoved();
        }
    }
    
    void SourcesManager::sourceSetCartesian(long index, double abscissa, double ordinate)
    {
        sourceSetAbscissa(index, abscissa);
        sourceSetOrdinate(index, ordinate);
    }
    
    void SourcesManager::sourceSetCartesian(long index, double abscissa, double ordinate, double height)
    {
        sourceSetAbscissa(index, abscissa);
        sourceSetOrdinate(index, ordinate);
        sourceSetHeight(index, height);
    }
    
    void SourcesManager::sourceSetAbscissa(long index, double abscissa)
    {
        if (index < 0)
            return;
        
        if(index >= m_sources.size())
        {
            for(int i = m_sources.size(); i < index; i++)
            {
                m_sources.push_back(new Source(0));
                m_sources[i]->setMaximumRadius(m_maximum_radius);
            }
            m_sources.push_back(new Source(1));
            m_sources[index]->setMaximumRadius(m_maximum_radius);
            m_sources[index]->setAbscissa(abscissa);
        }
        else if(index >= 0)
        {
            m_sources[index]->setAbscissa(abscissa);
            if(!m_sources[index]->getExistence())
                m_sources[index]->setExistence(1);
        }
        for(int i = 0; i < m_sources[index]->getNumberOfGroups(); i++)
        {
            int groupIndex = m_sources[index]->getGroupIndex(i);
            if(groupIndex >= 0 && groupIndex < m_groups.size())
                m_groups[groupIndex]->sourceHasMoved();
        }
    }
    
    void SourcesManager::sourceSetOrdinate(long index, double ordinate)
    {
        if (index < 0)
            return;
        
        if(index >= m_sources.size())
        {
            for(int i = m_sources.size(); i < index; i++)
            {
                m_sources.push_back(new Source(0));
                m_sources[i]->setMaximumRadius(m_maximum_radius);
            }
            m_sources.push_back(new Source(1));
            m_sources[index]->setMaximumRadius(m_maximum_radius);
            m_sources[index]->setOrdinate(ordinate);
        }
        else if(index >= 0)
        {
            m_sources[index]->setOrdinate(ordinate);
            if(!m_sources[index]->getExistence())
                m_sources[index]->setExistence(1);
        }
        for(int i = 0; i < m_sources[index]->getNumberOfGroups(); i++)
        {
            int groupIndex = m_sources[index]->getGroupIndex(i);
            if(groupIndex >= 0 && groupIndex < m_groups.size())
                m_groups[groupIndex]->sourceHasMoved();
        }
    }
    
    void SourcesManager::sourceSetHeight(long index, double height)
    {
        if (index < 0)
            return;
        
        if(index >= m_sources.size())
        {
            for(int i = m_sources.size(); i < index; i++)
            {
                m_sources.push_back(new Source(0));
                m_sources[i]->setMaximumRadius(m_maximum_radius);
            }
            m_sources.push_back(new Source(1));
            m_sources[index]->setMaximumRadius(m_maximum_radius);
            m_sources[index]->setHeight(height);
        }
        else if(index >= 0)
        {
            m_sources[index]->setHeight(height);
            if(!m_sources[index]->getExistence())
                m_sources[index]->setExistence(1);
        }
        for(int i = 0; i < m_sources[index]->getNumberOfGroups(); i++)
        {
            int groupIndex = m_sources[index]->getGroupIndex(i);
            if(groupIndex >= 0 && groupIndex < m_groups.size())
                m_groups[groupIndex]->sourceHasMoved();
        }
    }
    
    void SourcesManager::sourceSetColor(long index, double red, double green, double blue, double alpha)
    {
        if(index < m_sources.size() && index >= 0)
        {
            m_sources[index]->setColor(red, green, blue, alpha);
        }
    }
    
    void SourcesManager::sourceSetDescription(long index, std::string description)
    {
        if(index < m_sources.size() && index >= 0)
        {
            m_sources[index]->setDescription(description);
        }
    }
    void SourcesManager::checkMute()
    {
        for(int i = 0; i < m_groups.size(); i++)
        {
            m_groups[i]->setMute(1);
            for(int j = 0; j < m_groups[i]->getNumberOfSources(); j++)
            {
                int sourceIndex = m_groups[i]->getSourceIndex(j);
                if(sourceIndex >= 0 && sourceIndex < m_sources.size())
                {
                    if(m_sources[sourceIndex]->getMute() != 1)
                    {
                        m_groups[i]->setMute(0);
                    }
                }
            }
        }
    }
    
    void SourcesManager::sourceSetMute(long index, bool state)
    {
        if(index < m_sources.size() && index >= 0)
        {
            m_sources[index]->setMute(state);
            for(int i = 0; i < m_sources[index]->getNumberOfGroups(); i++)
            {
                int groupIndex = m_sources[index]->getGroupIndex(i);
                if(groupIndex >= 0 && groupIndex < m_groups.size())
                {
                    m_groups[groupIndex]->setMute(0);
                }
            }
            checkMute();
        }
    }
    
    /******************************************************************************/
    
    double SourcesManager::sourceGetRadius(long index)
    {
        if(index < m_sources.size() && index >= 0)
            return m_sources[index]->getRadius();
        
        return 0;
    }
    
    double SourcesManager::sourceGetAzimuth(long index)
    {
        if(index < m_sources.size() && index >= 0)
            return m_sources[index]->getAzimuth();
        
        return 0;
    }
    
    double SourcesManager::sourceGetElevation(long index)
    {
        if(index < m_sources.size() && index >= 0)
            return m_sources[index]->getElevation();
        
        return 0;
    }
    
    double SourcesManager::sourceGetAbscissa(long index)
    {
        if(index < m_sources.size() && index >= 0)
            return m_sources[index]->getAbscissa();
        
        return NULL;
    }
    
    double SourcesManager::sourceGetOrdinate(long index)
    {
        if(index < m_sources.size() && index >= 0)
            return m_sources[index]->getOrdinate();
        
        return 0;
    }
    
    double SourcesManager::sourceGetHeight(long index)
    {
        if(index < m_sources.size() && index >= 0)
            return m_sources[index]->getHeight();
        
        return 0;
    }
    
    double* SourcesManager::sourceGetColor(long index)
    {
        if(index < m_sources.size() && index >= 0)
            return m_sources[index]->getColor();
        
        return NULL;
    }
    
    std::string SourcesManager::sourceGetDescription(long index)
    {
        if(index < m_sources.size() && index >= 0)
            return m_sources[index]->getDescription();
        
        return NULL;
    }
    
    long SourcesManager::sourceGetExistence(long index)
    {
        if(index < m_sources.size() && index >= 0)
        {
            return m_sources[index]->getExistence();
        }
        return 0;
    }
    
    long SourcesManager::sourceGetNumberOfGroups(long index)
    {
        if(index < m_sources.size() && index >= 0)
        {
            return m_sources[index]->getNumberOfGroups();
        }
        return 0;
    }
    
    long SourcesManager::sourceGetGroupIndex(long sourceIndex, long groupIndex)
    {
        if(sourceIndex < m_sources.size() && sourceIndex >= 0)
        {
            return m_sources[sourceIndex]->getGroupIndex(groupIndex);
        }
        return 0;
    }
    
    long SourcesManager::sourceGetMute(long index)
    {
        if(index < m_sources.size() && index >= 0)
        {
            return m_sources[index]->getMute();
        }
        return 0;
    }
    
    /*******************************************************************************/
    /**********************************  GROUP  ************************************/
    /*******************************************************************************/
    
    void SourcesManager::groupSetSource(long groupIndex, long sourceIndex)
    {
        if (sourceIndex < 0 || groupIndex < 0 )
            return;
        
        if(groupIndex >= m_groups.size())
        {
            for(int i = m_groups.size(); i < groupIndex; i++)
            {
                m_groups.push_back(new SourcesGroup(this, 0));
                m_groups[i]->setMaximumRadius(m_maximum_radius);
                
            }
            m_groups.push_back(new SourcesGroup(this, 0));
            m_groups[groupIndex]->setMaximumRadius(m_maximum_radius);
            if(m_sources.size() > sourceIndex && m_sources[sourceIndex]->getExistence())
            {
                m_groups[groupIndex]->setExistence(1);
                m_groups[groupIndex]->addSource(sourceIndex);
                m_sources[sourceIndex]->setGroup(groupIndex);
                checkMute();
            }
        }
        else if(groupIndex >= 0)
        {
            if(m_sources.size() > sourceIndex && m_sources[sourceIndex]->getExistence())
            {
                m_groups[groupIndex]->addSource(sourceIndex);
                m_sources[sourceIndex]->setGroup(groupIndex);
                m_groups[groupIndex]->setExistence(1);
                checkMute();
            }
        }
    }
    
    void SourcesManager::groupRemoveSource(long groupIndex, long sourceIndex)
    {
        if (sourceIndex < 0 || groupIndex < 0 )
            return;
        
        if(groupIndex < m_groups.size() && groupIndex >= 0)
        {
            if(m_sources.size() > sourceIndex)
            {
                m_groups[groupIndex]->removeSource(sourceIndex);
                m_sources[sourceIndex]->removeGroup(groupIndex);
            }
        }
    }
    
    void SourcesManager::groupSetPolar(long groupIndex, double radius, double azimuth)
    {
        if(groupIndex < m_groups.size() && groupIndex >= 0)
        {
            m_groups[groupIndex]->setCoordinatesPolar(radius, azimuth);
        }
    }
    
    void SourcesManager::groupSetPolar(long groupIndex, double radius, double azimuth, double elevation)
    {
        if(groupIndex < m_groups.size() && groupIndex >= 0)
        {
            m_groups[groupIndex]->setCoordinatesPolar(radius, azimuth, elevation);
        }
    }
    
    void SourcesManager::groupSetRadius(long groupIndex, double radius)
    {
        if(groupIndex < m_groups.size() && groupIndex >= 0)
        {
            m_groups[groupIndex]->setRadius(radius);
        }
    }
    
    void SourcesManager::groupSetAzimuth(long groupIndex, double azimuth)
    {
        if(groupIndex < m_groups.size() && groupIndex >= 0)
        {
            m_groups[groupIndex]->setAzimuth(azimuth);
        }
    }
    
    void SourcesManager::groupSetElevation(long groupIndex, double elevation)
    {
        if(groupIndex < m_groups.size() && groupIndex >= 0)
        {
            m_groups[groupIndex]->setElevation(elevation);
        }
    }
    
    void SourcesManager::groupSetCartesian(long groupIndex, double abscissa, double ordinate)
    {
        if(groupIndex < m_groups.size() && groupIndex >= 0)
        {
            m_groups[groupIndex]->setCoordinatesCartesian(abscissa, ordinate);
        }
    }
    
    void SourcesManager::groupSetCartesian(long groupIndex, double abscissa, double ordinate, double height)
    {
        if(groupIndex < m_groups.size() && groupIndex >= 0)
        {
            m_groups[groupIndex]->setCoordinatesCartesian(abscissa, ordinate, height);
        }
    }
    
    void SourcesManager::groupSetAbscissa(long groupIndex, double abscissa)
    {
        if(groupIndex < m_groups.size() && groupIndex >= 0)
        {
            m_groups[groupIndex]->setAbscissa(abscissa);
        }
    }
    
    void SourcesManager::groupSetOrdinate(long groupIndex, double ordinate)
    {
        if(groupIndex < m_groups.size() && groupIndex >= 0)
        {
            m_groups[groupIndex]->setOrdinate(ordinate);
        }
    }
    
    void SourcesManager::groupSetHeight(long groupIndex, double height)
    {
        if(groupIndex < m_groups.size() && groupIndex >= 0)
        {
            m_groups[groupIndex]->setHeight(height);
        }
    }
    
    void SourcesManager::groupSetRelativePolar(long groupIndex, double radius, double azimuth)
    {
        if(groupIndex < m_groups.size() && groupIndex >= 0)
        {
            m_groups[groupIndex]->setRelativeCoordinatesPolar(radius, azimuth);
        }
    }
    
    void SourcesManager::groupSetRelativePolar(long groupIndex, double radius, double azimuth, double elevation)
    {
        if(groupIndex < m_groups.size() && groupIndex >= 0)
        {
            m_groups[groupIndex]->setRelativeCoordinatesPolar(radius, azimuth, elevation);
        }
    }
    
    void SourcesManager::groupSetRelativeRadius(long groupIndex, double radius)
    {
        if(groupIndex < m_groups.size() && groupIndex >= 0)
        {
            m_groups[groupIndex]->setRelativeRadius(radius);
        }
    }
    
    void SourcesManager::groupSetRelativeAzimuth(long groupIndex, double azimuth)
    {
        if(groupIndex < m_groups.size() && groupIndex >= 0)
        {
            m_groups[groupIndex]->setRelativeAzimuth(azimuth);
        }
    }
    
    void SourcesManager::groupSetRelativeElevation(long groupIndex, double elevation)
    {
        if(groupIndex < m_groups.size() && groupIndex >= 0)
        {
            m_groups[groupIndex]->setRelativeElevation(elevation);
        }
    }
    
    void SourcesManager::groupSetColor(long groupIndex, double red, double green, double blue, double alpha)
    {
        if(groupIndex < m_groups.size() && groupIndex >= 0)
        {
            m_groups[groupIndex]->setColor(red, green, blue, alpha);
        }
    }
    
    void SourcesManager::groupSetDescription(long groupIndex, std::string description)
    {
        if(groupIndex < m_groups.size() && groupIndex >= 0)
        {
            m_groups[groupIndex]->setDescription(description);
        }
    }
    
    void SourcesManager::groupRemove(long groupIndex)
    {
        if(groupIndex < m_groups.size() && groupIndex >= 0)
        {
            for(int i = 0; i < m_sources.size(); i++)
            {
                m_sources[i]->removeGroup(groupIndex);
                
            }
            for(int i = 0; i < m_sources.size(); i++)
            {
                m_groups[groupIndex]->removeSource(i);
            }
            m_groups[groupIndex]->setColor(0.2, 0.2, 0.2, 1.);
            m_groups[groupIndex]->setDescription("");
            m_groups[groupIndex]->setExistence(0);
            m_groups[groupIndex]->setMute(0);
        }
    }
    
    void SourcesManager::groupRemoveWithSources(long groupIndex)
    {
        for(int i = 0; i <= getMaximumIndexOfSource(); i++)
        {
            if (m_sources[i]->isOwnedByGroup(groupIndex)) {
                sourceRemove(i);
            }
        }
        groupRemove(groupIndex);
    }
    
    void SourcesManager::groupSetMute(long groupIndex, long aValue)
    {
        if(groupIndex < m_groups.size() && groupIndex >= 0)
        {
            m_groups[groupIndex]->setMute(aValue);
            for(int i = 0; i < m_groups[groupIndex]->getNumberOfSources(); i++)
            {
                int sourceIndex = m_groups[groupIndex]->getSourceIndex(i);
                if(sourceIndex >= 0 && sourceIndex < m_sources.size())
                    m_sources[sourceIndex]->setMute(aValue);
            }
        }
        checkMute();
    }
    
    
    /******************************************************************************/
    
    double SourcesManager::groupGetRadius(long index)
    {
        if(index < m_groups.size() && index >= 0)
            return m_groups[index]->getRadius();
        
        return 0;
    }
    
    double SourcesManager::groupGetAzimuth(long index)
    {
        if(index < m_groups.size() && index >= 0)
            return m_groups[index]->getAzimuth();
        
        return 0;
    }
    
    double SourcesManager::groupGetElevation(long index)
    {
        if(index < m_groups.size() && index >= 0)
            return m_groups[index]->getElevation();
        
        return 0;
    }
    
    double SourcesManager::groupGetAbscissa(long index)
    {
        if(index < m_groups.size() && index >= 0)
            return m_groups[index]->getAbscissa();
        
        return 0;
    }
    
    double SourcesManager::groupGetOrdinate(long index)
    {
        if(index < m_groups.size() && index >= 0)
            return m_groups[index]->getOrdinate();
        
        return 0;
    }
    
    double SourcesManager::groupGetHeight(long index)
    {
        if(index < m_groups.size() && index >= 0)
            return m_groups[index]->getHeight();
        
        return 0;
    }
    
    double* SourcesManager::groupGetColor(long index)
    {
        if(index < m_groups.size() && index >= 0)
            return m_groups[index]->getColor();
        
        return NULL;
    }
    
    std::string SourcesManager::groupGetDescription(long index)
    {
        if(index < m_groups.size() && index >= 0)
            return m_groups[index]->getDescription();
        
        return NULL;
    }
    
    long SourcesManager::groupGetExistence(long index)
    {
        if(index < m_groups.size() && index >= 0)
        {
            return m_groups[index]->getExistence();
        }
        return 0;
    }
    
    long SourcesManager::groupGetNumberOfSources(long index)
    {
        if(index < m_groups.size() && index >= 0)
        {
            return m_groups[index]->getNumberOfSources();
        }
        return 0;
    }
    
    long SourcesManager::groupGetSourceIndex(long groupIndex, long sourceIndex)
    {
        if(groupIndex < m_groups.size() && groupIndex >= 0)
        {
            return m_groups[groupIndex]->getSourceIndex(sourceIndex);
        }
        return 0;
    }
    
    long SourcesManager::groupGetMute(long index)
    {
        if(index < m_groups.size() && index >= 0)
        {
            return m_groups[index]->getMute();
        }
        return 0;
    }
    
    bool SourcesManager::groupGetIfSourceMuted(long index)
    {
        if(index < m_groups.size() && index >= 0)
        {
            for(long i = 0; i < groupGetNumberOfSources(index); i++)
            {
                if(sourceGetMute(groupGetSourceIndex(index, i)))
                {
                    return true;
                }
            }
        }
        return false;
    }
    
    long SourcesManager::groupGetNextIndex()
    {
        if(getNumberOfGroups() != 0)
        {
            for(int i = 0; i < getMaximumIndexOfGroup()+2; i++)
            {
                if(!groupGetExistence(i))
                {
                    return i;
                }
            }
        }
        return getMaximumIndexOfGroup();
    }
    
    void SourcesManager::groupClean()
    {
        for(int i = 0; i <= getNumberOfGroups(); i++)
        {
            if(groupGetExistence(i))
            {
                for(int j = 0; j <= getNumberOfGroups(); j++)
                {
                    if (i != j && groupGetExistence(j))
                    {
                        if(groupGetNumberOfSources(i) == groupGetNumberOfSources(j))
                        {
                            int check = 0;
                            for(int k = 0; k < groupGetNumberOfSources(i); k++)
                            {
                                for(int l = 0; l < groupGetNumberOfSources(i); l++)
                                {
                                    if(groupGetSourceIndex(i, k) == groupGetSourceIndex(j, l))
                                        check++;
                                }
                            }
                            if(check == groupGetNumberOfSources(j))
                                groupRemove(j);
                        }
                    }
                }
            }
        }
        
        for(int i = 0; i < getNumberOfGroups(); i++)
        {
            if(groupGetExistence(i))
            {
                if(groupGetNumberOfSources(i) < 2)
                {
                    groupRemove(i);
                }
            }
        }
    }
    
    /************************************/
    
    void SourcesManager::copyTo(SourcesManager* sourcesManagerDestination)
    {
        double* color;
        sourcesManagerDestination->setExistence(0);
        if(getExistence() == 1)
        {
            sourcesManagerDestination->setExistence(1);
            sourcesManagerDestination->setMaximumRadius(getLimitMaximum());
            
            for(long i = 0; i <= getMaximumIndexOfSource(); i++)
            {
                if(sourceGetExistence(i) == 1)
                {
                    color = sourceGetColor(i);
                    sourcesManagerDestination->sourceSetRadius(i, sourceGetRadius(i));
                    sourcesManagerDestination->sourceSetAzimuth(i, sourceGetAzimuth(i));
                    sourcesManagerDestination->sourceSetElevation(i, sourceGetElevation(i));
                    sourcesManagerDestination->sourceSetColor(i, color[0], color[1], color[2], color[3]);
                    sourcesManagerDestination->sourceSetDescription(i, sourceGetDescription(i));
                    sourcesManagerDestination->sourceSetMute(i, sourceGetMute(i));
                }
            }
            for(long i = 0; i <= getMaximumIndexOfGroup(); i++)
            {
                if(groupGetExistence(i) == 1)
                {
                    for(long j = 0; j < groupGetNumberOfSources(i); j++)
                    {
                        sourcesManagerDestination->groupSetSource(i, groupGetSourceIndex(i, j));
                    }
                    color = groupGetColor(i);
                    sourcesManagerDestination->groupSetColor(i, color[0], color[1], color[2], color[3]);
                    sourcesManagerDestination->groupSetDescription(i, groupGetDescription(i));
                }
            }
        }
    }
    
    /************************************/
    
    SourcesManager::~SourcesManager()
    {
        m_sources.clear();
        m_groups.clear();
    }
}
