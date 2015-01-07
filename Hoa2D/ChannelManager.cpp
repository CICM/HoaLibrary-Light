/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include "ChannelManager.h"

namespace Hoa2D
{
	ChannelManager::ChannelManager(unsigned int number_of_channels)
	{
		m_defaultAzimuths = NULL;
		setNumberOfChannels(number_of_channels);
		resetAzimuth(-1);
	}
	ChannelManager::~ChannelManager()
	{
		m_channels.clear();
		if(m_defaultAzimuths)
			free(m_defaultAzimuths);
	}
	
	void ChannelManager::setNumberOfChannels(unsigned int number_of_channels)
	{
		number_of_channels = clip_min(number_of_channels, 3);
		
		if(m_channels.size() > number_of_channels)
		{
			while (m_channels.size() != number_of_channels)
				m_channels.pop_back();
		}
		else if(m_channels.size() < number_of_channels)
		{
			while (m_channels.size() != number_of_channels)
				m_channels.push_back(new Channel());
		}
		
		for (int i = 0; i < m_channels.size(); i++)
			m_channels[i]->setSelected(0);
		
		setDefaultAzimuth();
	}
	
	void ChannelManager::setDefaultAzimuth()
	{
		if(m_defaultAzimuths)
			delete m_defaultAzimuths;
		
		m_defaultAzimuths = new double[m_channels.size()];
		
		for(int i = 0; i < m_channels.size(); i++)
			m_defaultAzimuths[i] = HOA_2PI / (double)m_channels.size() * (double)i;
	}
	
	void ChannelManager::resetAzimuth(const int index)
	{
		if (!isInside(index, -1, m_channels.size()))
			return;
			
		if(index == -1)
		{
			for(int i = 0; i < m_channels.size(); i++)
				m_channels[i]->setAzimuth(m_defaultAzimuths[i]);
		}
		else
			m_channels[index]->setAzimuth(m_defaultAzimuths[index]);
	}
	
	void ChannelManager::resetWideningValue(const int index)
	{
		if (!isInside(index, -1, m_channels.size()))
			return;
		
		if(index == -1)
		{
			for (int i = 0; i < m_channels.size(); i++)
				m_channels[i]->setWideningValue(1);
		}
		else
			m_channels[index]->setWideningValue(1);
	}
	
	void ChannelManager::setAzimuthList(double* angles, long size)
	{
		for(int i = 0; i < size && i < m_channels.size(); i++)
			m_channels[i]->setAzimuth(angles[i]);
	}
	
	void ChannelManager::setAzimuth(const int index, double angle)
	{
		if (!isInside(index, -1, m_channels.size()))
			return;
		
		if(index == -1)
		{
			for(int i = 0; i< m_channels.size(); i++)
				m_channels[i]->setAzimuth(angle);
		}
		else
			m_channels[index]->setAzimuth(angle);
	}
	
	void ChannelManager::setWideningValueList(double* wideningValues, long size)
	{
		for (int i=0; i < size && i < m_channels.size(); i++)
			m_channels[i]->setWideningValue(wideningValues[i]);
	}
	
	void ChannelManager::setWideningValue(const int index, double widerValue)
	{
		if (!isInside(index, -1, m_channels.size()))
			return;
		
		if (index == -1)
		{
			for (int i=0; i < m_channels.size(); i++)
				m_channels[i]->setWideningValue(widerValue);
		}
		else
			m_channels[index]->setWideningValue(widerValue);
	}
	
	void ChannelManager::setSelected(const int index, int _selectedState)
	{
		if (!isInside(index, -1, m_channels.size()))
			return;
		
		if (index == -1)
		{
			for (int i = 0; i < m_channels.size(); i++)
				m_channels[i]->setSelected(_selectedState);
		}
		else
			m_channels[index]->setSelected(_selectedState);
	}
	
	void ChannelManager::rotateSelectedChannels(double newRadian, int channelBeingDragged, int magnet)
	{
		if (!isInside(channelBeingDragged, 0, m_channels.size()))
			return;
		
		double deltaAngle;
		double oldAngle = m_channels[channelBeingDragged]->getAzimuth();
		deltaAngle = newRadian - m_channels[channelBeingDragged]->getAzimuth();
		m_channels[channelBeingDragged]->setAzimuth(newRadian);
		
		if (magnet > 0)
		{
			if ( getClosestDefChannelDistance(channelBeingDragged) < HOA_2PI / (double) m_channels.size() *0.1 )
			{
				m_channels[channelBeingDragged]->setAzimuth(getClosestDefChannelAzimuth(channelBeingDragged));
				deltaAngle =  m_channels[channelBeingDragged]->getAzimuth() - oldAngle;
			}
		}
		
		for (int i = 0; i < m_channels.size(); i++)
		{
			if (isSelected(i) && i != channelBeingDragged)
			{
				m_channels[i]->rotateAzimuth(deltaAngle);
			}
		}
	}
	
	/* ---- Fisheye ---- */
	
	void ChannelManager::setFisheyeDestAzimuth(double azimuth)
	{
		m_fisheyeDestAzimuth = wrap_twopi(azimuth);
		m_fisheyeStep = 0;
	}
	
	void ChannelManager::setFisheyeStepWithDelta(const int index, double delta)
	{
		m_fisheyeStep = clip_minmax(m_fisheyeStep + delta, 0, 1);
		setFisheyeStepDirect(index, m_fisheyeStep);
	}
	
	void ChannelManager::setFisheyeStepDirect(const int index, double fisheyeStep)
	{
		if (!isInside(index, -2, m_channels.size()))
			return;
				
		m_fisheyeStep = clip_minmax(fisheyeStep, 0, 1);
		
		if (index == -2)
		{
			for (int i=0; i < m_channels.size(); i++)
				if (isSelected(i))
					m_channels[i]->setAzimuth(radianInterp(m_fisheyeStep, m_channels[i]->getFisheyeStartAzimuth(), m_fisheyeDestAzimuth));
		}
		else if (index == -1)
		{
			for (int i = 0; i < m_channels.size(); i++)
				m_channels[i]->setAzimuth(radianInterp(m_fisheyeStep, m_channels[i]->getFisheyeStartAzimuth(), m_fisheyeDestAzimuth));
		}
		else
			m_channels[index]->setAzimuth(radianInterp(m_fisheyeStep, m_channels[index]->getFisheyeStartAzimuth(), m_fisheyeDestAzimuth));
	}
	
	void ChannelManager::setFisheyeStartAzimuth(const int index)
	{
		if (!isInside(index, -2, m_channels.size()))
			return;
		
		if (index == -2)
		{
			for (int i=0; i < m_channels.size(); i++)
				if (isSelected(i))
					m_channels[i]->setFisheyeStartAzimuth();
		}
		else if (index == -1)
		{
			for (int i=0; i < m_channels.size(); i++)
				m_channels[i]->setFisheyeStartAzimuth();
		}
		else
			m_channels[index]->setFisheyeStartAzimuth();
	}
	
	void ChannelManager::setFisheyeStartAzimuth(const int index, double azimuth)
	{
		if (!isInside(index, -2, m_channels.size()))
			return;
		
		if (index == -2) {
			for (int i=0; i < m_channels.size(); i++)
				if (isSelected(i))
					m_channels[i]->setFisheyeStartAzimuth(azimuth);
		}
		else if (index == -1)
			for (int i=0; i<m_channels.size(); i++)
				m_channels[i]->setFisheyeStartAzimuth(azimuth);
		else
			m_channels[index]->setFisheyeStartAzimuth(azimuth);
	}
	
	/* --- Utility --- */
	
	double ChannelManager::radianInterp(double step, double startRad, double endRad)
	{
		double start = wrap_twopi(startRad);
		double end   = wrap_twopi(endRad);
		
		// anti-clockwise
		
		if ( wrap_twopi(end - start) <= HOA_PI )
		{
			if (end - start >= 0)
				return wrap_twopi( start + step*(end - start) );
			else
				return wrap_twopi( start + step*( (end + HOA_2PI) - start) );
		}
		
		// clockwise
		
		else
		{
			if (end - start <= 0)
				return wrap_twopi( start + step*(end - start) );
			else
				return wrap_twopi( start - step*( (start + HOA_2PI) - end) );
		}
	}
	
	double ChannelManager::getClosestDefChannelAzimuth(const int index)
	{
		return getClosestDefChannelAzimuth(getAzimuth(index));
	}
	
	double ChannelManager::getClosestDefChannelAzimuth(double azimuth)
	{
		double angle = wrap_twopi(azimuth);
		double distance = HOA_2PI;
		double tempDistance, closestAngle;
		closestAngle = angle;
		
		for (int i = 0; i < m_channels.size(); i++)
		{
			tempDistance = radianClosestDistance(angle, m_defaultAzimuths[i]);
			if (tempDistance < distance) {
				distance = tempDistance;
				closestAngle = m_defaultAzimuths[i];
			}
		}
		
		return closestAngle;
	}
	
	double ChannelManager::getClosestDefChannelDistance(const int index)
	{
		return getClosestDefChannelDistance(getAzimuth(index));
	}
	
	double ChannelManager::getClosestDefChannelDistance(double azimuth)
	{
		double angle = wrap_twopi(azimuth);
		double distance = HOA_2PI;
		
		for (int i = 0; i < m_channels.size(); i++)
			distance = min(distance, radianClosestDistance(angle, m_defaultAzimuths[i]));
		
		return distance;
	}
	
	void ChannelManager::setAzimuthToClosestDefChannelAzimuth(const int index)
	{
		if (!isInside(index, 0, m_channels.size())) return;
		m_channels[index]->setAzimuth(getClosestDefChannelAzimuth(index));
	}
	
	long ChannelManager::getNumberOfSelectedChannels()
	{
		long numberOfSelectedChannels = 0;
		for (int i = 0; i < m_channels.size(); i++)
			if (m_channels[i]->isSelected()) numberOfSelectedChannels++;
		return numberOfSelectedChannels;
	}
}
