/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef __DEF_HOA_2D_CHANNELMANAGER__
#define __DEF_HOA_2D_CHANNELMANAGER__

#include "../Hoa.h"

namespace Hoa2D
{
	//! The channel manager
    /** The channel manager should be used to manage a group of channel.
     */
	class ChannelManager
	{
	private:

		class Channel
		{
		private:
			
			double  m_azimuth;
			double  m_wideningValue;
			bool    m_isSelected;
			double  m_fisheyeStartAzimuth;
			double  m_fisheyeEndAzimuth;
			
		public:

			Channel(double azimuth = 0.0f, double wideningValue = 1.0f)
			{
				setAzimuth(azimuth);
				m_fisheyeEndAzimuth = m_fisheyeStartAzimuth = m_azimuth;
				setWideningValue(wideningValue);
				m_isSelected = false;
			}
			~Channel() {}
			
			bool   isSelected()			const {return m_isSelected;}
			double getWideningValue()		const {return m_wideningValue;}
			double getFisheyeStartAzimuth() const { return m_fisheyeStartAzimuth; }
			double getAbscissa()		const {return abscissa(1.0f, m_azimuth);}
			double getOrdinate()		const {return ordinate(1.0f, m_azimuth);}
			double getAzimuth()			const {return m_azimuth;}
			
			void setAzimuth(const double angle) { m_azimuth = wrap_twopi(angle); }
			void rotateAzimuth(const double deltaAzimuth) { m_azimuth = wrap_twopi(m_azimuth + deltaAzimuth); }
			void setWideningValue(const double widerValue) { m_wideningValue = clip_minmax(widerValue, 0, 1); }
			void setFisheyeStartAzimuth() { m_fisheyeStartAzimuth = m_azimuth; }
			void setFisheyeStartAzimuth(const double radian) { m_fisheyeStartAzimuth = wrap_twopi(radian); }
			void setFisheyeEndAngle() { m_fisheyeEndAzimuth = m_azimuth; }
			void setFisheyeEndAngle(const double radian) { m_fisheyeEndAzimuth = wrap_twopi(radian); }
			void setSelected(const int state) { m_isSelected = (state == -1) ? !m_isSelected : (state != 0); }
		};
		
		std::vector <Channel*> m_channels;
		double* m_defaultAzimuths;
		double  m_fisheyeStep;
		double  m_fisheyeDestAzimuth;

		double radianInterp(double step, double startRad, double endRad);
		void setDefaultAzimuth();
		
	public:
		
		//! The channel manager constructor.
		/**	The channel manager constructor allocates and initialize the member values for each channel.
		 
		 @param     number_of_channels		The number of channels you need to manage, default is 1.
		 */
		ChannelManager(unsigned int number_of_channels = 1);
		
		//! The channel manager destructor.
        /**	The channel manager destructor free the memory.
         */
		~ChannelManager();
		
		//! Set the number of channels you want to manage.
		/**
		 @param     number_of_channels		The number of channels you need to manage.
         */
		void setNumberOfChannels(unsigned int number_of_channels);
		
		//! Set one or each channel azimuth value.
		/**
		 If you want to reset all of the channels azimuth value, pass -1 in the index parameter.
         
		 @param     index		The index of the channel (pass -1 to reset all of them).
		 @param     azimuth		The new azimuth value in radians.
		 @see resetAzimuth, setAzimuthList
         */
		void setAzimuth(const int index, double azimuth);
		
		//! Set several channel azimuth with a list.
		/**
		 @param     azimuths	An array of azimuth values in radians.
		 @param     size		The size of the array.
		 @see setAzimuth, setWideningValue
         */
		void setAzimuthList(double* azimuths, long size);
		
		//! Reset one or each channel azimuth value.
		/**	Reset one or each channel to the default azimuth value.
		 If you want to reset all of the channels azimuth value, pass -1 in the index parameter.
         
		 @param     index		The index of the channel (pass -1 to reset all of them).
		 @see setAzimuth, setAzimuthList
         */
		void resetAzimuth(const int index = -1);
		
		//! Set one or each channel wideningValue value.
		/**
		 If you want to reset all of the channels azimuth value, pass -1 in the index parameter.
         
		 @param     index				The index of the channel (pass -1 to reset all of them).
		 @param     wideningValue		The new widening value (between 0 and 1).
		 @see resetWideningValue, setWideningValueList
         */
		void setWideningValue(const int index, const double wideningValue);
		
		//! Set several channel wideningValue value with a list.
		/**
		 @param     wideningValues	An array of widening values (between 0 and 1).
		 @param     size			The size of the array.
		 @see setWideningValue, resetWideningValue
         */
		void setWideningValueList(double* wideningValues, long size);
		
		//! Reset one or each channel wideningValue to 1.
		/**
		 If you want to reset all of the channels wideningValue value, pass -1 in the index parameter.
         
		 @param     index		The index of the channel (pass -1 to reset all of them).
		 @see setWideningValue, setWideningValueList
         */
		void resetWideningValue(const int index = -1);
		
		//! Set the selected state of one or each channel.
		/**
		 If you want to reset all of the channels selected state, pass -1 in the index parameter.
         
		 @param     index				The index of the channel (pass -1 to reset all of them).
		 @param     selectedState		The selected state, pass 0 to unselect, 1 to select and -1 to toggle the state.
		 @see getNumberOfSelectedChannels
         */
		void setSelected(const int index, int selectedState);
		
		//! Rotate all of the currently selected channels relative to one channel azimuth rotation.
		/**
		 @param     newRadian				The new azimuth of the channel which is being dragged
		 @param     channelBeingDragged		The index of the channel which is being dragged.
		 @param     magnet					Pass 1 to magnetize the channel to one of the default azimuth positions
         */
		void rotateSelectedChannels(double newRadian, int channelBeingDragged, int magnet = 0);
		
		//! Set fisheye destination azimuth.
		/**
		 @param     azimuth		The azimuth value of the fisheye destination in radians.
		 @see getFisheyeDestAzimuth, setFisheyeStartAzimuth
         */
		void setFisheyeDestAzimuth(double azimuth);
		
		//! Set the fisheye start azimuth of one or several channel with current azimuth.
		/**
		 @param     index				The index of the channel if greather than -1, pass -2 to affect all selected channels, pass -1 to affect all channels
		 @see
         */
		void setFisheyeStartAzimuth(const int index);
		
		//! Set the fisheye start azimuth of one or several channel.
		/**
		 @param     index				The index of the channel if greather than -1, pass -2 to affect all selected channels, pass -1 to affect all channels
         */
		void setFisheyeStartAzimuth(const int index, double radian);
		
		//! Increment or decrement the fisheye step by a value.
		/**
		 @param     index				The index of the channel if greather than -1, pass -2 to affect all selected channels, pass -1 to affect all channels
		 @see setFisheyeStepDirect
         */
		void setFisheyeStepWithDelta(const int index, double delta);
		
		//! Set the fisheye step directly.
		/**
		 @param     index				The index of the channel if greather than -1, pass -2 to affect all selected channels, pass -1 to affect all channels
		 @param		fisheyeStep			The fisheye step (between 0 and 1).
		 @see setFisheyeStepWithDelta
         */
		void setFisheyeStepDirect(const int index, double fisheyeStep);
		
		//! Set channel azimuth value to the closest default azimuth position value.
		/**
		 @param     index		The index of the channel.
		 @see resetAzimuth, setAzimuth
         */
		void setAzimuthToClosestDefChannelAzimuth(const int index);
		
		//! Retrieve the number of selected channels.
		/**
		 @return     The number of selected channels.
		 @see setSelected
         */
		long          getNumberOfSelectedChannels();
		
		//! Retrieve the closest default azimuth position of a channel.
		/**
		 @param     index		The index of the channel.
		 @return    closest default azimuth position of the channel.
		 @see getClosestDefChannelDistance
         */
		double        getClosestDefChannelAzimuth(const int index);
		
		//! Retrieve the closest default azimuth position relative to an azimuth value.
		/**
		 @param     azimuth		The relative azimuth.
		 @return    The closest default azimuth position of the channel relative to the azimuth value you passed in.
		 @see getClosestDefChannelDistance
         */
		double        getClosestDefChannelAzimuth(double azimuth);
		
		//! Retrieve the distance to the closest default azimuth position relative to the azimuth of a channel.
		/**
		 @param     index		The index of the channel.
		 @return    The distance to the closest default azimuth position.
		 @see getClosestDefChannelAzimuth
         */
		double        getClosestDefChannelDistance(const int index);
		
		//! Retrieve the distance to the closest default azimuth position relative to an azimuth value.
		/**
		 @param     azimuth		The relative azimuth.
		 @return    The distance to the closest default azimuth position of the channel relative to the azimuth value you passed in.
		 @see getClosestDefChannelAzimuth
         */
		double        getClosestDefChannelDistance(double azimuth);
		
		//! Retrieve the fisheye destination azimuth value.
		/**
		 @return    The fisheye destination azimuth value.
		 @see setFisheyeDestAzimuth
         */
		inline double getFisheyeDestAzimuth()			const { return m_fisheyeDestAzimuth; }
		
		//! Retrieve the number of channels currently managed.
		/**
		 @return    The number of channels currently managed.
		 @see setNumberOfChannels
         */
		inline long   getNumberOfChannels()					const { return m_channels.size(); }
		
		//! Retrieve the current fisheye step.
		/**
		 @return    The current fisheye step.
		 @see setFisheyeStepDirect, setFisheyeStepWithDelta
         */
		inline double getFisheyeStep()					const { return m_fisheyeStep; }
		
		//! Retrieve the selection state of a channel.
		/**
		 @param     index		The index of the channel.
		 @return    The selection state of the channel.
		 @see setSelected
         */
		inline bool   isSelected(const int index)	const
		{
			if (!isInside(index, 0, m_channels.size()))
				return 0;
			return m_channels[index]->isSelected();
		}
		
		//! Retrieve the wideningValue value of a channel.
		/**
		 @param     index		The index of the channel.
		 @return    The wideningValue value of the channel.
		 @see setWideningValue
         */
		inline double getWideningValue(const int index)	const
		{
			if (!isInside(index, 0, m_channels.size()))
				return 0;
			return m_channels[index]->getWideningValue();
		}
		
		//! Retrieve the fisheye start azimuth value of a channel.
		/**
		 @param     index		The index of the channel.
		 @return    The fisheye start azimuth value of the channel.
		 @see setFisheyeStartAzimuth
         */
		inline double getFisheyeStartAzimuth(const int index) const
		{
			if (!isInside(index, 0, m_channels.size()))
				return 0;
			return m_channels[index]->getFisheyeStartAzimuth();
		}
		
		//! Retrieve the abscissa of a channel.
		/**
		 @param     index		The index of the channel.
		 @return    The abscissa of the channel.
         */
		inline double getAbscissa(const int index)		const
		{
			if (!isInside(index, 0, m_channels.size()))
				return 0;
			return m_channels[index]->getAbscissa();
		}
		
		//! Retrieve the ordinate of a channel.
		/**
		 @param     index		The index of the channel.
		 @return    The ordinate of the channel.
         */
		inline double getOrdinate(const int index)		const
		{
			if (!isInside(index, 0, m_channels.size()))
				return 0;
			return m_channels[index]->getOrdinate();
		}
		
		//! Retrieve the azimuth of a channel.
		/**
		 @param     index		The index of the channel.
		 @return    The azimuth of the channel.
		 @see setAzimuth
         */
		inline double getAzimuth(const int index)		const
		{
			if (!isInside(index, 0, m_channels.size()))
				return 0;
			return m_channels[index]->getAzimuth();
		}
	};
}

#endif
