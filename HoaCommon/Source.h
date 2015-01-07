/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_SOURCE
#define DEF_SOURCE

#include "../Hoa.h"

namespace HoaCommon
{
	//! The source
    /** The source store and manage source informations like the position, color, mute...
     */
	class Source
	{
	private:
		
		double					m_radius;
		double					m_azimuth;
        double                  m_elevation;
		double*					m_color;
		std::string             m_description;
		bool                    m_exist;
		std::vector <long>      m_groups;
		double                  m_maximum_radius;
		bool                    m_mute;
		
	public:
		
		//! The source constructor.
		/**	The source constructor allocates and initialize the member values for a source.
		 
			@param     existence		The existence state of the source.
			@param     radius			The radius of the source.
			@param     azimuth			The azimuth of the source.
            @param     elevation		The elevation of the source.
		 */
		Source(bool existence = true, double radius = 0., double azimuth = 0., double elevation = 0.);
		
		//! The source destructor.
        /**	The source destructor free the memory.
         */
		~Source();
		
		//! Set the existence state of the source.
		/**
			@param     state		The existence state of the source.
			@see getExistence
         */
		void setExistence(bool state);
		
        //! Set the position of the source with polar coordinates.
		/**
            @param     radius			The radius of the source.
            @param     azimuth			The azimuth of the source.
            @see setCoordinatesCartesian
         */
		void setCoordinatesPolar(double radius, double azimuth);
        
		//! Set the position of the source with polar coordinates.
		/**
			@param     radius			The radius of the source.
			@param     azimuth			The azimuth of the source.
            @param     elevation        The elevation of the source.
			@see setCoordinatesCartesian
         */
		void setCoordinatesPolar(double radius, double azimuth, double elevation);
		
		//! Set the radius of the source.
		/**
			@param     radius			The radius of the source.
			@see getRadius
         */
		void setRadius(double radius);
		
		//! Set the azimuth of the source.
		/**
			@param     azimuth			The azimuth of the source.
         */
		void setAzimuth(double azimuth);
        
        //! Set the elevation of the source.
		/**
         @param     elevation			The elevation of the source.
         */
		void setElevation(double elevation);
		
		//! Set the position of the source with cartesians coordinates.
		/**
			@param     abscissa		The abscissa of the source.
			@param     ordinate		The ordinate of the source.
         */
		void setCoordinatesCartesian(double abscissa, double ordinate);
        
        //! Set the position of the source with cartesians coordinates.
		/**
            @param     abscissa		The abscissa of the source.
            @param     ordinate		The ordinate of the source.
            @param     height		The height of the source.
         */
		void setCoordinatesCartesian(double abscissa, double ordinate, double height);
		
		//! Set the abscissa of the source.
		/**
			@param     abscissa		The abscissa of the source.
         */
		void setAbscissa(double abscissa);
		
		//! Set the ordinate of the source.
		/**
			@param     ordinate		The ordinate of the source.
         */
		void setOrdinate(double ordinate);
        
        //! Set the height of the source.
		/**
         @param     height		The height of the source.
         */
		void setHeight(double height);
		
		//! Set the color of the source.
		/**
		 * @param     red		The red component of the color.
		 * @param     green		The green component of the color
		 * @param     blue		The blue component of the color
		 * @param     alpha		The alpha component of the color
         */
		void setColor(double red, double green, double blue, double alpha);
		
		//! Set the description of the source.
		/**
		 * @param     description		The text description of the source.
         */
		void setDescription(std::string description);
		
		//! Add source to an indexed group.
		/**
			@param     groupIndex		The index of the group.
         */
		void setGroup(long groupIndex);
		
		//! Remove source from an indexed group.
		/**
			@param     groupIndex		The index of the group.
         */
		void removeGroup(long groupIndex);
		
		//! Set the maximum radius of the source
		/**
			@param     limitValue		The radius limit value.
         */
		void setMaximumRadius(double limitValue);
		
		//! Set the mute state of the source
		/**
			@param     state		The mute state of the source.
			@see getMute
         */
		void setMute(bool state);
		
		//! Get the existence state of the source.
		/**
			@return		The existence state of the source.
			@see setExistence
         */
		bool			getExistence()		const {return m_exist;}
		
		//! Get the radius of the source.
		/**
			@return		The radius of the source.
			@see setRadius, setCoordinatesPolar
         */
		double			getRadius()			const {return m_radius;}
		
		//! Get the azimuth of the source.
		/**
			@return		The azimuth of the source.
			@see setAzimuth, setCoordinatesPolar
         */
		double			getAzimuth()		const {return m_azimuth;}
        
        //! Get the elevation of the source.
		/**
            @return		The elevation of the source.
            @see setElevation, setCoordinatesPolar
         */
		double			getElevation()		const {return m_elevation;}
		
		//! Get the abscissa of the source.
		/**
			@return		The abscissa of the source.
			@see setAbscissa, setCoordinatesCartesian
         */
		double			getAbscissa()		const {return abscissa(m_radius, m_azimuth, m_elevation);}
		
		
		//! Get the ordinate of the source.
		/**
			@return		The ordinate of the source.
			@see setOrdinate, setCoordinatesCartesian
         */
		double			getOrdinate()		const {return ordinate(m_radius, m_azimuth, m_elevation);}
        
        //! Get the height of the source.
		/**
         @return		The height of the source.
         @see setHeight, setCoordinatesCartesian
         */
		double			getHeight()		const {return height(m_radius, m_azimuth, m_elevation);}
		
		//! Get the color of the source.
		/**
			@return		The rgba color of the source as an array of 4 double.
			@see setColor
         */
		double*			getColor()			const {return m_color;}
		
		//! Get the description of the source.
		/**
			@return		The description of the source.
			@see setDescription
         */
		std::string		getDescription()	const {return m_description;}
		
		//! Get the number of group the source is owned by.
		/**
			@return		The number of group.
			@see setDescription
         */
		long			getNumberOfGroups() const {return m_groups.size();}
		
		//! Get the the group index the source is owned by at a particular index.
		/**
			@param			index			The index of the group
			@return		The group index.
         */
		long			getGroupIndex(long index);
		
		//! Determine if the source is owned by a particular group
		/**
			@param			groupIndex		The index of the group
			@return		true if the source is in the group, false otherwise.
         */
		bool			isOwnedByGroup(long groupIndex);
		
		//! Get the mute state of the source
		/**
			@return		The mute state of the source.
			@see setMute
         */
		bool			getMute()			const {return m_mute;}
	};
}
#endif