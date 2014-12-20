/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_SOURCESGROUP
#define DEF_SOURCESGROUP

#include "Source.h"
#include "SourcesManager.h"

namespace HoaCommon
{
	class SourcesManager;
	
	//! The sources group
    /** The SourcesGroup should be used to store and manage multiple Source
     */
	class SourcesGroup
	{
		
	private:
		
		std::vector <long>      m_sources;
		std::string             m_description;
		long                    m_exist;
		double*					m_color;
		double					m_centroid_x;
		double					m_centroid_y;
        double					m_centroid_z;
		SourcesManager*         m_source_manager;
		double                  m_maximum_radius;
		long                    m_mute;
		
		void computeCentroid();
		void shiftPolar(double radius, double azimuth);
        void shiftPolar(double radius, double azimuth, double elevation);
		void shiftRadius(double radius);
		void shiftAzimuth(double azimuth);
        void shiftElevation(double elevation);
		void shiftCartesian(double abscissa, double ordinate);
        void shiftCartesian(double abscissa, double ordinate, double height);
		void shiftAbscissa(double abscissa);
		void shiftOrdinate(double ordinate);
        void shiftHeight(double ordinate);
		
	public:
		
		//! The source group constructor.
		/**	The source group constructor allocates and initialize the member values for a source group.
		 
			@param     sourcesManager		A SourceManager object's pointer.
			@param     existence			The existence state of the source.
		 */
		SourcesGroup(SourcesManager* sourcesManager, bool existence);
		
		//! The source group destructor.
        /**	The source group destructor free the memory.
         */
		~SourcesGroup();
		
		//! Set the existence state of the group.
		/**
			@param     state		The existence state of the group.
			@see getExistence
         */
		void setExistence(bool state);
		
		//! Set the position of the group with polar coordinates.
		/**
			@param     radius			The radius of the group.
			@param     azimuth			The azimuth of the group.
			@see setRelativeCoordinatesPolar, setCoordinatesCartesian
         */
		void setCoordinatesPolar(double radius, double azimuth);
        
        //! Set the position of the group with polar coordinates.
		/**
         @param     radius			The radius of the group.
         @param     azimuth			The azimuth of the group.
         @param     elevation			The elevation of the group.
         @see setRelativeCoordinatesPolar, setCoordinatesCartesian
         */
		void setCoordinatesPolar(double radius, double azimuth, double elevation);
		
		//! Set the radius of the group.
		/**
			@param     radius			The radius of the group.
			@see getRadius
         */
		void setRadius(double radius);
		
		//! Set the azimuth of the group.
		/**
		 @param     azimuth			The azimuth of the group.
         */
		void setAzimuth(double azimuth);
        
        //! Set the elevation of the group.
		/**
		 @param     elevation			The elevation of the group.
         */
		void setElevation(double elevation);
		
		//! Set the position of the group with cartesians coordinates.
		/**
			@param     abscissa		The abscissa of the group.
			@param     ordinate		The ordinate of the group.
         */
		void setCoordinatesCartesian(double abscissa, double ordinate);
        
        //! Set the position of the group with cartesians coordinates.
		/**
         @param     abscissa		The abscissa of the group.
         @param     ordinate		The ordinate of the group.
         @param     height		The height of the group.
         */
		void setCoordinatesCartesian(double abscissa, double ordinate, double height);
		
		//! Set the abscissa of the group.
		/**
			@param     abscissa		The abscissa of the group.
         */
		void setAbscissa(double abscissa);
		
		//! Set the ordinate of the group.
		/**
			@param     ordinate		The ordinate of the group.
         */
		void setOrdinate(double ordinate);
        
        //! Set the height of the group.
		/**
         @param     height		The height of the group.
         */
		void setHeight(double height);
		
		//! Set the color of the group.
		/**
			@param     red			The red component of the color.
			@param     green		The green component of the color
			@param     blue			The blue component of the color
			@param     alpha		The alpha component of the color
         */
		void setColor(double red, double green, double blue, double alpha);
		
		//! Set the description of the group.
		/**
			@param     description		The text description of the group.
         */
		void setDescription(std::string description);
		
		//! Set the maximum radius of the group
		/**
			@param     limitValue		The radius limit value.
         */
		void setMaximumRadius(double limitValue);
		
		//! Store a new source in this group
		/**
			@param     sourceIndex		The index of the source to store.
         */
		void addSource(long sourceIndex);
		
		//! Remove a source from this group
		/**
			@param     sourceIndex		The index of the source to store.
         */
		void removeSource(long sourceIndex);
		
		//! Notify the group that a source has moved
		/** You need to call this function whenever a source has moved to update group information.
			@param     sourceIndex		The index of the source to store.
         */
		void sourceHasMoved();
		
		//! Set the mute state of the group
		/**
			@param     state		The mute state of the group.
			@see getMute
         */
		void setMute(long aValue);
		
		//! Set the position of the group with relative polar coordinates.
		/**
			@param     radius			The relative radius value.
			@param     azimuth			The relative azimuth value.
			@see setCoordinatesPolar, setRelativeRadius
         */
		void setRelativeCoordinatesPolar(double radius, double azimuth);
        
        //! Set the position of the group with relative polar coordinates.
		/**
         @param     radius			The relative radius value.
         @param     azimuth			The relative azimuth value.
         @param     elevation		The relative elevation value.
         @see setCoordinatesPolar, setRelativeRadius
         */
		void setRelativeCoordinatesPolar(double radius, double azimuth, double elevation);
		
		//! Set the radius of the group with a relative value.
		/**
			@param     radius			The relative radius value.
			@see setCoordinatesPolar, setRelativeRadius
         */
		void setRelativeRadius(double radius);
		
		//! Set the azimuth of the group with a relative value.
		/**
			@param     azimuth			The relative azimuth value.
			@see setCoordinatesPolar, setRelativeRadius
         */
		void setRelativeAzimuth(double azimuth);
        
        //! Set the elevation of the group with a relative value.
		/**
         @param     elevation			The relative elevation value.
         @see setCoordinatesPolar, setRelativeRadius
         */
		void setRelativeElevation(double elevation);
		
		//! Get the existence state of the group.
		/**
			@return		The existence state of the group.
			@see setExistence
         */
		bool			getExistence()		const {return m_exist;}
		
		//! Get the radius of the group.
		/**
			@return		The radius of the group.
			@see setRadius, setCoordinatesPolar
         */
		double			getRadius()			const { return radius(m_centroid_x, m_centroid_y, m_centroid_z);}
		
		//! Get the azimuth of the group.
		/**
			@return		The azimuth of the group.
			@see setAzimuth, setCoordinatesPolar
         */
		double			getAzimuth()		const {return azimuth(m_centroid_x, m_centroid_y, m_centroid_z) + HOA_PI2;}
        
        //! Get the elevation of the group.
		/**
         @return		The elevation of the group.
         @see setAzimuth, setCoordinatesPolar
         */
		double			getElevation()		const {return elevation(m_centroid_x, m_centroid_y, m_centroid_z);}
		
		//! Get the abscissa of the group.
		/**
			@return		The abscissa of the group.
			@see setAbscissa, setCoordinatesCartesian
         */
		double			getAbscissa()		const {return m_centroid_x;}
		
		//! Get the ordinate of the group.
		/**
			@return		The ordinate of the group.
			@see setOrdinate, setCoordinatesCartesian
         */
		double			getOrdinate()		const {return m_centroid_y;}
        
        //! Get the height of the group.
		/**
         @return		The height of the group.
         @see setOrdinate, setCoordinatesCartesian
         */
		double			getHeight()		const {return m_centroid_z;}
		
		//! Get the color of the group.
		/**
			@return		The rgba color of the group as an array of 4 double.
			@see setColor
         */
		double*			getColor()			const {return m_color;}
		
		//! Get the description of the group.
		/**
			@return		The description of the group.
			@see setDescription
         */
		std::string		getDescription()	const {return m_description;}
		
		//! Get the number of sources owned by this group.
		/**
			@return		The number of sources.
			@see setDescription
         */
		long			getNumberOfSources() const {return m_sources.size();}
		
		//! Get the mute state of the group
		/**
			@return		The mute state of the group.
			@see setMute
         */
		bool			getMute()			const {return m_mute;}
		
		//! Get the the index of a source stored at a particular index by the group.
		/**
			@param			index			The index of the source.
			@return			The index of the source if it exists, -1 otherwise.
         */
		long   getSourceIndex(long index)	const
		{
			if(index < m_sources.size() && index >= 0)
				return m_sources[index];
			else
				return -1;
		}
	};
}
#endif
