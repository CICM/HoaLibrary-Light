/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_SOURCE_LIGHT
#define DEF_HOA_SOURCE_LIGHT

#include "Math.hpp"

namespace hoa
{
    class SourcesManager;
    class SourcesGroup;
    
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
		double			getAbscissa()		const {return Math<double>::abscissa(m_radius, m_azimuth, m_elevation);}
		
		
		//! Get the ordinate of the source.
		/**
			@return		The ordinate of the source.
			@see setOrdinate, setCoordinatesCartesian
         */
		double			getOrdinate()		const {return Math<double>::ordinate(m_radius, m_azimuth, m_elevation);}
        
        //! Get the height of the source.
		/**
         @return		The height of the source.
         @see setHeight, setCoordinatesCartesian
         */
		double			getHeight()		const {return Math<double>::height(m_radius, m_azimuth, m_elevation);}
		
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
        double			getRadius()			const { return Math<double>::radius(m_centroid_x, m_centroid_y, m_centroid_z);}
        
        //! Get the azimuth of the group.
        /**
         @return		The azimuth of the group.
         @see setAzimuth, setCoordinatesPolar
         */
        double			getAzimuth()		const {return Math<double>::azimuth(m_centroid_x, m_centroid_y, m_centroid_z) + HOA_PI2;}
        
        //! Get the elevation of the group.
        /**
         @return		The elevation of the group.
         @see setAzimuth, setCoordinatesPolar
         */
        double			getElevation()		const {return Math<double>::elevation(m_centroid_x, m_centroid_y, m_centroid_z);}
        
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
    
    class SourcesManager
    {
        
    private:
        
        double                      m_maximum_radius;
        std::vector <Source*>       m_sources;
        std::vector <SourcesGroup*> m_groups;
        bool                        m_exist;
        double                      m_zoom;
        
        void checkMute();
        
    public:
        
        //! The source manager constructor.
        /**	The source manager constructor allocates and initialize the member values.
         *
         * @param     maximumRadius		The maximum radius the sources or groups in the source manager could have
         * @param     existence			The existence state of the source manager.
         */
        SourcesManager(double maximumRadius = 1., bool existence = true);
        
        //! The source manager destructor free the memory.
        ~SourcesManager();
        
        //! Make a copy of this SourcesManager into an other
        void copyTo(SourcesManager* sourcesManagerDestination);
        
        //! Clear all the sources and groups.
        void clearAll();
        
        //! Set the maximum radius the sources and groups can have
        /**
         * @param     limitValue		The radius limit value.
         */
        void setMaximumRadius(double limitValue);
        
        //! Set the existence state of the sources manager
        /** If the existence state is false it will delete all sources and groups already stored.
         *
         * @param     state		The existence state.
         */
        void setExistence(bool state);
        
        //! Set the zoom factor.
        /** Set the zoom factor between 0 and 1.
         
         @param     zoom		The zoom factor.
         */
        void setZoom(double zoom);
        
        //! Get the maximum index of sources.
        /**
         @return		The maximum index of source.
         */
        long getMaximumIndexOfSource();
        
        //! Get the number of sources actually managed.
        /**
         @return		The number of sources.
         */
        long getNumberOfSources();
        
        //! Get the maximum index of the sources actually managed.
        /**
         @return		The maximum index.
         */
        long getMaximumIndexOfGroup();
        
        //! Get the number of groups actually managed.
        /**
         @return		The number of group.
         */
        long getNumberOfGroups();
        
        //! Get the maximum radius of the sources and groups.
        /**
         @return		The maximum radius.
         */
        double getLimitMaximum();
        
        //! Get the existence state of the source manager.
        /**
         @return		The existence state.
         */
        bool getExistence();
        
        //! Get the zoom factor value.
        /**
         @return		The zoom factor value.
         */
        double getZoom();
        
        /* ------------------------------------------------------------------------ */
        /* ------------------------------- Sources -------------------------------- */
        /* ------------------------------------------------------------------------ */
        
        //! Add a new source with polar coordinates.
        /**
         * @param     radius			The radius of the source.
         * @param     azimuth			The azimuth of the source.
         * @see sourceNewCartesian
         */
        void sourceNewPolar(double radius, double azimuth);
        
        //! Add a new source with polar coordinates.
        /**
         * @param     radius			The radius of the source.
         * @param     azimuth			The azimuth of the source.
         * @param     elevation			The elevation of the source.
         * @see sourceNewCartesian
         */
        void sourceNewPolar(double radius, double azimuth, double elevation);
        
        //! Add a new source with cartesian coordinates.
        /**
         * @param     abscissa			The abscissa of the source.
         * @param     ordinate			The ordinate of the source.
         * @see sourceNewPolar
         */
        void sourceNewCartesian(double abscissa, double ordinate);
        
        //! Add a new source with cartesian coordinates.
        /**
         * @param     abscissa			The abscissa of the source.
         * @param     ordinate			The ordinate of the source.
         * @param     height			The height of the source.
         * @see sourceNewPolar
         */
        void sourceNewCartesian(double abscissa, double ordinate, double height);
        
        //! Set position of a source with polar coordinates.
        /**
         * @param     index				The index of the source.
         * @param     radius			The radius of the source.
         * @param     azimuth			The azimuth of the source.
         * @see sourceSetRadius, sourceSetAzimuth, sourceSetCartesian
         */
        void sourceSetPolar(long index, double radius, double azimuth);
        
        //! Set position of a source with polar coordinates.
        /**
         * @param     index				The index of the source.
         * @param     radius			The radius of the source.
         * @param     azimuth			The azimuth of the source.
         * @param     elevation			The elevation of the source.
         * @see sourceSetRadius, sourceSetAzimuth, sourceSetCartesian
         */
        void sourceSetPolar(long index, double radius, double azimuth, double elevation);
        
        //! Set radius of a source.
        /**
         * @param     index				The index of the source.
         * @param     radius			The radius of the source.
         * @see sourceSetPolar, sourceSetAzimuth
         */
        void sourceSetRadius(long index, double radius);
        
        //! Set azimuth of a source.
        /**
         * @param     index				The index of the source.
         * @param     azimuth			The azimuth of the source.
         * @see sourceSetRadius, sourceSetAzimuth
         */
        void sourceSetAzimuth(long index, double azimuth);
        
        //! Set elevation of a source.
        /**
         * @param     index				The index of the source.
         * @param     elevation			The elevation of the source.
         * @see sourceSetRadius, sourceSetAzimuth
         */
        void sourceSetElevation(long index, double elevation);
        
        //! Set position of a source with cartesian coordinates.
        /**
         * @param     index				The index of the source.
         * @param     abscissa			The abscissa of the source.
         * @param     ordinate			The ordinate of the source.
         * @see sourceSetAbscissa, sourceSetOrdinate, sourceSetPolar
         */
        void sourceSetCartesian(long index, double abscissa, double ordinate);
        
        //! Set position of a source with cartesian coordinates.
        /**
         * @param     index				The index of the source.
         * @param     abscissa			The abscissa of the source.
         * @param     ordinate			The ordinate of the source.
         * @param     height			The height of the source.
         * @see sourceSetAbscissa, sourceSetOrdinate, sourceSetPolar
         */
        void sourceSetCartesian(long index, double abscissa, double ordinate, double height);
        
        //! Set abscissa of a source.
        /**
         * @param     index				The index of the source.
         * @param     abscissa			The abscissa of the source.
         * @see sourceSetOrdinate, sourceSetHeight
         */
        void sourceSetAbscissa(long index, double abscissa);
        
        //! Set ordinate of a source.
        /**
         * @param     index				The index of the source.
         * @param     ordinate			The ordinate of the source.
         * @see sourceSetAbscissa, sourceSetHeight
         */
        void sourceSetOrdinate(long index, double ordinate);
        
        //! Set height of a source.
        /**
         * @param     index				The index of the source.
         * @param     height			The height of the source.
         * @see sourceSetOrdinate, sourceSetAbscissa
         */
        void sourceSetHeight(long index, double height);
        
        //! Set the rgba color of a source.
        /** All values are clipped between 0 and 1.
         * @param     index				The index of the source.
         * @param     red				The red component of the color.
         * @param     green				The green component of the color
         * @param     blue				The blue component of the color
         * @param     alpha				The alpha component of the color
         */
        void sourceSetColor(long index, double red, double green, double blue, double alpha);
        
        //! Add a description to a given source.
        /**
         * @param     index				The index of the source.
         * @param     description		The text description of the source.
         */
        void sourceSetDescription(long index, std::string description);
        
        //! Remove a source.
        /** This will also remove the source from all the group that the source is a part of.
         * @param     index				The index of the source to remove.
         */
        void sourceRemove(long index);
        
        //! Set the mute state of a source.
        /**
         * @param     index				The index of the source.
         * @param     state				The mute state of the source.
         */
        void sourceSetMute(long index, bool state);
        
        //! Retrieve the existence state of a source.
        /**
         * @param     index				The index of the source.
         * @return						The existence state of the source.
         */
        long sourceGetExistence(long index);
        
        //! Get the radius of a source.
        /**
         * @param     index				The index of the source.
         * @return						The radius of the source.
         */
        double sourceGetRadius(long index);
        
        //! Get the azimuth of a source.
        /**
         * @param     index				The index of the source.
         * @return						The azimuth of the source.
         */
        double sourceGetAzimuth(long index);
        
        //! Get the elevation of a source.
        /**
         * @param     index				The index of the source.
         * @return						The elevation of the source.
         */
        double sourceGetElevation(long index);
        
        //! Get the abscissa of a source.
        /**
         * @param     index				The index of the source.
         * @return						The abscissa of the source.
         */
        double sourceGetAbscissa(long index);
        
        //! Get the ordinate of a source.
        /**
         * @param     index				The index of the source.
         * @return						The ordinate of the source.
         */
        double sourceGetOrdinate(long index);
        
        //! Get the height of a source.
        /**
         * @param     index				The index of the source.
         * @return						The height of the source.
         */
        double sourceGetHeight(long index);
        
        //! Get the rgba color of a source.
        /**
         * @param     index				The index of the source.
         * @return						The rgba color of the source as an array of 4 values (red, green, blue, alpha).
         */
        double* sourceGetColor(long index);
        
        //! Get the text description of a source.
        /**
         * @param     index				The index of the source.
         * @return						The text description.
         */
        std::string sourceGetDescription(long index);
        
        //! Get the number of group a source is owned by.
        /**
         * @param     index				The index of the source.
         * @return		The number of group.
         */
        long sourceGetNumberOfGroups(long index);
        
        //! Get the the group index the source is owned by at a particular index.
        /**
         * @param		sourceIndex			The index of the source.
         * @param		groupIndex			The index of the group
         * @return		The group index.
         */
        long sourceGetGroupIndex(long sourceIndex, long groupIndex);
        
        //! Retrieve the mute state of a source.
        /**
         * @param     index				The index of the source.
         * @return						The mute state of the source.
         */
        long sourceGetMute(long index);
        
        /* ------------------------------------------------------------------------ */
        /* -------------------------------- Groups -------------------------------- */
        /* ------------------------------------------------------------------------ */
        
        //! Add a source to a group.
        /**
         * @param     groupIndex		The index of the group.
         * @param     sourceIndex		The index of the source.
         * @see groupRemoveSource
         */
        void groupSetSource(long groupIndex, long sourceIndex);
        
        //! Remove source from a group.
        /**
         * @param     groupIndex		The index of the group.
         * @param     sourceIndex		The index of the source.
         * @see groupSetSource
         */
        void groupRemoveSource(long groupIndex, long sourceIndex);
        
        //! Set position of a group with polar coordinates.
        /**
         * @param     index				The index of the group.
         * @param     radius			The radius of the group.
         * @param     azimuth			The azimuth of the group.
         * @see groupSetRadius, groupSetAzimuth, groupSetCartesian
         */
        void groupSetPolar(long index, double radius, double azimuth);
        
        //! Set position of a group with polar coordinates.
        /**
         * @param     index				The index of the group.
         * @param     radius			The radius of the group.
         * @param     azimuth			The azimuth of the group.
         * @param     elevation			The elevation of the group.
         * @see groupSetRadius, groupSetAzimuth, groupSetCartesian
         */
        void groupSetPolar(long index, double radius, double azimuth, double elevation);
        
        //! Set radius of a group.
        /**
         * @param     index				The index of the group.
         * @param     radius			The radius of the group.
         * @see groupSetPolar, groupSetAzimuth
         */
        void groupSetRadius(long index, double radius);
        
        //! Set azimuth of a group.
        /**
         * @param     index				The index of the group.
         * @param     azimuth			The azimuth of the group.
         * @see groupSetRadius, groupSetAzimuth
         */
        void groupSetAzimuth(long index, double azimuth);
        
        //! Set elevation of a group.
        /**
         * @param     index				The index of the group.
         * @param     elevation			The elevation of the group.
         * @see groupSetRadius, groupSetAzimuth
         */
        void groupSetElevation(long index, double elevation);
        
        //! Set position of a group with cartesian coordinates.
        /**
         * @param     index				The index of the group.
         * @param     abscissa			The abscissa of the group.
         * @param     ordinate			The ordinate of the group.
         * @see groupSetAbscissa, groupSetOrdinate, groupSetPolar
         */
        void groupSetCartesian(long index, double abscissa, double ordinate);
        
        //! Set position of a group with cartesian coordinates.
        /**
         * @param     index				The index of the group.
         * @param     abscissa			The abscissa of the group.
         * @param     ordinate			The ordinate of the group.
         * @param     height			The height of the group.
         * @see groupSetAbscissa, groupSetOrdinate, groupSetPolar
         */
        void groupSetCartesian(long index, double abscissa, double ordinate, double height);
        
        //! Set abscissa of a group.
        /**
         * @param     index				The index of the group.
         * @param     abscissa			The abscissa of the group.
         * @see groupSetOrdinate
         */
        void groupSetAbscissa(long index, double abscissa);
        
        //! Set ordinate of a group.
        /**
         * @param     index				The index of the group.
         * @param     ordinate			The ordinate of the group.
         * @see groupSetAbscissa
         */
        void groupSetOrdinate(long index, double ordinate);
        
        //! Set height of a group.
        /**
         * @param     index				The index of the group.
         * @param     height			The height of the group.
         * @see groupSetAbscissa
         */
        void groupSetHeight(long index, double height);
        
        //! Set position of a group with relative polar coordinates.
        /**
         * @param     index				The index of the group.
         * @param     radius			The relative radius of the group.
         * @param     azimuth			The relative azimuth of the group.
         * @see groupSetRadius, groupSetAzimuth, groupSetCartesian
         */
        void groupSetRelativePolar(long groupIndex, double radius, double azimuth);
        
        //! Set position of a group with relative polar coordinates.
        /**
         * @param     index				The index of the group.
         * @param     radius			The relative radius of the group.
         * @param     azimuth			The relative azimuth of the group.
         * @param     elevation			The relative elevation of the group.
         * @see groupSetRadius, groupSetAzimuth, groupSetCartesian
         */
        void groupSetRelativePolar(long groupIndex, double radius, double azimuth, double elevation);
        
        //! Set radius of a group with relative value.
        /**
         * @param     index				The index of the group.
         * @param     radius			The relative radius of the group.
         * @see groupSetRadius, groupSetAzimuth, groupSetPolar
         */
        void groupSetRelativeRadius(long groupIndex, double radius);
        
        //! Set azimuth of a group with relative value.
        /**
         * @param     index				The index of the group.
         * @param     azimuth			The relative azimuth of the group.
         * @see groupSetAzimuth, groupSetRadius, groupSetPolar
         */
        void groupSetRelativeAzimuth(long groupIndex, double azimuth);
        
        //! Set elevation of a group with relative value.
        /**
         * @param     index				The index of the group.
         * @param     azimuth			The relative elevation of the group.
         * @see groupSetAzimuth, groupSetRadius, groupSetPolar
         */
        void groupSetRelativeElevation(long groupIndex, double elevation);
        
        //! Set the rgba color of a group.
        /** All values are clipped between 0 and 1.
         * @param     index				The index of the group.
         * @param     red				The red component of the color.
         * @param     green				The green component of the color
         * @param     blue				The blue component of the color
         * @param     alpha				The alpha component of the color
         */
        void groupSetColor(long index, double red, double green, double blue, double alpha);
        
        //! Add a description to a given group.
        /**
         * @param     index				The index of the group.
         * @param     description		The text description of the group.
         */
        void groupSetDescription(long index, std::string description);
        
        //! Remove group.
        /**
         * @param     groupIndex		The index of the group.
         * @see groupRemoveWithSources
         */
        void groupRemove(long groupIndex);
        
        //! Remove group and sources it contains.
        /**
         * @param     groupIndex		The index of the group.
         * @see groupRemove
         */
        void groupRemoveWithSources(long groupIndex);
        
        //! Get the number of sources a group contains.
        /**
         * @param     groupIndex		The index of the group.
         */
        long groupGetNumberOfSources(long groupIndex);
        
        //! Get the the index of a source stored at a particular index by a group.
        /**
         * @param			groupIndex			The index of the group.
         * @param			sourceIndex			The index of the source.
         * @return			The index of the source if it exists, -1 otherwise.
         */
        long groupGetSourceIndex(long groupIndex, long sourceIndex);
        
        
        //! Set the mute state of a group.
        /**
         * @param     index				The index of the group.
         * @param     state				The mute state of the group.
         */
        void groupSetMute(long index, long state);
        
        //! Clean all groups
        void groupClean();
        
        //! Retrieve the existence state of a group.
        /**
         * @param     index				The index of the group.
         * @return						The existence state of the group.
         */
        long groupGetExistence(long index);
        
        //! Get the radius of a group.
        /**
         * @param     index				The index of the group.
         * @return						The radius of the group.
         */
        double groupGetRadius(long index);
        
        //! Get the azimuth of a group.
        /**
         * @param     index				The index of the group.
         * @return						The azimuth of the group.
         */
        double groupGetAzimuth(long index);
        
        //! Get the elevation of a group.
        /**
         * @param     index				The index of the group.
         * @return						The elevation of the group.
         */
        double groupGetElevation(long index);
        
        //! Get the abscissa of a group.
        /**
         * @param     index				The index of the group.
         * @return						The abscissa of the group.
         */
        double groupGetAbscissa(long index);
        
        //! Get the ordinate of a group.
        /**
         * @param     index				The index of the group.
         * @return						The ordinate of the group.
         */
        double groupGetOrdinate(long index);
        
        //! Get the height of a group.
        /**
         * @param     index				The index of the group.
         * @return						The height of the group.
         */
        double groupGetHeight(long index);
        
        //! Get the rgba color of a group.
        /**
         * @param     index				The index of the group.
         * @return						The rgba color of the group as an array of 4 values (red, green, blue, alpha).
         */
        double* groupGetColor(long index);
        
        //! Get the text description of a group.
        /**
         * @param     index				The index of the group.
         * @return						The text description.
         */
        std::string groupGetDescription(long index);
        
        //! Retrieve the mute state of a group.
        /**
         * @param     index				The index of the group.
         * @return						The mute state of the group.
         */
        long groupGetMute(long index);
        
        //! Returns true if a source is is muted in a group.
        /**
         * @param     index				The index of the group.
         * @return						True if a source is is muted in the group, false otherwise.
         */
        bool groupGetIfSourceMuted(long index);
        
        //! Retrieve the next free group index.
        /**
         * @return						The next free group index.
         */
        long groupGetNextIndex();
    };
   
}
#endif