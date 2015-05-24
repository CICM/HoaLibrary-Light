/*
// Copyright (c) 2012-2015 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_SOURCE_LIGHT
#define DEF_HOA_SOURCE_LIGHT

#include "Math.hpp"

namespace hoa
{
    //! Source class is used to simulate punctual sources.
    /** Source class is used to simulate punctual sources and control its.
     */
    class Source
    {
    public:
        class Group;

        typedef  map<ulong, Source*>::const_iterator  const_source_iterator;
        typedef  map<ulong, Source*>::iterator source_iterator;
        typedef  map<ulong, Group*>::const_iterator  const_group_iterator;
        typedef  map<ulong, Group*>::iterator group_iterator;

        //! Manager class is used to control punctual sources and group of sources.
        /** Manager class is used to control punctual sources and group of sources.
         */
        class Manager
        {
        private:
            const double         m_maximum_radius;
            map<ulong, Source*>  m_sources;
            map<ulong, Group*>   m_groups;
            double               m_zoom;

        public:

            //! The manager constructor.
            /**	The manager constructor allocates and initialize the member values.
             *
             * @param     maximumRadius		The maximum radius the sources or groups in the source manager could have
             */
            Manager(const double maximumRadius = 1.) :
                m_maximum_radius(maximumRadius)
            {
                m_zoom = 1;
            }

            //! The manager constructor by copy.
            /**	The manager constructor allocates and initialize the member values.
             *
             * @param     other		It's a contructor by copy an 'other' manager
             */
            Manager(const Manager& other) :
                m_maximum_radius(other.getMaximumRadius())
            {
                m_zoom = other.getZoom();

                for(const_source_iterator it = other.m_sources.begin() ; it != other.m_sources.end() ; it ++)
                {
                    m_sources[it->first] = new Source(*it->second);
                }
                for(const_group_iterator it = other.m_groups.begin() ; it != other.m_groups.end() ; it ++)
                {
                    m_groups[it->first] = new Group(*it->second);

                    map<ulong, Source*>& tmp = it->second->getSources();
                    for (const_source_iterator ti = tmp.begin() ; ti != tmp.end() ; ti ++)
                    {
                        m_groups[it->first]->addSource(m_sources[ti->first]);
                    }
                }
            }

            //! The source manager destructor.
            /** The source manager destructor free the memory.
             */
            ~Manager() noexcept
            {
                clear();
            }

            //! Clear and free the memory
            /** Clear and free the memory
             */
            inline void clear()
            {
                for(const_group_iterator it = m_groups.begin() ; it != m_groups.end() ; ++it)
                {
                    delete it->second;
                }
                for(const_source_iterator it = m_sources.begin() ; it != m_sources.end() ; ++it)
                {
                    delete it->second;
                }
                m_groups.clear();
                m_sources.clear();
            }

            //! Set the zoom factor.
            /** Set the zoom factor between 0 and 1.
             @param     zoom		The zoom factor.
             */
            inline void setZoom(const double zoom)
            {
                m_zoom = Math<double>::clip(zoom, 1. / m_maximum_radius, 1.);
            }

            //! Get the maximum radius of the sources and groups.
            /** Get the maximum radius of the sources and groups.
             @return		The maximum radius.
             */
            inline const double getMaximumRadius() const noexcept
            {
                return m_maximum_radius;
            }

            //! Get the zoom factor value.
            /** Get the zoom factor value.
             @return		The zoom factor value.
             */
            inline const double getZoom() const noexcept
            {
                return m_zoom;
            }

            //! Add a new Source.
            /** Add a new Source.
             @param     index       The index of the new source.
             @param     radius      The radius of the new source.
             @param     azimuth     The azimuth of the new source.
             @param     elevation   The elevation of the new source.
             @return	            The created source.
             */
            inline Source* newSource (const ulong index, const double radius = 0., const double azimuth = 0., const double elevation = 0.) noexcept
            {
                const_source_iterator it = m_sources.find(index);
                if(it == m_sources.end())
                {
                    m_sources[index] = new Source(m_maximum_radius, index, radius, azimuth, elevation);
                    return m_sources[index];
                }
                return it->second;
            }

            //! Remove a Source.
            /** Remove a Source.
             @param     index   The index of the source.
             */
            inline void removeSource (const ulong index) noexcept
            {
                const_source_iterator it = m_sources.find(index);
                if(it != m_sources.end())
                {
                    delete it->second;
                    m_sources.erase(index);
                }
            }

            //! Get the Sources map size.
            /** Get the Sources map size.
             @return     The sources map size.
             */
            inline ulong getSourcesSize() const noexcept
            {
                return m_sources.size();
            }

            //! Check if the Sources map is empty.
            /** Check if the Sources map is empty.
             @return    The state of the source map content.
             */
            inline bool isSourcesEmpty() const noexcept
            {
                return m_sources.empty();
            }

            //! Get the Groups map size.
            /** Get the Groups map size.
             @return     The groups map size.
             */
            inline ulong getGroupsSize() const noexcept
            {
                return m_groups.size();
            }

            //! Check if the Groups map is empty.
            /** Check if the Groups map is empty.
             @return    The state of the groups map content.
             */
            inline bool isGroupsEmpty() const noexcept
            {
                return m_groups.empty();
            }

            //! Add a new Group.
            /** Add a new Group.
             @param     index       The index of the new group.
             @return	            The created group.
             */
            inline Group* newGroup (const ulong index) noexcept
            {
                const_group_iterator it = m_groups.find(index);
                if(it == m_groups.end())
                {
                    m_groups[index] = new Group(this, index);
                    return m_groups[index];
                }
                return it->second;
            }

            //! Remove a Group.
            /** Remove a Group.
             @param     index   The index of the group.
             */
            inline void removeGroup (const ulong index) noexcept
            {
                const_group_iterator it = m_groups.find(index);
                if(it != m_groups.end())
                {
                    map<ulong, Source*>& sources = it->second->getSources();
                    for (source_iterator ti = sources.begin() ; ti != sources.end() ; ++ti)
                    {
                        Source* src = ti->second;
                        if(src)
                        {
                            src->removeGroup(index);
                        }
                    }
                    
                    delete it->second;
                    m_groups.erase(index);
                }
            }

            //! Remove a Group with its sources.
            /** Remove a Group with its sources.
             @param     index   The index of the group.
             */
            void removeGroupWithSources (const ulong index) noexcept
            {
                group_iterator it = m_groups.find(index);
                if(it != m_groups.end())
                {
                    map<ulong, Source*> sources = it->second->getSources();
                    for (source_iterator ti = sources.begin() ; ti != sources.end() ; ++ti)
                    {
                        delete ti->second;
                        m_sources.erase(ti->first);
                    }
                    
                    delete it->second;
                    m_groups.erase(it);
                }
            }

            //! Get one source.
            /** Get one source.
             @param     index   The index of the source.
             @return            A pointer on the source.
             */
            inline Source* getSource(const ulong index)
            {
                const_source_iterator it = m_sources.find(index);
                if(it != m_sources.end())
                {
                    return it->second;
                }
                return NULL;
            }

            //! Get the first const iterator of the Sources map.
            /** Get the first const iterator of the Sources map.
             @return        The first const iterator of the sources map.
             */
            inline const_source_iterator getFirstSource() const noexcept
            {
                return m_sources.begin();
            }

            //! Get the first iterator of the Sources map.
            /** Get the first iterator of the Sources map.
             @return        The const iterator of the sources map.
             */
            inline source_iterator getFirstSource() noexcept
            {
                return m_sources.begin();
            }

            //! Get the last const iterator of the Sources map.
            /** Get the last const iterator of the Sources map.
             @return        The last const iterator of the sources map.
             */
            inline const_source_iterator getLastSource() const noexcept
            {
                return m_sources.end();
            }

            //! Get the last iterator of the Sources map.
            /** Get the last iterator of the Sources map.
             @return        The last iterator of the sources map.
             */
            inline source_iterator getLastSource() noexcept
            {
                return m_sources.end();
            }

            //! Removes the groups which have less than 2 sources
            /** Removes the groups which have less than 2 sources
             */
            inline void cleanEmptyGroup() noexcept
            {
                group_iterator it = m_groups.begin();
                while (it != m_groups.end())
                {
                    if (it->second->m_sources.size() < 2)
                    {
                        group_iterator to = it;
                        ++to;
                        delete it->second;
                        m_groups.erase(it);
                        it  = to;
                    }
                    else
                    {
                        ++it;
                    }
                }
            }

            //! Removes the groups which have exactly the same sources.
            /** Removes the groups which have exactly the same sources.
             */
            inline void cleanDuplicatedGroup() noexcept
            {
                for(group_iterator it = m_groups.begin() ; it != m_groups.end(); ++it)
                {
                    group_iterator ti = it;
                    while(ti != m_groups.end())
                    {
                        if(it->first != ti->first && *it->second == *ti->second)
                        {
                            group_iterator to = ti;
                            ++to;
                            delete ti->second;
                            m_groups.erase(ti);
                            ti = to;
                        }
                        else
                        {
                            ++ti;
                        }
                    }
                }
            }

            //! Get one group.
            /** Get one group.
             @param     index   The index of the group.
             @return            A pointer on the group.
             */
            inline Group* getGroup(const ulong index)
            {
                const_group_iterator it = m_groups.find(index);
                if(it != m_groups.end())
                {
                    return it->second;
                }
                return NULL;
            }

            //! Get the first const iterator of the Groups map.
            /** Get the first const iterator of the Groups map.
             @return        The first const iterator of the groups map.
             */
            inline const_group_iterator getFirstGroup() const noexcept
            {
                return m_groups.begin();
            }

            //! Get the first iterator of the Groups map.
            /** Get the first iterator of the Groups map.
             @return        The first iterator of the groups map.
             */
            inline group_iterator getFirstGroup() noexcept
            {
                return m_groups.begin();
            }

            //! Get the last const iterator of the Groups map.
            /** Get the last const iterator of the Groups map.
             @return        The last const iterator of the groups map.
             */
            inline const_group_iterator getLastGroup() const noexcept
            {
                return m_groups.end();
            }

            //! Get the last iterator of the Groups map.
            /** Get the last iterator of the Groups map.
             @return        The last iterator of the groups map.
             */
            inline group_iterator getLastGroup() noexcept
            {
                return m_groups.end();
            }
        };

		//! Set the position of the source with polar coordinates.
		/** Set the position of the source with polar coordinates.
			@param     radius			The radius of the source.
			@param     azimuth			The azimuth of the source.
            @param     elevation        The elevation of the source.
			@see setCoordinatesCartesian
         */
		inline void setCoordinatesPolar(const double radius, const double azimuth, const double elevation = 0.)
		{
			setRadius(radius);
		    setAzimuth(azimuth);
		    setElevation(elevation);
       	}

		//! Set the radius of the source.
		/** Set the radius of the source.
			@param     radius			The radius of the source.
			@see getRadius
         */
		inline void setRadius(const double radius)
		{
			if(m_maximum_radius >= 0)
		    {
		        if(radius < -m_maximum_radius || radius > m_maximum_radius)
		            return;
		    }
		    m_radius = max(radius, (double)0.);
		    notifyCoordinates();
		}

		//! Set the azimuth of the source.
		/** Set the azimuth of the source.
			@param     azimuth			The azimuth of the source.
         */
		inline void setAzimuth(const double azimuth)
		{
			m_azimuth = Math<double>::wrap_twopi(azimuth);
			notifyCoordinates();
		}

        //! Set the elevation of the source.
		/** Set the elevation of the source.
         @param     elevation			The elevation of the source.
         */
		inline void setElevation(const double elevation)
		{
			m_elevation = Math<double>::wrap_pi(elevation);
		    if(m_elevation > HOA_PI2 || m_elevation < -HOA_PI2)
		    {
		        m_azimuth = Math<double>::wrap_twopi(m_azimuth + HOA_PI);
		        m_elevation = -elevation;
		    }
		    notifyCoordinates();
		}

        //! Set the position of the source with cartesian coordinates.
		/** Set the position of the source with cartesian coordinates.
            @param     abscissa		The abscissa of the source.
            @param     ordinate		The ordinate of the source.
            @param     height		The height of the source.
         */
		inline void setCoordinatesCartesian(const double abscissa, const double ordinate, const double height = 0.)
		{
			setRadius(Math<double>::radius(abscissa, ordinate, height));
        	setAzimuth(Math<double>::azimuth(abscissa, ordinate, height));
        	setElevation(Math<double>::elevation(abscissa, ordinate, height));
        }

		//! Set the abscissa of the source.
		/** Set the abscissa of the source.
			@param     abscissa		The abscissa of the source.
         */
		inline void setAbscissa(const double abscissa)
		{
            setCoordinatesCartesian(abscissa, getOrdinate(), getHeight());
		}

		//! Set the ordinate of the source.
		/** Set the ordinate of the source.
			@param     ordinate		The ordinate of the source.
         */
		inline void setOrdinate(const double ordinate)
		{
            setCoordinatesCartesian(getAbscissa(), ordinate, getHeight());
        }

        //! Set the height of the source.
		/** Set the height of the source.
         @param     height		The height of the source.
         */
		inline void setHeight(const double height)
		{
            setCoordinatesCartesian(getAbscissa(), getOrdinate(), height);
		}

		//! Set the color of the source.
		/** Set the color of the source.
		 * @param     red		The red component of the color.
		 * @param     green		The green component of the color.
		 * @param     blue		The blue component of the color.
		 * @param     alpha		The alpha component of the color.
         */
		inline void setColor(const double red, const double green, const double blue, const double alpha)
		{
			m_color[0]	=  Math<double>::clip(red, 0., 1.);
		    m_color[1]	=  Math<double>::clip(green, 0., 1.);
		    m_color[2]	=  Math<double>::clip(blue, 0., 1.);
		    m_color[3]	=  Math<double>::clip(alpha, 0., 1.);
        }

		//! Set the description of the source.
		/** Set the description of the source.
		 * @param     description		The text description of the source.
         */
		inline void setDescription(const string description)
		{
			m_description = description;
		}

		//! Set the mute state of the source.
		/** Set the mute state of the source.
			@param     state		The mute state of the source.
			@see getMute
         */
		inline void setMute(const bool state)
		{
			m_mute = state;
			notifyMute();
		}

        //! Get the maximum radius of the source.
		/** Get the maximum radius of the source.
			@return		The maximum radius of the source.
         */
        inline const double getMaximumRadius() const noexcept
        {
            return m_maximum_radius;
        }

        //! Get the index of the source.
		/** Get the index of the source.
			@return		The index of the source.
         */
        inline const ulong getIndex() const noexcept
        {
            return m_index;
        }

		//! Get the radius of the source.
		/** Get the radius of the source.
			@return		The radius of the source.
			@see setRadius, setCoordinatesPolar
         */
		inline const double	getRadius()	const noexcept
		{
			return m_radius;
		}

		//! Get the azimuth of the source.
		/** Get the azimuth of the source.
			@return		The azimuth of the source.
			@see setAzimuth, setCoordinatesPolar
         */
		inline const double	getAzimuth() const noexcept
		{
			return m_azimuth;
		}

        //! Get the elevation of the source.
		/** Get the elevation of the source.
            @return		The elevation of the source.
            @see setElevation, setCoordinatesPolar
         */
		inline const double	getElevation() const noexcept
		{
			return m_elevation;
		}

		//! Get the abscissa of the source.
		/** Get the abscissa of the source.
			@return		The abscissa of the source.
			@see setAbscissa, setCoordinatesCartesian
         */
		inline const double	getAbscissa() const
		{
			return Math<double>::abscissa(m_radius, m_azimuth, m_elevation);
		}

		//! Get the ordinate of the source.
		/** Get the ordinate of the source.
			@return		The ordinate of the source.
			@see setOrdinate, setCoordinatesCartesian
         */
		inline const double	getOrdinate() const
		{
			return Math<double>::ordinate(m_radius, m_azimuth, m_elevation);
		}

        //! Get the height of the source.
		/** Get the height of the source.
         @return		The height of the source.
         @see setHeight, setCoordinatesCartesian
         */
		inline const double	getHeight() const
		{
			return Math<double>::height(m_radius, m_azimuth, m_elevation);
		}

		//! Get the color of the source.
		/** Get the color of the source.
			@return		The rgba color of the source as an array of 4 double.
			@see setColor
         */
		inline const double* getColor() const noexcept
		{
			return m_color;
		}

		//! Get the description of the source.
		/** Get the description of the source.
			@return		The description of the source.
			@see setDescription
         */
		inline const string	getDescription() const noexcept
		{
			return m_description;
		}

		//! Get the mute state of the source.
		/** Get the mute state of the source.
			@return		The mute state of the source.
			@see setMute
         */
		inline const bool getMute() const noexcept
		{
			return m_mute;
		}

        //! Get the size of the Groups map.
        /** Get the size of the Groups map.
         @return    The group map size
         */
        inline ulong getGroupsSize() const noexcept
        {
            return m_groups.size();
        }

        //! Check if the Groups map is empty.
        /** Check if the Groups map is empty.
         @return    The state of the groups map content.
         */
        inline bool isGroupsEmpty() const noexcept
        {
            return m_groups.empty();
        }

        //! Get the Groups map.
        /** Get the Groups map.
         @return    A reference of the groups map.
         */
        inline map<ulong, Group*>& getGroups() noexcept
        {
            return m_groups;
        }

        //! Group class is used to control punctual sources.
        /** Group class is used to control punctual sources.
         */
        class Group
        {
        friend class Source;
        friend class Manager;

        private:
            ulong                   m_index;
            map<ulong, Source*>     m_sources;
            string                  m_description;
            double			        m_color[4];
            double			        m_centroid_x;
            double			        m_centroid_y;
            double			        m_centroid_z;
            double                  m_maximum_radius;
            bool                    m_mute;
            bool                    m_subMute;
            const Manager*          m_manager;

            //! The group constructor.
            /**	The group constructor allocates and initialize the member values for a source group.
             @param     manager		A pointer on a manager object
             @param     index       The index of the group
             */
            Group(const Manager* manager, const ulong index) : m_manager(manager)
            {
                m_maximum_radius = m_manager->getMaximumRadius();
                m_index = index;
                setColor(0.2, 0.2, 0.2, 1.);
                m_description = "";
                computeCentroid();
                m_mute = false;
            }

            //! The group constructor by copy.
            /**	The group constructor allocates and initialize the member values for a source group.
             @param     other		It's a contructor by copy an 'other' group
             */
            Group(const Group& other) : m_manager(other.m_manager)
            {
                m_maximum_radius = m_manager->getMaximumRadius();
                m_index = other.getIndex();
                const double* color = other.getColor();
                setColor(color[0], color[1], color[2], color[3]);
                m_description = other.getDescription();
                m_mute = other.getMute();
                m_sources = other.m_sources;
                computeCentroid();
            }

            //! The source group destructor.
            /**	The source group destructor free the memory.
             */
            ~Group() noexcept
            {
                for (const_source_iterator it = m_sources.begin() ; it != m_sources.end() ; it ++)
                {
                    m_sources[it->first]->removeGroup(m_index);
                }
                m_sources.clear();
            }

            //! Compute the group position for each moving of its sources.
            /** Compute the group position for each moving of its sources.
             */
            inline void notifyCoordinates() noexcept
            {
                computeCentroid();
            }

            //! Check the group mute state for each change of mute state of its sources.
            /** Check the group mute state for each change of mute state of its sources.
             */
            inline void notifyMute() noexcept
            {
                ulong numberOfMutedSources = 0;
                for (const_source_iterator it = m_sources.begin() ; it != m_sources.end() ; it ++)
                {
                    if (it->second->getMute())
                        numberOfMutedSources ++;
                }
                if (numberOfMutedSources)
                    m_subMute = true;
                else
                    m_subMute = false;
                if (numberOfMutedSources == m_sources.size())
                    m_mute = true;
                else
                    m_mute = false;
            }

            //! Compute the group position.
            /** Compute the group position.
             */
            inline void computeCentroid()
            {
                m_centroid_x = 0.;
                m_centroid_y = 0.;
                m_centroid_z = 0.;
                if(m_sources.size())
                {
                    for (const_source_iterator it = m_sources.begin() ; it != m_sources.end() ; it ++)
                    {
                        m_centroid_x += it->second->getAbscissa();
                        m_centroid_y += it->second->getOrdinate();
                        m_centroid_z += it->second->getHeight();
                    }
                    m_centroid_x /= m_sources.size();
                    m_centroid_y /= m_sources.size();
                    m_centroid_z /= m_sources.size();
                }
            }

            //! Compute the new polar coordinates of the Group.
            /** Compute the new polar coordinates of the Group.
             @param     radius      The radius factor of shifting.
             @param     azimuth     The azimuth factor of shifting.
             @param     elevation   The elevation factor of shifting.
             */
            inline void shiftPolar(const double radius, const double azimuth, const double elevation = 0.)
            {
                shiftRadius(radius);
                shiftAzimuth(azimuth);
                shiftElevation(elevation);
            }

            //! Compute the new radius of the Group.
            /** Compute the new radius of the Group.
             @param     radius      The radius factor of shifting.
             */
            void shiftRadius(double radius)
            {
                for (const_source_iterator it = m_sources.begin() ; it != m_sources.end() ; it ++)
                {
                    it->second->setRadius(radius + it->second->getRadius());
                }
            }

            //! Compute the new azimuth of the Group.
            /** Compute the new azimuth of the Group.
             @param     azimuth     The azimuth factor of shifting.
             */
            inline void shiftAzimuth(double azimuth)
            {
                for (const_source_iterator it = m_sources.begin() ; it != m_sources.end() ; it ++)
                {
                    it->second->setAzimuth(azimuth + it->second->getAzimuth());
                }
            }

            //! Compute the elevation of the Group.
            /** Compute the elevation of the Group.
             @param     elevation   The elevation factor of shifting.
             */
            inline void shiftElevation(double elevation)
            {
                for (const_source_iterator it = m_sources.begin() ; it != m_sources.end() ; it ++)
                {
                    it->second->setElevation(elevation + it->second->getElevation());
                }
            }

            //! Compute the new cartesian coordinates of the Group.
            /** Compute the new cartesian coordinates of the Group.
             @param     abscissa    The abscissa factor of shifting.
             @param     ordinate    The ordinate factor of shifting.
             @param     height      The height factor of shifting.
             */
            inline void shiftCartesian(const double abscissa, const double ordinate, const double height = 0.)
            {
                shiftAbscissa(abscissa);
                shiftOrdinate(ordinate);
                shiftHeight(height);
            }

            //! Compute the new abscissa of the Group.
            /** Compute the new abscissa of the Group.
             @param     abscissa    The abscissa factor of shifting.
             */
            void shiftAbscissa(double abscissa)
            {
                if(m_maximum_radius >= 0)
                {
                    if(abscissa < 0.)
                    {
                        double refValue = -m_maximum_radius * 2.;
                        for (const_source_iterator it = m_sources.begin() ; it != m_sources.end() ; it ++)
                        {
                            double circleValue = -sqrt(m_maximum_radius * m_maximum_radius - it->second->getOrdinate() * it->second->getOrdinate());
                            if(circleValue - it->second->getAbscissa() > refValue)
                                refValue = circleValue - it->second->getAbscissa();
                        }
                        if(abscissa < refValue)
                        {
                            abscissa = refValue;
                        }
                    }
                    else if(abscissa >= 0.)
                    {
                        double refValue = m_maximum_radius * 2.;
                        for (const_source_iterator it = m_sources.begin() ; it != m_sources.end() ; it ++)
                        {
                            double circleValue = sqrt(m_maximum_radius * m_maximum_radius - it->second->getOrdinate() * it->second->getOrdinate());
                            if(circleValue - it->second->getAbscissa() < refValue)
                                refValue = circleValue - it->second->getAbscissa();
                        }
                    }
                }
                for (const_source_iterator it = m_sources.begin() ; it != m_sources.end() ; it ++)
                {
                    it->second->setAbscissa(abscissa + it->second->getAbscissa());
                }
            }

            //! Compute the new ordinate of the Group.
            /** Compute the new ordinate of the Group.
             @param     ordinate    The ordinate factor of shifting.
             */
            void shiftOrdinate(double ordinate)
            {
                if(m_maximum_radius >= 0)
                {
                    if(ordinate < 0.)
                    {
                        double refValue = -m_maximum_radius * 2.;
                        for (const_source_iterator it = m_sources.begin() ; it != m_sources.end() ; it ++)
                        {
                            double circleValue = -sqrt(m_maximum_radius * m_maximum_radius - it->second->getAbscissa() * it->second->getAbscissa());
                            if(circleValue - it->second->getOrdinate() > refValue)
                                refValue = circleValue - it->second->getOrdinate();
                        }
                        if(ordinate < refValue)
                        {
                            ordinate = refValue;
                        }
                    }
                    else if(ordinate >= 0.)
                    {
                        double refValue = m_maximum_radius * 2.;
                        for (const_source_iterator it = m_sources.begin() ; it != m_sources.end() ; it ++)
                        {
                            double circleValue = sqrt(m_maximum_radius * m_maximum_radius - it->second->getAbscissa() * it->second->getAbscissa());
                            if(circleValue - it->second->getOrdinate() < refValue)
                                refValue = circleValue - it->second->getOrdinate();
                        }
                    }
                }
                for (const_source_iterator it = m_sources.begin() ; it != m_sources.end() ; it ++)
                {
                    it->second->setOrdinate(ordinate + it->second->getOrdinate());
                }
            }

            //! Compute the new height of the Group.
            /** Compute the new height of the Group.
             @param     height      The height factor of shifting.
             */
            void shiftHeight(double height)
            {
                if(m_maximum_radius >= 0)
                {
                    if(height < 0.)
                    {
                        double refValue = -m_maximum_radius * 2.;
                        for (const_source_iterator it = m_sources.begin() ; it != m_sources.end() ; it ++)
                        {
                            double circleValue = -sqrt(m_maximum_radius * m_maximum_radius - it->second->getAbscissa() * it->second->getAbscissa());
                            if(circleValue - it->second->getHeight() > refValue)
                                refValue = circleValue - it->second->getHeight();
                        }
                        if(height < refValue)
                        {
                            height = refValue;
                        }
                    }
                    else if(height >= 0.)
                    {
                        double refValue = m_maximum_radius * 2.;
                        for (const_source_iterator it = m_sources.begin() ; it != m_sources.end() ; it ++)
                        {
                            double circleValue = sqrt(m_maximum_radius * m_maximum_radius - it->second->getAbscissa() * it->second->getAbscissa());
                            if(circleValue - it->second->getHeight() < refValue)
                                refValue = circleValue - it->second->getHeight();
                        }
                    }
                }
                for (const_source_iterator it = m_sources.begin() ; it != m_sources.end() ; it ++)
                {
                    it->second->setHeight(height + it->second->getHeight());
                }
            }

        public:

            //! Add a new Source.
            /** Add a new Source (and add the group to the source).
             @param     source  The source to add.
             @return            The state of the addition of the source.
             */
            inline bool addSource(Source* source) noexcept
            {
                if(source)
                {
                    const_source_iterator it = m_sources.find(source->getIndex());
                    if(it == m_sources.end())
                    {
                        source->addGroup(this);
                        m_sources[source->getIndex()] = source;
                        computeCentroid();
                        notifyMute();
                        return true;
                    }
                }
                return false;
            }

            //! Remove a Source.
            /** Remove a Source.
             @param     index   The index of the source.
             */
            inline void removeSource(const ulong index) noexcept
            {
                m_sources.erase(index);
                computeCentroid();
            }

            //! Set the position of the group with polar coordinates.
            /** Set the position of the group with polar coordinates.
             @param     radius			The radius of the group.
             @param     azimuth			The azimuth of the group.
             @param     elevation			The elevation of the group.
             @see setRelativeCoordinatesPolar, setCoordinatesCartesian
             */
            inline void setCoordinatesPolar(const double radius, const double azimuth, const double elevation = 0.)
            {
                setCoordinatesCartesian(Math<double>::abscissa(radius, azimuth, elevation), Math<double>::ordinate(radius, azimuth, elevation));
            }

            //! Set the radius of the group.
            /** Set the radius of the group.
             @param     radius			The radius of the group.
             @see getRadius
             */
            inline void setRadius(const double radius)
            {
                setCoordinatesCartesian(Math<double>::abscissa(radius, getAzimuth(), getElevation()), Math<double>::ordinate(radius, getAzimuth(), getElevation()), Math<double>::height(radius, getAzimuth(), getElevation()));
            }

            //! Set the azimuth of the group.
            /** Set the azimuth of the group.
             @param     azimuth			The azimuth of the group.
             */
            inline void setAzimuth(const double azimuth)
            {
                setCoordinatesCartesian(Math<double>::abscissa(getRadius(), azimuth, getElevation()), Math<double>::ordinate(getRadius(), azimuth, getElevation()), Math<double>::height(getRadius(), azimuth, getElevation()));
            }

            //! Set the elevation of the group.
            /** Set the elevation of the group.
             @param     elevation			The elevation of the group.
             */
            inline void setElevation(const double elevation)
            {
                setCoordinatesCartesian(Math<double>::abscissa(getRadius(), getAzimuth(), elevation), Math<double>::ordinate(getRadius(), getAzimuth(), elevation), Math<double>::height(getRadius(), getAzimuth(), elevation));
            }

            //! Set the position of the group with cartesian coordinates.
            /** Set the position of the group with cartesian coordinates.
             @param     abscissa		The abscissa of the group.
             @param     ordinate		The ordinate of the group.
             @param     height		The height of the group.
             */
            inline void setCoordinatesCartesian(double abscissa, double ordinate, double height = 0.)
            {
                abscissa = abscissa - getAbscissa();
                ordinate = ordinate - getOrdinate();
                height = height - getHeight();
                shiftCartesian(abscissa, ordinate, height);
                computeCentroid();
            }

            //! Set the abscissa of the group.
            /** Set the abscissa of the group.
             @param     abscissa		The abscissa of the group.
             */
            inline void setAbscissa(const double abscissa)
            {
                double aAbscissaOffset = abscissa - getAbscissa();
                shiftAbscissa(aAbscissaOffset);
                computeCentroid();
            }

            //! Set the ordinate of the group.
            /** Set the ordinate of the group.
             @param     ordinate		The ordinate of the group.
             */
            inline void setOrdinate(const double ordinate)
            {
                double aOrdinateOffset = ordinate - getOrdinate();
                shiftOrdinate(aOrdinateOffset);
                computeCentroid();
            }

            //! Set the height of the group.
            /** Set the height of the group.
             @param     height		The height of the group.
             */
            inline void setHeight(const double height)
            {
                double aHeightOffset = height - getHeight();
                shiftHeight(aHeightOffset);
                computeCentroid();
            }

            //! Set the color of the group.
            /** Set the color of the group.
             @param     red			The red component of the color.
             @param     green		The green component of the color.
             @param     blue			The blue component of the color.
             @param     alpha		The alpha component of the color.
             */
            inline void setColor(const double red, const double green, const double blue, const double alpha)
            {
                m_color[0]	=  Math<double>::clip(red, 0., 1.);
                m_color[1]	=  Math<double>::clip(green, 0., 1.);
                m_color[2]	=  Math<double>::clip(blue, 0., 1.);
                m_color[3]	=  Math<double>::clip(alpha, 0., 1.);
            }

            //! Set the description of the group.
            /** Set the description of the group.
             @param     description		The text description of the group.
             */
            inline void setDescription(const string description)
            {
                m_description = description;
            }

            //! Set the mute state of the group.
            /** Set the mute state of the group.
             @param     state		The mute state of the group.
             @see getMute
             */
            inline void setMute(const bool state)
            {
                for (const_source_iterator it = m_sources.begin() ; it != m_sources.end() ; it ++)
                {
                    it->second->setMute(state);
                }
                m_mute = state;
            }

            //! Set the position of the group with relative polar coordinates.
            /** Set the position of the group with relative polar coordinates.
             @param     radius			The relative radius value.
             @param     azimuth			The relative azimuth value.
             @param     elevation		The relative elevation value.
             @see setCoordinatesPolar, setRelativeRadius
             */
            inline void setRelativeCoordinatesPolar(const double radius, const double azimuth, const double elevation = -HOA_PI2)
            {
                setRelativeRadius(radius);
                setRelativeAzimuth(azimuth);
                setRelativeElevation(elevation);
            }

            //! Set the radius of the group with a relative value.
            /** Set the radius of the group with a relative value.
             @param     radius			The relative radius value.
             @see setCoordinatesPolar, setRelativeRadius
             */

            inline void setRelativeRadius(const double radius)
            {
                double aRadiusOffset = max(radius, (double)0.) - getRadius();
                shiftRadius(aRadiusOffset);
                computeCentroid();
            }

            //! Set the azimuth of the group with a relative value.
            /** Set the azimuth of the group with a relative value.
             @param     azimuth			The relative azimuth value.
             @see setCoordinatesPolar, setRelativeRadius
             */
            inline void setRelativeAzimuth(double azimuth)
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

            //! Set the elevation of the group with a relative value.
            /** Set the elevation of the group with a relative value.
             @param     elevation			The relative elevation value.
             @see setCoordinatesPolar, setRelativeRadius
             */
            inline void setRelativeElevation(double elevation)
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

            //! Get the manager of the Group.
            /** Get the manager of the Group.
                @return     A pointer on the manager of the group.
             */
            inline const Manager* getManager() const noexcept
            {
                return m_manager;
            }

            //! Get the maximum radius of the Group.
            /** Get the maximum radius of the Group.
             @return    The maximum radius of the group.
             */
            inline const double getMaximumRadius() const noexcept
            {
                return m_maximum_radius;
            }

            //! Get the index of the Group.
            /** Get the index of the Group.
             @return    The index of the group.
             */
            inline const ulong getIndex() const noexcept
            {
                return m_index;
            }

            //! Get the radius of the group.
            /** Get the radius of the group.
             @return		The radius of the group.
             @see setRadius, setCoordinatesPolar
             */
            inline const double	getRadius()	const noexcept
            {
                return Math<double>::radius(m_centroid_x, m_centroid_y, m_centroid_z);
            }

            //! Get the azimuth of the group.
            /** Get the azimuth of the group.
             @return		The azimuth of the group.
             @see setAzimuth, setCoordinatesPolar
             */
            inline const double	getAzimuth() const noexcept
            {
                return Math<double>::azimuth(m_centroid_x, m_centroid_y, m_centroid_z) + HOA_PI2;
            }

            //! Get the elevation of the group.
            /** Get the elevation of the group.
             @return		The elevation of the group.
             @see setAzimuth, setCoordinatesPolar
             */
            inline const double	getElevation() const noexcept
            {
                return Math<double>::elevation(m_centroid_x, m_centroid_y, m_centroid_z);
            }

            //! Get the abscissa of the group.
            /** Get the abscissa of the group.
             @return		The abscissa of the group.
             @see setAbscissa, setCoordinatesCartesian
             */
            inline const double	getAbscissa() const noexcept
            {
                return m_centroid_x;
            }

            //! Get the ordinate of the group.
            /** Get the ordinate of the group.
             @return		The ordinate of the group.
             @see setOrdinate, setCoordinatesCartesian
             */
            inline const double	getOrdinate() const noexcept
            {
                return m_centroid_y;
            }

            //! Get the height of the group.
            /** Get the height of the group.
             @return		The height of the group.
             @see setOrdinate, setCoordinatesCartesian
             */
            inline const double	getHeight() const noexcept
            {
                return m_centroid_z;
            }

            //! Get the color of the group.
            /** Get the color of the group.
             @return		The rgba color of the group as an array of 4 double numbers.
             @see setColor
             */
            inline const double*	getColor()	const noexcept
            {
                return m_color;
            }

            //! Get the description of the group.
            /** Get the description of the group.
             @return		The description of the group.
             @see setDescription
             */
            inline const string	getDescription() const noexcept
            {
                return m_description;
            }

            //! Get the mute state of the group.
            /** Get the mute state of the group.
             @return		The mute state of the group.
             @see setMute
             */
            inline const bool getMute()	const noexcept
            {
                return m_mute;
            }

            //! Get the general mute state of sources of the group.
            /** Get the general mute state of sources of the group.
             @return        The sub general mute state of sources of the group.
             */
            inline const bool getSubMute()	const noexcept
            {
                return m_subMute;
            }

            //! Get the size of the Sources map.
            /** Get the size of the Sources map.
             @return    The sources map size.
             */
            inline ulong getSourcesSize() const noexcept
            {
                return m_sources.size();
            }

            //! Check if the Sources map is empty.
            /** Check if the Sources map is empty.
             @return    The state of the sources map content.
             */
            inline bool isSourcesEmpty() const noexcept
            {
                return m_sources.empty();
            }

            //! Get the Sources map.
            /** Get the Sources map.
             @return    A reference of the sources map.
             */
            inline map<ulong, Source*>& getSources() noexcept
            {
                return m_sources;
            }

            inline bool operator== (Group& other)
            {
                if (m_sources.size() == other.m_sources.size())
                {
                    ulong clones = 0;
                    for (const_source_iterator it = m_sources.begin() ; it != m_sources.end() ; it ++)
                    {
                        const_source_iterator ti = other.m_sources.find(it->first);
                        if (ti != other.m_sources.end())
                            clones ++;
                    }
                    if (clones == m_sources.size())
                        return true;
                }
                return false;
            }
        };

    private:
        ulong                m_index;
        double		         m_radius;
		double		         m_azimuth;
        double               m_elevation;
		double		         m_color[4];
		string               m_description;
		map<ulong, Group*>   m_groups;
		double               m_maximum_radius;
		bool                 m_mute;

		//! The source constructor.
		/**	The source constructor allocates and initialize the member values for a source.
            @param     maximumRadius    The maximum radius of the source.
            @param     index            The index of the source.
			@param     radius			The radius of the source.
			@param     azimuth			The azimuth of the source.
            @param     elevation		The elevation of the source.
		 */
		Source(const double maximumRadius, const ulong index, const double radius = 0., const double azimuth = 0., const double elevation = 0.)
		{
            m_maximum_radius = maximumRadius;
            m_index = index;
            m_radius = radius;
            m_azimuth = azimuth;
            m_elevation = elevation;
            setColor(0.2, 0.2, 0.2, 1.);
            m_description = "";
            m_mute = false;
      	}

        //! The source constructor by copy.
		/**	The source constructor allocates and initialize the member values for a source.
            @param     othe    It's a constructor by copy an 'other' source.
		 */
      	Source(const Source& other)
		{
            m_maximum_radius = other.getIndex();
            m_index = other.getIndex();
			m_radius = other.getRadius();
			m_azimuth = other.getAzimuth();
			m_elevation = other.getElevation();
			const double* color = other.getColor();
			setColor(color[0], color[1], color[2], color[3]);
			m_description = other.getDescription();
			m_mute = other.getMute();
            m_groups = other.m_groups;
      	}

      	//! The source destructor.
        /**	The source destructor free the memory.
         */
		~Source() noexcept
		{
        	for (const_group_iterator it = m_groups.begin() ; it != m_groups.end() ; it ++)
		    {
                m_groups[it->first]->removeSource(m_index);
            }
        	m_groups.clear();
        }

		//! Add a new group.
		/** Add a new group.
		 @param     group   The group to add.
		 @return            The state of the addition of the group.
		 */
        inline bool addGroup(Group* group) noexcept
        {
            const_group_iterator it = m_groups.find(group->getIndex());
            if(it == m_groups.end() && group)
            {
                m_groups[group->getIndex()] = group;
                return true;
            }
            return false;
        }

        //! Remove a group.
        /** Remove a group.
         @param     index   The index of the group.
         */
        inline void removeGroup(const ulong index) noexcept
        {
            m_groups.erase(index);
        }

        //! Call the groups of the source for each moving to compute their new position.
        /** Call the groups of the source for each moving to compute their new position.
         */
        inline void notifyCoordinates() noexcept
        {
            for (const_group_iterator it = m_groups.begin() ; it != m_groups.end() ; it ++)
		    {
                it->second->notifyCoordinates();
            }
        }

        //! Call the groups of the source for each change of its mute state to check their mute state.
        /** Call the groups of the source for each change of its mute state to check their mute state.
         */
        inline void notifyMute() noexcept
        {
            for (const_group_iterator it = m_groups.begin() ; it != m_groups.end() ; it ++)
		    {
                it->second->notifyMute();
            }
        }
    };

}
#endif
