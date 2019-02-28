/*
// Copyright (c) 2012-2017 CICM - Universite Paris 8 - Labex Arts H2H.
// Authors :
// 2012: Pierre Guillot, Eliott Paris & Julien Colafrancesco.
// 2012-2015: Pierre Guillot & Eliott Paris.
// 2015: Pierre Guillot & Eliott Paris & Thomas Le Meur (Light version)
// 2016-2017: Pierre Guillot.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#pragma once

#include "Hoa_Math.hpp"

//! @cond

namespace hoa
{
    //! The source class is used to simulate punctual sources.
    /** The source class is used to simulate punctual sources and control its.
     */
    class Source
    {
    public:
        class Group;

        typedef  std::map<size_t, Source*>::iterator          source_iterator;
        typedef  std::map<size_t, Source*>::const_iterator    const_source_iterator;
        typedef  std::map<size_t, Group*>::iterator           group_iterator;
        typedef  std::map<size_t, Group*>::const_iterator     const_group_iterator;

        //! The manager class is used to control punctual sources and group of sources.
        /** The manager class is used to control punctual sources and group of sources.
         */
        class Manager
        {
        private:
            const double        m_maximum_radius;
            std::map<size_t, Source*> m_sources;
            std::map<size_t, Group*>  m_groups;
            double              m_zoom;

        public:

            //! The manager constructor.
            /**	The manager constructor allocates and initialize the member values.
             *
             * @param     maximumRadius		The maximum radius the sources or groups in the source manager could have
             */
            Manager(const double maximumRadius = 1.) : m_maximum_radius(maximumRadius), m_zoom(1)
            {
                ;
            }

            //! The manager constructor by copy.
            /**	The manager constructor allocates and initialize the member values.
             *
             * @param     other		It's a contructor by copy an 'other' manager
             */
            Manager(const Manager& other) : m_maximum_radius(other.m_maximum_radius), m_zoom(other.m_zoom)
            {
                for(const_source_iterator it = other.m_sources.begin() ; it != other.m_sources.end() ; it ++)
                {
                    m_sources[it->first] = new Source(*it->second);
                }
                for(const_group_iterator it = other.m_groups.begin() ; it != other.m_groups.end() ; it ++)
                {
                    Group* grp = new Group(*it->second);
                    m_groups[it->first] = grp;

                    std::map<size_t, Source*>& tmp = it->second->getSources();
                    for (source_iterator ti = tmp.begin() ; ti != tmp.end() ; ti ++)
                    {
                        grp->addSource(m_sources[ti->first]);
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
                for(group_iterator it = m_groups.begin() ; it != m_groups.end() ; ++it)
                {
                    delete it->second;
                }
                for(source_iterator it = m_sources.begin() ; it != m_sources.end() ; ++it)
                {
                    delete it->second;
                }
                m_groups.clear();
                m_sources.clear();
            }

            //! Removes all groups.
            /** Removes all groups.
             */
            inline void clearGroups()
            {
                for(group_iterator it = m_groups.begin() ; it != m_groups.end() ; ++it)
                {
                    delete it->second;
                }
                m_groups.clear();
            }

            //! Set the zoom factor.
            /** Set the zoom factor between 0 and 1.
             @param     zoom		The zoom factor.
             */
            inline void setZoom(const double zoom)
            {
                m_zoom = clip(zoom, 1. / m_maximum_radius, 1.);
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

            //! Add a new Source to the manager.
            /** Add a new Source to the map container of the manager.
             @param     index       The index of the new source.
             @param     radius      The radius of the new source.
             @param     azimuth     The azimuth of the new source.
             @param     elevation   The elevation of the new source.
             @return	            The created source.
             */
            inline Source* newSource (const size_t index, const double radius = 0., const double azimuth = 0., const double elevation = 0.) noexcept
            {
                source_iterator it = m_sources.find(index);
                if(it == m_sources.end())
                {
                    Source* src = new Source(m_maximum_radius, index, radius, azimuth, elevation);
                    m_sources[index] = src;
                    return src;
                }
                return it->second;
            }

            //! Remove a Source from the manager.
            /** Remove a Source from the map container of the manager.
             @param     index   The index of the source.
             */
            inline void removeSource (const size_t index) noexcept
            {
                source_iterator it = m_sources.find(index);
                if(it != m_sources.end())
                {
                    delete it->second;
                    m_sources.erase(index);
                    cleanDuplicatedGroup();
                    cleanEmptyGroup();
                }
            }

            //! Get the Sources map size of the manager.
            /** Get the Sources map size of the manager.
             @return     The sources map size.
             */
            inline size_t getNumberOfSources() const noexcept
            {
                return m_sources.size();
            }

            //! Check if the Sources map of the manager is empty.
            /** Check if the Sources map of the manager is empty.
             @return    The state of the source map content.
             */
            inline bool isSourcesEmpty() const noexcept
            {
                return m_sources.empty();
            }

            //! Get the Groups map size of the manager.
            /** Get the Groups map size of the manager.
             @return     The groups map size.
             */
            inline size_t getNumberOfGroups() const noexcept
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

            //! The group constructor.
            /**	The group constructor allocates and initialize the member values for a source group.
             @param     index       The index of the group
             */
            Source::Group* createGroup (const size_t index)
            {
                return new Source::Group(this, index);
            }

            //! Add a Group to the manager.
            /** Add a Group to the map container of the manager.
             @param     group       The group to add.
             @return True if success, otherwise false.
             */
            inline bool addGroup (Group* group) noexcept
            {
                if (group && m_groups.find(group->m_index) == m_groups.end())
                {
                    if(group->getNumberOfSources() >= 2)
                    {
                        for(group_iterator it = m_groups.begin() ; it != m_groups.end(); ++it)
                        {
                            if(*group == *it->second)
                            {
                                return false;
                            }
                        }

                        m_groups[group->m_index] = group;
                        return true;
                    }
                }
                return false;
            }

            //! Remove a Group from the manager.
            /** Remove a Group from the map container of the manager.
             @param     index   The index of the group.
             */
            inline void removeGroup (const size_t index) noexcept
            {
                group_iterator it = m_groups.find(index);
                if(it != m_groups.end())
                {
                    std::map<size_t, Source*>& sources = it->second->getSources();
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

            //! Remove a Group from the manager with its sources.
            /** Remove a Group from the map container of the manager with its sources.
             @param     index   The index of the group.
             */
            void removeGroupWithSources (const size_t index) noexcept
            {
                group_iterator it = m_groups.find(index);
                if(it != m_groups.end())
                {
                    std::map<size_t, Source*>& sources = it->second->getSources();
                    source_iterator si = sources.begin();
                    while (si != sources.end())
                    {
                        source_iterator to = si;
                        ++to;
                        delete si->second;
                        sources.erase(si->first);
                        m_sources.erase(si->first);
                        si = to;
                    }

                    m_groups.erase(it->first);
                    delete it->second;
                    cleanEmptyGroup();
                }
            }

            //! Get one source.
            /** Get one source.
             @param     index   The index of the source.
             @return            A pointer on the source.
             */
            inline Source* getSource(const size_t index)
            {
                source_iterator it = m_sources.find(index);
                if(it != m_sources.end())
                {
                    return it->second;
                }
                return NULL;
            }

            //! Get the first const iterator of the Sources map of the manager.
            /** Get the first const iterator of the Sources map of the manager.
             @return        The first const iterator of the sources map.
             */
            inline const_source_iterator getFirstSource() const noexcept
            {
                return m_sources.begin();
            }

            //! Get the first iterator of the Sources map of the manager.
            /** Get the first iterator of the Sources map of the manager.
             @return        The const iterator of the sources map.
             */
            inline source_iterator getFirstSource() noexcept
            {
                return m_sources.begin();
            }

            //! Get the last const iterator of the Sources map of the manager.
            /** Get the last const iterator of the Sources map of the manager.
             @return        The last const iterator of the sources map.
             */
            inline const_source_iterator getLastSource() const noexcept
            {
                return m_sources.end();
            }

            //! Get the last iterator of the Sources map of the manager.
            /** Get the last iterator of the Sources map of the manager.
             @return        The last iterator of the sources map.
             */
            inline source_iterator getLastSource() noexcept
            {
                return m_sources.end();
            }

            //! Remove the groups which have less than 2 sources from the manager.
            /** Remove the groups which have less than 2 sources from the map container of the manager.
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

            //! Remove the group which have exactly the same sources from the manager.
            /** Remove the group which have exactly the same sources from the map container of the manager.
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
            inline Group* getGroup(const size_t index)
            {
                group_iterator it = m_groups.find(index);
                if(it != m_groups.end())
                {
                    return it->second;
                }
                return NULL;
            }

            //! Get the first const iterator of the Groups map of the manager.
            /** Get the first const iterator of the Groups map of the manager.
             @return        The first const iterator of the groups map.
             */
            inline const_group_iterator getFirstGroup() const noexcept
            {
                return m_groups.begin();
            }

            //! Get the first iterator of the Groups map of the manager.
            /** Get the first iterator of the Groups map of the manager.
             @return        The first iterator of the groups map.
             */
            inline group_iterator getFirstGroup() noexcept
            {
                return m_groups.begin();
            }

            //! Get the last const iterator of the Groups map of the manager.
            /** Get the last const iterator of the Groups map of the manager.
             @return        The last const iterator of the groups map.
             */
            inline const_group_iterator getLastGroup() const noexcept
            {
                return m_groups.end();
            }

            //! Get the last iterator of the Groups map of the manager.
            /** Get the last iterator of the Groups map of the manager.
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
			@see setCoordinatesCartesian
         */
        inline void setCoordinatesPolar(const double radius, const double azimuth)
		{
			setRadius(radius);
		    setAzimuth(azimuth);
       	}

		//! Set the position of the source with polar coordinates.
		/** Set the position of the source with polar coordinates.
			@param     radius			The radius of the source.
			@param     azimuth			The azimuth of the source.
            @param     elevation        The elevation of the source.
			@see setCoordinatesCartesian
         */
		inline void setCoordinatesPolar(const double radius, const double azimuth, const double elevation)
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
            m_radius = std::max(radius, (double)0.);
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
		    if(m_elevation > HOA_PI2)
		    {
		        m_azimuth = Math<double>::wrap_twopi(m_azimuth + HOA_PI);
                m_elevation = HOA_PI2 - (m_elevation - HOA_PI2);
		    }
		    else if(m_elevation < -HOA_PI2)
		    {
		        m_azimuth = Math<double>::wrap_twopi(m_azimuth + HOA_PI);
                m_elevation = -HOA_PI2 + (-m_elevation - HOA_PI2);
		    }
		    notifyCoordinates();
		}

		//! Set the position of the source with cartesian coordinates.
		/** Set the position of the source with cartesian coordinates.
            @param     abscissa		The abscissa of the source.
            @param     ordinate		The ordinate of the source.
         */
		inline void setCoordinatesCartesian(const double abscissa, const double ordinate)
		{
		    const double height = getHeight();
			setRadius(Math<double>::radius(abscissa, ordinate, height));
        	setAzimuth(Math<double>::azimuth(abscissa, ordinate, height));
        	setElevation(Math<double>::elevation(abscissa, ordinate, height));
        }

        //! Set the position of the source with cartesian coordinates.
		/** Set the position of the source with cartesian coordinates.
            @param     abscissa		The abscissa of the source.
            @param     ordinate		The ordinate of the source.
            @param     height		The height of the source.
         */
		inline void setCoordinatesCartesian(const double abscissa, const double ordinate, const double height)
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
			m_color[0]	=  clip(red, 0., 1.);
		    m_color[1]	=  clip(green, 0., 1.);
		    m_color[2]	=  clip(blue, 0., 1.);
		    m_color[3]	=  clip(alpha, 0., 1.);
        }

		//! Set the description of the source.
		/** Set the description of the source.
		 * @param     description		The text description of the source.
         */
		inline void setDescription(const std::string description)
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
        inline const size_t getIndex() const noexcept
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
		inline const std::string	getDescription() const noexcept
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

        //! Get the size of the Groups map of the source.
        /** Get the size of the Groups map of the source.
         @return    The group map size
         */
        inline size_t getNumberOfGroups() const noexcept
        {
            return m_groups.size();
        }

        //! Check if the Groups map of the source is empty.
        /** Check if the Groups map of the source is empty.
         @return    The state of the groups map content.
         */
        inline bool isGroupsEmpty() const noexcept
        {
            return m_groups.empty();
        }

        //! Get the Groups map of the source.
        /** Get the Groups map of the source.
         @return    A reference of the groups map.
         */
        inline std::map<size_t, Group*>& getGroups() noexcept
        {
            return m_groups;
        }

        //! The group class is used to control punctual sources.
        /** The group class is used to control punctual sources.
         */
        class Group
        {
        friend class Source;
        friend class Manager;

        private:
            const Manager*          m_manager;
            size_t                   m_index;
            std::map<size_t, Source*>     m_sources;
            std::string                  m_description;
            double			        m_color[4];
            double			        m_centroid_x;
            double			        m_centroid_y;
            double			        m_centroid_z;
            double                  m_maximum_radius;
            bool                    m_mute;
            bool                    m_subMute;

            //! The group constructor.
            /**	The group constructor allocates and initialize the member values for a source group.
             @param     manager		A pointer on a manager object
             @param     index       The index of the group
             */
            Group(const Manager* manager, const size_t index) : m_manager(manager)
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
            Group(const Group& other) :
            m_manager(other.m_manager),
            m_index(other.m_index),
            m_description(other.m_description),
            m_centroid_x(other.m_centroid_x),
            m_centroid_y(other.m_centroid_y),
            m_centroid_z(other.m_centroid_z),
            m_maximum_radius(other.m_maximum_radius),
            m_mute(other.m_mute),
            m_subMute(other.m_subMute)
            {
                memcpy(m_color, other.m_color, 4 * sizeof(double));
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
                size_t numberOfMutedSources = 0;
                for (source_iterator it = m_sources.begin() ; it != m_sources.end() ; it ++)
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
                    for (source_iterator it = m_sources.begin() ; it != m_sources.end() ; it ++)
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
                for (source_iterator it = m_sources.begin() ; it != m_sources.end() ; it ++)
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
                for (source_iterator it = m_sources.begin() ; it != m_sources.end() ; it ++)
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
                for (source_iterator it = m_sources.begin() ; it != m_sources.end() ; it ++)
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
                        for (source_iterator it = m_sources.begin() ; it != m_sources.end() ; it ++)
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
                        for (source_iterator it = m_sources.begin() ; it != m_sources.end() ; it ++)
                        {
                            double circleValue = sqrt(m_maximum_radius * m_maximum_radius - it->second->getOrdinate() * it->second->getOrdinate());
                            if(circleValue - it->second->getAbscissa() < refValue)
                                refValue = circleValue - it->second->getAbscissa();
                        }
                        if(abscissa > refValue)
                        {
                            abscissa = refValue;
                        }
                    }
                }
                for (source_iterator it = m_sources.begin() ; it != m_sources.end() ; it ++)
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
                        for (source_iterator it = m_sources.begin() ; it != m_sources.end() ; it ++)
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
                        for (source_iterator it = m_sources.begin() ; it != m_sources.end() ; it ++)
                        {
                            double circleValue = sqrt(m_maximum_radius * m_maximum_radius - it->second->getAbscissa() * it->second->getAbscissa());
                            if(circleValue - it->second->getOrdinate() < refValue)
                                refValue = circleValue - it->second->getOrdinate();
                        }
                        if(ordinate > refValue)
                        {
                            ordinate = refValue;
                        }
                    }
                }
                for (source_iterator it = m_sources.begin() ; it != m_sources.end() ; it ++)
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
                        for (source_iterator it = m_sources.begin() ; it != m_sources.end() ; it ++)
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
                        for (source_iterator it = m_sources.begin() ; it != m_sources.end() ; it ++)
                        {
                            double circleValue = sqrt(m_maximum_radius * m_maximum_radius - it->second->getAbscissa() * it->second->getAbscissa());
                            if(circleValue - it->second->getHeight() < refValue)
                                refValue = circleValue - it->second->getHeight();
                        }
                        if(height > refValue)
                        {
                            height = refValue;
                        }
                    }
                }
                for (source_iterator it = m_sources.begin() ; it != m_sources.end() ; it ++)
                {
                    it->second->setHeight(height + it->second->getHeight());
                }
            }

        public:

            //! The group destructor.
            /**	The group destructor free the memory.
             */
            ~Group() noexcept
            {
                for (source_iterator it = m_sources.begin() ; it != m_sources.end() ; it ++)
                {
                    m_sources[it->first]->removeGroup(m_index);
                }
                m_sources.clear();
            }

            //! Add a new Source to the group.
            /** Add a new Source to the map container of the group (and add the group to the source).
             @param     source  The source to add.
             @return            The state of the addition of the source.
             */
            inline bool addSource(Source* source) noexcept
            {
                if(source)
                {
                    source_iterator it = m_sources.find(source->getIndex());
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

            //! Remove a Source from the group.
            /** Remove a Source from the map container of the group.
             @param     index   The index of the source.
             */
            inline void removeSource(const size_t index) noexcept
            {
                m_sources.erase(index);
                computeCentroid();
            }

            //! Set the position of the group with polar coordinates.
            /** Set the position of the group with polar coordinates.
             @param     radius			The radius of the group.
             @param     azimuth			The azimuth of the group.
             @see setRelativeCoordinatesPolar, setCoordinatesCartesian
             */
            inline void setCoordinatesPolar(const double radius, const double azimuth)
            {
                if(m_maximum_radius >= 0)
                {
                    if(radius < -m_maximum_radius || radius > m_maximum_radius)
                        return;
                }
                const double elevation = getElevation();
                setCoordinatesCartesian(Math<double>::abscissa(radius, azimuth, elevation), Math<double>::ordinate(radius, azimuth, elevation));
            }

            //! Set the position of the group with polar coordinates.
            /** Set the position of the group with polar coordinates.
             @param     radius			The radius of the group.
             @param     azimuth			The azimuth of the group.
             @param     elevation			The elevation of the group.
             @see setRelativeCoordinatesPolar, setCoordinatesCartesian
             */
            inline void setCoordinatesPolar(const double radius, const double azimuth, const double elevation)
            {
                if(m_maximum_radius >= 0)
                {
                    if(radius < -m_maximum_radius || radius > m_maximum_radius)
                        return;
                }
                setCoordinatesCartesian(Math<double>::abscissa(radius, azimuth, elevation), Math<double>::ordinate(radius, azimuth, elevation), Math<double>::height(radius, azimuth, elevation));
            }

            //! Set the radius of the group.
            /** Set the radius of the group.
             @param     radius			The radius of the group.
             @see getRadius
             */
            inline void setRadius(const double radius)
            {
                if(m_maximum_radius >= 0)
                {
                    if(radius < -m_maximum_radius || radius > m_maximum_radius)
                        return;
                }
                const double azimuth = getAzimuth();
                const double elevation = getElevation();
                setCoordinatesCartesian(Math<double>::abscissa(radius, azimuth, elevation), Math<double>::ordinate(radius, azimuth, elevation), Math<double>::height(radius, azimuth, elevation));
            }

            //! Set the azimuth of the group.
            /** Set the azimuth of the group.
             @param     azimuth			The azimuth of the group.
             */
            inline void setAzimuth(const double azimuth)
            {
                const double radius = getRadius();
                const double elevation = getElevation();
                setCoordinatesCartesian(Math<double>::abscissa(radius, azimuth, elevation), Math<double>::ordinate(radius, azimuth, elevation), Math<double>::height(radius, azimuth, elevation));
            }

            //! Set the elevation of the group.
            /** Set the elevation of the group.
             @param     elevation			The elevation of the group.
             */
            inline void setElevation(const double elevation)
            {
                const double radius = getRadius();
                const double azimuth = getAzimuth();
                setCoordinatesCartesian(Math<double>::abscissa(radius, azimuth, elevation), Math<double>::ordinate(radius, azimuth, elevation), Math<double>::height(radius, azimuth, elevation));
            }

            //! Set the position of the group with cartesian coordinates.
            /** Set the position of the group with cartesian coordinates.
             @param     abscissa		The abscissa of the group.
             @param     ordinate		The ordinate of the group.
             */
            inline void setCoordinatesCartesian(double abscissa, double ordinate)
            {
                abscissa = abscissa - getAbscissa();
                ordinate = ordinate - getOrdinate();
                shiftCartesian(abscissa, ordinate);
                computeCentroid();
            }

            //! Set the position of the group with cartesian coordinates.
            /** Set the position of the group with cartesian coordinates.
             @param     abscissa		The abscissa of the group.
             @param     ordinate		The ordinate of the group.
             @param     height		The height of the group.
             */
            inline void setCoordinatesCartesian(double abscissa, double ordinate, double height)
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
                m_color[0]	=  clip(red, 0., 1.);
                m_color[1]	=  clip(green, 0., 1.);
                m_color[2]	=  clip(blue, 0., 1.);
                m_color[3]	=  clip(alpha, 0., 1.);
            }

            //! Set the description of the group.
            /** Set the description of the group.
             @param     description		The text description of the group.
             */
            inline void setDescription(const std::string description)
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
                for (source_iterator it = m_sources.begin() ; it != m_sources.end() ; it ++)
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
                double aRadiusOffset = radius - getRadius();
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
                azimuth = Math<double>::wrap_twopi(azimuth);
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
                elevation = Math<double>::wrap_twopi(elevation + HOA_PI2);
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
            inline const size_t getIndex() const noexcept
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
                return Math<double>::azimuth(m_centroid_x, m_centroid_y, m_centroid_z);
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
            inline const std::string	getDescription() const noexcept
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

            //! Get the size of the Sources map of the group.
            /** Get the size of the Sources map of the group.
             @return    The sources map size.
             */
            inline size_t getNumberOfSources() const noexcept
            {
                return m_sources.size();
            }

            //! Check if the Sources map of the group is empty.
            /** Check if the Sources map of the group is empty.
             @return    The state of the sources map content.
             */
            inline bool isSourcesEmpty() const noexcept
            {
                return m_sources.empty();
            }

            //! Get the Sources map of the group.
            /** Get the Sources map of the group.
             @return    A reference of the sources map.
             */
            inline std::map<size_t, Source*>& getSources() noexcept
            {
                return m_sources;
            }

            inline bool operator== (Group& other)
            {
                if (m_sources.size() == other.m_sources.size())
                {
                    size_t clones = 0;
                    for (source_iterator it = m_sources.begin() ; it != m_sources.end() ; it ++)
                    {
                        source_iterator ti = other.m_sources.find(it->first);
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
        size_t                m_index;
        double		         m_radius;
		double		         m_azimuth;
        double               m_elevation;
		double		         m_color[4];
        std::string          m_description;
        std::map<size_t, Group*>   m_groups;
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
		Source(const double maximumRadius, const size_t index, const double radius = 0., const double azimuth = 0., const double elevation = 0.)
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
            @param     other    The other source.
		 */
        Source(const Source& other) :
        m_index(other.m_index),
        m_radius(other.m_radius),
        m_azimuth(other.m_azimuth),
        m_elevation(other.m_elevation),
        m_description(other.m_description),
        m_maximum_radius(other.m_maximum_radius),
        m_mute(other.m_mute)
		{
            memcpy(m_color, other.m_color, 4 * sizeof(double));
      	}

      	//! The source destructor.
        /**	The source destructor free the memory.
         */
		~Source() noexcept
		{
        	for (group_iterator it = m_groups.begin() ; it != m_groups.end() ; it ++)
		    {
                m_groups[it->first]->removeSource(m_index);
            }
        	m_groups.clear();
        }

		//! Add a new group to the source.
		/** Add a new group to the map container of the source.
		 @param     group   The group to add.
		 @return            The state of the addition of the group.
		 */
        inline bool addGroup(Group* group) noexcept
        {
            group_iterator it = m_groups.find(group->getIndex());
            if(it == m_groups.end() && group)
            {
                m_groups[group->getIndex()] = group;
                return true;
            }
            return false;
        }

        //! Remove a group from the source.
        /** Remove a group from the map container of the source.
         @param     index   The index of the group.
         */
        inline void removeGroup(const size_t index) noexcept
        {
            m_groups.erase(index);
        }

        //! Call the groups of the source for each moving to compute their new position.
        /** Call the groups of the source for each moving to compute their new position.
         */
        inline void notifyCoordinates() noexcept
        {
            for (group_iterator it = m_groups.begin() ; it != m_groups.end() ; it ++)
		    {
                it->second->notifyCoordinates();
            }
        }

        //! Call the groups of the source for each change of its mute state to check their mute state.
        /** Call the groups of the source for each change of its mute state to check their mute state.
         */
        inline void notifyMute() noexcept
        {
            for (group_iterator it = m_groups.begin() ; it != m_groups.end() ; it ++)
		    {
                it->second->notifyMute();
            }
        }
        
        static inline T clip(const T& n, const T& lower, const T& upper)
        {
            return std::max(lower, std::min(n, upper));
        }
    };

}
