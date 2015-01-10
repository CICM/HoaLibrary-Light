/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_SOURCE
#define DEF_SOURCE

#include "HoaMath.hpp"

namespace hoa
{
	//! The source
    /** The source store and manage source informations like the position, color, mute...
     */
	template <typename T> class Source
	{
    public:
        struct Color
        {
            T red;
            T green;
            T blue;
            T alpha;
        };
        
        class Group;
        class Manager;
        static T  m_maximum_radius;
	private:
        const ulong m_index;
		T		m_radius;
		T		m_azimuth;
        T       m_elevation;
        bool    m_mute;
		Color   m_color;
        string  m_description;
		
        vector<weak_ptr<Group>> m_groups;
        
        inline void clean()
        {
            for(auto it = m_groups.begin(); it != m_groups.end(); )
            {
                shared_ptr<Group> group = (*it).lock();
                if(group)
                {
                    ++it;
                }
                else
                {
                    it = m_groups.erase(it);
                }
            }
        }
        
        inline void notify()
        {
            for(auto it = m_groups.begin(); it != m_groups.end(); )
            {
                shared_ptr<Group> group = (*it).lock();
                if(group)
                {
                    group->computeCentroid();
                    ++it;
                }
                else
                {
                    it = m_groups.erase(it);
                }
            }
        }
        
	public:

        Source(const ulong index, const T radius, const T  azimuth, const T  elevation) noexcept :
        m_index(index),
        m_mute(false),
        m_color({0.2, 0.2, 0.2, 1.}),
        m_description("")
        {
            setRadius(radius);
            setAzimuth(azimuth);
            setElevation(elevation);
        }
		
		~Source()
        {
            m_groups.clear();
        }
        
        inline ulong getIndex() const noexcept
        {
            return m_index;
        }
        
        inline void setMute(const bool state) noexcept
        {
            m_mute = state;
        }

		inline void setPolar(const T radius, const T azimuth, const T elevation) noexcept
        {
            if(!m_maximum_radius || radius <= m_maximum_radius)
            {
                m_radius = max(radius, (T)0.);
            }
            m_azimuth = wrap_twopi(azimuth);
            m_elevation = wrap(elevation, -HOA_PI, HOA_PI);
            if(m_elevation > HOA_PI2)
            {
                m_azimuth = wrap_twopi(m_azimuth + HOA_PI);
                m_elevation = HOA_PI2 - (elevation - HOA_PI2);
            }
            else if(m_elevation < -HOA_PI2)
            {
                m_azimuth = wrap_twopi(m_azimuth + HOA_PI);
                m_elevation = -HOA_PI2 + (-elevation + HOA_PI2);
            }
            notify();
        }

		inline void setRadius(const T radius) noexcept
        {
            if(!m_maximum_radius || radius <= m_maximum_radius)
            {
                m_radius = max(radius, (T)0.);
            }
            notify();
        }

		inline void setAzimuth(const T azimuth) noexcept
        {
            m_azimuth = wrap_twopi(azimuth);
            notify();
        }

		inline void setElevation(const T elevation) noexcept
        {
            m_elevation = wrap(elevation, -HOA_PI, HOA_PI);
            if(m_elevation > HOA_PI2)
            {
                m_azimuth = wrap_twopi(m_azimuth + HOA_PI);
                m_elevation = HOA_PI2 - (elevation - HOA_PI2);
            }
            else if(m_elevation < -HOA_PI2)
            {
                m_azimuth = wrap_twopi(m_azimuth + HOA_PI);
                m_elevation = -HOA_PI2 + (-elevation + HOA_PI2);
            }
            notify();
        }

		inline void setCartesian(const T  abscissa, const T  ordinate, const T  height) noexcept
        {
            setPolar(radius(abscissa, ordinate, height), azimuth(abscissa, ordinate, height), elevation(abscissa, ordinate, height));
        }

		inline void setAbscissa(const T abscissa) noexcept
        {
            const T ordinate = getOrdinate();
            const T height = getHeight();
            setPolar(radius(abscissa, ordinate, height), azimuth(abscissa, ordinate, height), elevation(abscissa, ordinate, height));
        }

		inline void setOrdinate(const T ordinate) noexcept
        {
            const T abscissa = getAbscissa();
            const T height = getHeight();
            setPolar(radius(abscissa, ordinate, height), azimuth(abscissa, ordinate, height), elevation(abscissa, ordinate, height));
        }

		inline void setHeight(const T height) noexcept
        {
            const T abscissa = getAbscissa();
            const T ordinate = getOrdinate();
            setPolar(radius(abscissa, ordinate, height), azimuth(abscissa, ordinate, height), elevation(abscissa, ordinate, height));
        }

		inline void setColor(Color const& color) noexcept
        {
            m_color = color;
        }

		inline void setDescription(string const& description) noexcept
        {
            m_description = description;
        }
		
		void setGroup(shared_ptr<Group> group) noexcept
        {
            clean();
            if(group)
            {
                auto it = find(m_groups.begin(), m_groups.end(), group);
                if(it == m_groups.end())
                {
                    m_groups.push_back(group);
                }
            }
        }

		void removeGroup(shared_ptr<Group> group)
        {
            clean();
            if(group)
            {
                auto it = find(m_groups.begin(), m_groups.end(), group);
                if(it != m_groups.end())
                {
                    m_groups.erase(it);
                }
            }
        }

		inline T getRadius() const noexcept
        {
            return m_radius;
        }

		inline T getAzimuth() const noexcept
        {
            return m_azimuth;
        }

		inline T getElevation() const noexcept
        {
            return m_elevation;
        }

		inline T getAbscissa() const noexcept
        {
            return abscissa(m_radius, m_azimuth, m_elevation);
        }

		inline T getOrdinate() const noexcept
        {
            return ordinate(m_radius, m_azimuth, m_elevation);
        }

		inline T getHeight() const noexcept
        {
            return height(m_radius, m_azimuth, m_elevation);
        }

		inline Color getColor() const noexcept
        {
            return m_color;
        }

		inline string getDescription()	const noexcept
        {
            return m_description;
        }

		inline ulong getNumberOfGroups() const noexcept
        {
            clean();
            return m_groups.size();
        }

		inline shared_ptr<Group> getGroup(const ulong index) const noexcept
        {
            clean();
            return m_groups[index].lock();
        }

		inline bool isOwned(shared_ptr<Group> group) const noexcept
        {
            clean();
            if(group)
            {
                return m_groups.find(group) != m_groups.end();
            }
            else
            {
                return false;
            }
        }
    
		inline bool getMute() const noexcept
        {
            return m_mute;
        }
	};
    template <typename T> T Source<T>::m_maximum_radius = 0.;
    
    template<typename T> class Source<T>::Group
    {
    private:
        T				m_centroid_x;
        T				m_centroid_y;
        T				m_centroid_z;
        bool            m_mute;
        Color			m_color;
        string          m_description;
        vector<weak_ptr<Source<T>>> m_sources;
        
        inline void computeCentroid() noexcept
        {
            m_centroid_x = m_centroid_y = m_centroid_z = 0.;
            if(m_sources.size())
            {
                for(auto it = m_sources.begin(); it != m_sources.end(); ++it)
                {
                    shared_ptr<Source<T>> source = it.lock();
                    if(source)
                    {
                        m_centroid_x += source->getAbscissa();
                        m_centroid_y += source->getOrdinate();
                        m_centroid_z += source->getHeight();
                    }
                }
                m_centroid_x /= m_sources.size();
                m_centroid_y /= m_sources.size();
                m_centroid_z /= m_sources.size();
            }
        }
        
        inline void shiftPolar(const T radius, const T azimuth, const T elevation) noexcept
        {
            shiftRadius(radius);
            shiftAzimuth(azimuth);
            shiftElevation(elevation);
        }
        
        inline void shiftRadius(T radius)
        {
            if(m_maximum_radius)
            {
                if(radius < 0.)
                {
                    T refRadius = m_maximum_radius;
                    for(auto it = m_sources.begin(); it != m_sources.end(); ++it)
                    {
                        shared_ptr<Source<T>> source = it.lock();
                        if(source)
                        {
                            if(source->getRadius() < refRadius)
                            {
                                refRadius = source->getRadius();
                            }
                        }
                        
                    }
                    if(radius + refRadius < 0.)
                    {
                        radius = -refRadius;
                    }
                }
                else
                {
                    T refRadius = -m_maximum_radius;
                    for(auto it = m_sources.begin(); it != m_sources.end(); ++it)
                    {
                        shared_ptr<Source<T>> source = it.lock();
                        if(source)
                        {
                            if(source->getRadius() > refRadius)
                            {
                                refRadius = source->getRadius();
                            }
                        }
                        
                    }
                    if(radius + refRadius > m_maximum_radius)
                    {
                        radius = m_maximum_radius - refRadius;
                    }
                    
                }
            }
            for(auto it = m_sources.begin(); it != m_sources.end(); ++it)
            {
                shared_ptr<Source<T>> source = it.lock();
                if(source)
                {
                    source->setRadius(radius + source->getRadius());
                }
            }
        }
        
        inline void shiftAzimuth(const T azimuth)
        {
            for(auto it = m_sources.begin(); it != m_sources.end(); ++it)
            {
                shared_ptr<Source<T>> source = it.lock();
                if(source)
                {
                    source->setAzimuth(azimuth + source->getAzimuth());
                }
            }
        }
        
        inline void shiftElevation(const T elevation)
        {
            for(auto it = m_sources.begin(); it != m_sources.end(); ++it)
            {
                shared_ptr<Source<T>> source = it.lock();
                if(source)
                {
                    source->setElevation(azimuth + source->getElevation());
                }
            }
        }
        
        inline void shiftCartesian(double abscissa, double ordinate, double height);
        inline void shiftAbscissa(double abscissa);
        inline void shiftOrdinate(double ordinate);
        inline void shiftHeight(double ordinate);
        
    public:
        
        //! The source group constructor.
        /**	The source group constructor allocates and initialize the member values for a source group.
         */
        Group() noexcept :
        m_mute(false),
        m_color({0.2, 0.2, 0.2, 1.}),
        m_description("")
        {
            computeCentroid();
        }
        
        //! The source group destructor.
        /**	The source group destructor free the memory.
         */
        ~Group()
        {
            m_sources.clear();
        }
        
        //! Set the position of the group with polar coordinates.
        /**
         @param     radius			The radius of the group.
         @param     azimuth			The azimuth of the group.
         @param     elevation			The elevation of the group.
         */
        inline void setCoordinatesPolar(const T radius, const T azimuth, const T elevation) noexcept
        {
            setCoordinatesCartesian(abscissa(radius, azimuth, elevation), ordinate(radius, azimuth, elevation));
        }
        
        //! Set the radius of the group.
        /**
         @param     radius			The radius of the group.
         */
        inline void setRadius(const T radius) noexcept
        {
            setCoordinatesCartesian(abscissa(radius, getAzimuth(), getElevation()), ordinate(radius, getAzimuth(), getElevation()), height(radius, getAzimuth(), getElevation()));
        }
        
        //! Set the azimuth of the group.
        /**
         @param     azimuth			The azimuth of the group.
         */
        inline void setAzimuth(const T azimuth) noexcept
        {
            setCoordinatesCartesian(abscissa(getRadius(), azimuth, getElevation()), ordinate(getRadius(), azimuth, getElevation()), height(getRadius(), azimuth, getElevation()));
        }
        
        //! Set the elevation of the group.
        /**
         @param     elevation			The elevation of the group.
         */
        inline void setElevation(const T elevation) noexcept
        {
            setCoordinatesCartesian(abscissa(getRadius(), getAzimuth(), elevation), ordinate(getRadius(), getAzimuth(), elevation), height(getRadius(), getAzimuth(), elevation));

        }
        
        //! Set the position of the group with cartesians coordinates.
        /**
         @param     abscissa		The abscissa of the group.
         @param     ordinate		The ordinate of the group.
         @param     height		The height of the group.
         */
        inline void setCoordinatesCartesian(const T abscissa, const T ordinate, const T height) noexcept
        {
            const T nabs = abscissa - getAbscissa();
            const T nord = ordinate - getOrdinate();
            const T nhei = height - getHeight();
            shiftAbscissa(nabs);
            shiftOrdinate(nord);
            shiftHeight(nhei);
        }
        
        //! Set the abscissa of the group.
        /**
         @param     abscissa		The abscissa of the group.
         */
        inline void setAbscissa(const T abscissa) noexcept
        {
            const T nabs = abscissa - getAbscissa();
            shiftAbscissa(nabs);
        }
        
        //! Set the ordinate of the group.
        /**
         @param     ordinate		The ordinate of the group.
         */
        inline void setOrdinate(const T ordinate) noexcept
        {
            const T nord = ordinate - getOrdinate();
            shiftOrdinate(nord);
        }
        
        //! Set the height of the group.
        /**
         @param     height		The height of the group.
         */
        inline void setHeight(const T height) noexcept
        {
            const T nhei = height - getHeight();
            shiftHeight(nhei);
        }
        
        //! Set the color of the group.
        /**
         @param     color The color
         */
        inline void setColor(Color const& color) noexcept
        {
            m_color = color;
        }
        
        //! Set the description of the group.
        /**
         @param     description		The text description of the group.
         */
        inline void setDescription(string const& description) noexcept
        {
            m_description = description;
        }
        
        //! Store a new source in this group
        /**
         @param     sourceIndex		The index of the source to store.
         */
        inline void addSource(shared_ptr<Source<T>> source)
        {
            
        }
        
        //! Remove a source from this group
        /**
         @param     sourceIndex		The index of the source to store.
         */
        inline void removeSource(long sourceIndex);
        
        //! Set the mute state of the group
        /**
         @param     state		The mute state of the group.
         @see getMute
         */
        inline void setMute(const bool state) noexcept
        {
            m_mute = state;
            for(auto it = m_sources.begin(); it != m_sources.end(); ++it)
            {
                shared_ptr<Source<T>> source = it.lock();
                if(source)
                {
                    source->setMute(m_mute);
                }
            }
        }
        
        //! Set the position of the group with relative polar coordinates.
        /**
         @param     radius			The relative radius value.
         @param     azimuth			The relative azimuth value.
         @param     elevation		The relative elevation value.
         */
        inline void setRelativeCoordinatesPolar(const T radius, const T azimuth, const T elevation) noexcept
        {
            shiftRadius(max(radius, 0.) - getRadius());
            shiftAzimuth(azimuth + HOA_PI2 - getAzimuth());
            shiftElevation(elevation + HOA_PI2 - getElevation());
        }
        
        //! Set the radius of the group with a relative value.
        /**
         @param     radius			The relative radius value.
         @see setCoordinatesPolar, setRelativeRadius
         */
        inline void setRelativeRadius(const T radius) noexcept
        {
            shiftRadius(max(radius, 0.) - getRadius());
        }
        
        //! Set the azimuth of the group with a relative value.
        /**
         @param     azimuth			The relative azimuth value.
         @see setCoordinatesPolar, setRelativeRadius
         */
        inline void setRelativeAzimuth(const T azimuth) noexcept
        {
            shiftAzimuth(azimuth + HOA_PI2 - getAzimuth());
        }
        
        //! Set the elevation of the group with a relative value.
        /**
         @param     elevation			The relative elevation value.
         @see setCoordinatesPolar, setRelativeRadius
         */
        inline void setRelativeElevation(const T elevation) noexcept
        {
            shiftElevation(elevation + HOA_PI2 - getElevation());
        }
        
        //! Get the radius of the group.
        /**
         @return		The radius of the group.
         @see setRadius, setCoordinatesPolar
         */
        inline T getRadius() const noexcept
        {
            return radius(m_centroid_x, m_centroid_y, m_centroid_z);
        }
        
        //! Get the azimuth of the group.
        /**
         @return		The azimuth of the group.
         @see setAzimuth, setCoordinatesPolar
         */
        inline T getAzimuth() const noexcept
        {
            return azimuth(m_centroid_x, m_centroid_y, m_centroid_z) + HOA_PI2;
        }
        
        //! Get the elevation of the group.
        /**
         @return		The elevation of the group.
         @see setAzimuth, setCoordinatesPolar
         */
        inline T getElevation() const noexcept
        {
            return elevation(m_centroid_x, m_centroid_y, m_centroid_z);
        }
        
        //! Get the abscissa of the group.
        /**
         @return		The abscissa of the group.
         @see setAbscissa, setCoordinatesCartesian
         */
        inline T getAbscissa() const noexcept
        {
            return m_centroid_x;
        }
        
        //! Get the ordinate of the group.
        /**
         @return		The ordinate of the group.
         @see setOrdinate, setCoordinatesCartesian
         */
        inline T getOrdinate()		const noexcept
        {
            return m_centroid_y;
        }
        
        //! Get the height of the group.
        /**
         @return		The height of the group.
         @see setOrdinate, setCoordinatesCartesian
         */
        inline T getHeight() const noexcept
        {
            return m_centroid_z;
        }
        
        //! Get the color of the group.
        /**
         @return		The rgba color of the group as an array of 4 double.
         @see setColor
         */
        inline Color getColor() const noexcept
        {
            return m_color;
        }
        
        //! Get the description of the group.
        /**
         @return		The description of the group.
         @see setDescription
         */
        inline string getDescription() const noexcept
        {
            return m_description;
        }
        
        //! Get the number of sources owned by this group.
        /**
         @return		The number of sources.
         @see setDescription
         */
        inline ulong getNumberOfSources() const noexcept
        {
            return m_sources.size();
        }
        
        //! Get the mute state of the group
        /**
         @return		The mute state of the group.
         @see setMute
         */
        inline bool getMute() const noexcept
        {
            return m_mute;
        }
        
        //! Get the the index of a source stored at a particular index by the group.
        /**
         @param			index			The index of the source.
         @return			The index of the source if it exists, -1 otherwise.
         */
        shared_ptr<Source<T>> getSource(const ulong index)	const noexcept
        {
            return m_sources[index];
        }
    };
    
    template<typename T> class Source<T>::Manager
    {
        
    private:
        map<ulong, shared_ptr<Source>> m_sources;
        map<ulong, shared_ptr<Group>>  m_groups;
        T                              m_zoom;
    public:
        
        //! The source manager constructor.
        /**	The source manager constructor allocates and initialize the member values.
         *
         * @param     maximumRadius		The maximum radius the sources or groups in the source manager could have
         * @param     existence			The existence state of the source manager.
         */
        Manager() noexcept :
        m_zoom(1.)
        {
            
        }
        
        //! The source manager destructor free the memory.
        ~Manager()
        {
            m_sources.clear();
            m_groups.clear();
        }
        
        //! Make a copy of this SourcesManager into an other
        void copyTo(Manager* other);
        
        //! Clear all the sources and groups.
        void clearAll() noexcept
        {
            m_sources.clear();
            m_groups.clear();
        }
        
        //! Set the zoom factor.
        /** Set the zoom factor between 0 and 1.
         
         @param     zoom		The zoom factor.
         */
        inline void setZoom(const T zoom) noexcept
        {
            m_zoom = clip(zoom, 1. / m_maximum_radius, 1.);;
        }
        
        //! Get the maximum index of sources.
        /**
         @return		The maximum index of source.
         */
        inline long getMaximumIndexOfSource() const noexcept
        {
            if(!m_sources.empty())
            {
                m_sources.rbegin()->first;
            }
            else
            {
                return -1;
            }
        }
        
        //! Get the number of sources actually managed.
        /**
         @return		The number of sources.
         */
        inline ulong getNumberOfSources() const noexcept
        {
            m_sources.size();
        }
        
        //! Get the maximum index of the sources actually managed.
        /**
         @return		The maximum index.
         */
        inline long getMaximumIndexOfGroup() const noexcept
        {
            if(!m_groups.empty())
            {
                m_groups.rbegin()->first;
            }
            else
            {
                return -1;
            }
        }
        
        //! Get the number of groups actually managed.
        /**
         @return		The number of group.
         */
        inline ulong getNumberOfGroups()  const noexcept
        {
            m_groups.size();
        }
        
        //! Get the zoom factor value.
        /**
         @return		The zoom factor value.
         */
        T getZoom() const noexcept
        {
            return m_zoom;
        }
        
        /* ------------------------------------------------------------------------ */
        /* ------------------------------- Sources -------------------------------- */
        /* ------------------------------------------------------------------------ */
        
        //! Add a new source with polar coordinates.
        /**
         * @param     radius			The radius of the source.
         * @param     azimuth			The azimuth of the source.
         * @param     elevation			The elevation of the source.
         * @see sourceNewCartesian
         */
        shared_ptr<Source<T>> addSourcePolar(const T radius, const T azimuth, const T elevation)
        {
            ulong index = 0;
            for(auto it = m_sources.begin(); it != m_sources.end(); ++it)
            {
                if(index != it->first)
                {
                    shared_ptr<Source<T>> source = make_shared<Source<T>>(index, radius, azimuth, elevation);
                    m_sources[index] = source;
                    return source;
                }
                ++index;
            }
            shared_ptr<Source<T>> source = make_shared<Source<T>>(m_sources.size(), radius, azimuth, elevation);
            m_sources[m_sources.size()] = source;
            return source;
            
        }
        
        //! Add a new source with cartesian coordinates.
        /**
         * @param     abscissa			The abscissa of the source.
         * @param     ordinate			The ordinate of the source.
         * @param     height			The height of the source.
         * @see sourceNewPolar
         */
        shared_ptr<Source<T>> addSourceCartesian(const T abscissa, const T ordinate, const T height)
        {
            return addSourcePolar(radius(abscissa, ordinate, height), azimuth(abscissa, ordinate, height), elevation(abscissa, ordinate, height));
        }
        
        //! Remove a source.
        /** Remove a source.
         * @param     index				The index of the source.
         */
        inline void removeSource(const ulong index) noexcept
        {
            auto it = m_sources.find(index);
            if(it != m_sources.end())
            {
                m_sources.erase(it);
            }
        }
        
        //! Retrieve a source.
        /** Retrieve a source.
         * @param     index				The index of the source.
         */
        inline shared_ptr<Source<T>> getSource(const ulong index) const noexcept
        {
            auto it = m_sources.find(index);
            if(it != m_sources.end())
            {
                return it->second;
            }
            else
            {
                return nullptr;
            }
        }
        
        //! Remove a source.
        /** Remove a source.
         * @param     index				The index of the source.
         */
        inline void removeGroup(const ulong index) noexcept
        {
            auto it = m_groups.find(index);
            if(it != m_groups.end())
            {
                m_groups.erase(it);
            }
        }
        
        //! Retrieve a group.
        /** Retrieve a group.
         * @param     index				The index of the group.
         */
        inline shared_ptr<Source<T>::Group> getGroup(const ulong index) const noexcept
        {
            auto it = m_groups.find(index);
            if(it != m_groups.end())
            {
                return it->second;
            }
            else
            {
                return nullptr;
            }
        }
    };
}
#endif