/*
// Copyright (c) 2012-2015 Pierre Guillot, Eliott Paris & Thomas Le Meur CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_TOOLS_LIGHT
#define DEF_HOA_TOOLS_LIGHT

#include "Math.hpp"
#include "Signal.hpp"

namespace hoa
{
    template <typename T> class Line
    {
    private:
        T       m_value_old;
        T       m_value_new;
        T       m_value_step;
        ulong   m_counter;
        ulong   m_ramp;

    public:
        //! The line constructor.
        /**	The line constructor allocates and initialize the base classes.
         */
        Line() noexcept
        {
            ;
        }

        //! The destructor.
        /** The destructor free the memory.
         */
        ~Line()
        {
            ;
        }

        //! Get the ramp value.
        /** Get the ramp value.
        @return The ramp value.
         */
        inline ulong getRamp() const noexcept
        {
            return m_ramp;
        }

        //! Get the current value.
        /** Get the current value.
        @return The current value.
         */
        inline T getValue() const noexcept
        {
            return m_value_new;
        }

        //! Set the ramp value.
        /** Set the ramp value.
        @param ramp    The new value of the ramp.
         */
        inline void setRamp(const ulong ramp) noexcept
        {
            m_ramp = max(ramp, (ulong)1);
        }

        //! Set, linearly, the current value.
        /** Set, linearly, the current value.
        @param value    The new value of the current value.
         */
        inline void setValue(const T value) noexcept
        {
            m_value_new = value;
            m_value_step = (m_value_new - m_value_old) / (T)m_ramp;
            m_counter = 0;
        }

        //! Set, directly, the current value.
        /** Set, directly, the current value.
        @param value    The new value of the current value.
         */
        inline void setValueDirect(const T value) noexcept
        {
            m_value_old = m_value_new = value;
            m_value_step = 0.;
            m_counter = 0;
        }

        //! This method performs the count of the virtual points of the line.
        /** This method performs the count of the virtual points of the line.
        @return The old value of the counter.
         */
        inline T process() noexcept
        {
            m_value_old += m_value_step;
            if(m_counter++ >= m_ramp)
            {
                m_value_old  = m_value_new;
                m_value_step = 0.;
                m_counter    = 0;
            }
            return m_value_old;
        }
    };

    template <Dimension D, typename T> class PolarLines;

    template <typename T> class PolarLines<Hoa2d, T>
    {

    private:
        const ulong m_number_of_sources;
        T*      m_values_old;
        T*      m_values_new;
        T*      m_values_step;
        ulong   m_counter;
        ulong   m_ramp;

    public:
        //! The line constructor.
        /**	The line constructor allocates and initialize the base classes.
        @param numberOfSources  The number of sources.
         */
        PolarLines(ulong numberOfSources) noexcept :
        m_number_of_sources(numberOfSources)
        {
            m_values_old    = Signal<T>::alloc(m_number_of_sources * 2);
            m_values_new    = Signal<T>::alloc(m_number_of_sources * 2);
            m_values_step   = Signal<T>::alloc(m_number_of_sources * 2);
        }

        //! The destructor.
        /** The destructor free the memory.
         */
        ~PolarLines()
        {
            Signal<T>::free(m_values_old);
            Signal<T>::free(m_values_new);
            Signal<T>::free(m_values_step);
        }

        //! Get the number of sources.
        /** Get the number of sources.
        @return The number of sources.
         */
        inline ulong getNumberOfSources() const noexcept
        {
            return m_number_of_sources;
        }

        //! Get the ramp value.
        /** Get the ramp value.
        @return The ramp value.
         */
        inline ulong getRamp() const noexcept
        {
            return m_ramp;
        }

        //! Get the radius of a source.
        /** Get the radius of a source.
        @param index    The index of the source.
        @return The radius of a source.
         */
        inline T getRadius(const ulong index) const noexcept
        {
            return m_values_new[index];
        }

        //! Get the azimuth of a source.
        /** Get the azimuth of a source.
        @param index    The index of the source.
        @return The azimuth of a source.
         */
        inline T getAzimuth(const ulong index) const noexcept
        {
            return m_values_new[m_number_of_sources +index];
        }

        //! Set the ramp value.
        /** Set the ramp value.
        @param ramp    The new value of the ramp.
         */
        inline void setRamp(const ulong ramp) noexcept
        {
            m_ramp = max(ramp, (ulong)1);
        }

        //! Set, linearly, the radius of a source.
        /** Set, linearly, the radius of a source.
        @param index    The index of the source.
        @param radi    The new value of the radius.
         */
        inline void setRadius(const ulong index, const T radi) noexcept
        {
            m_values_new[index]  = radi;
            m_values_step[index] = (m_values_new[index] - m_values_old[index]) / (T)m_ramp;
            m_counter = 0;
        }

        //! Set, linearly, the azimuth of a source.
        /** Set, linearly, the azimuth of a source.
        @param index    The index of the source.
        @param azim     The new value of the azimuth.
         */
        inline void setAzimuth(const ulong index, const T azim) noexcept
        {
            m_values_new[index + m_number_of_sources] = Math<T>::wrap_twopi(azim);
            m_values_old[index + m_number_of_sources] = Math<T>::wrap_twopi(m_values_old[index + m_number_of_sources]);

            T distance;
            if(m_values_old[index + m_number_of_sources] > m_values_new[index + m_number_of_sources])
                distance = (m_values_old[index + m_number_of_sources] - m_values_new[index + m_number_of_sources]);
            else
                distance = (m_values_new[index + m_number_of_sources] - m_values_old[index + m_number_of_sources]);

            if(distance <= HOA_PI)
            {
                m_values_step[index + m_number_of_sources] = (m_values_new[index + m_number_of_sources] - m_values_old[index + m_number_of_sources]) / (T)m_ramp;
            }
            else
            {
                if(m_values_new[index + m_number_of_sources] > m_values_old[index + m_number_of_sources])
                {
                    m_values_step[index + m_number_of_sources] = ((m_values_new[index + m_number_of_sources] - HOA_2PI) - m_values_old[index + m_number_of_sources]) / (T)m_ramp;
                }
                else
                {
                    m_values_step[index + m_number_of_sources] = ((m_values_new[index + m_number_of_sources] + HOA_2PI) - m_values_old[index + m_number_of_sources]) / (T)m_ramp;
                }
            }
            m_counter = 0;
        }

        //! Set, directly, the radius of a source.
        /** Set, directly, the radius of a source.
        @param index    The index of the source.
        @param radi     The new value of the radius.
         */
        inline void setRadiusDirect(const ulong index, const T radi) noexcept
        {
            m_values_old[index] = m_values_new[index] = radi;
            m_values_step[index] = 0.;
            m_counter = 0;
        }

        //! Set, directly, the azimuth of a source.
        /** Set, directly, the azimuth of a source.
        @param index    The index of the source.
        @param azim     The new value of the azimuth.
         */
        inline void setAzimuthDirect(ulong index, const T azim) noexcept
        {
            m_values_old[index + m_number_of_sources] = m_values_new[index + m_number_of_sources] = azim;
            m_values_step[index + m_number_of_sources] = 0.;
            m_counter = 0;
        }

        //! This method performs the count of the virtual points of the line.
        /** This method performs the count of the virtual points of the line.
         */
        void process(T* vector) noexcept
        {
            Signal<T>::add(m_number_of_sources * 2, m_values_step, m_values_old);
            if(m_counter++ >= m_ramp)
            {
                Signal<T>::copy(m_number_of_sources * 2, m_values_new, m_values_old);
                Signal<T>::clear(m_number_of_sources * 2, m_values_step);
                m_counter    = 0;
            }
            Signal<T>::copy(m_number_of_sources * 2, m_values_old, vector);
        }
    };

    template <typename T> class PolarLines<Hoa3d, T>
    {

    private:
        const ulong m_number_of_sources;
        T*      m_values_old;
        T*      m_values_new;
        T*      m_values_step;
        ulong   m_counter;
        ulong   m_ramp;

    public:
        //! The line constructor.
        /**	The line constructor allocates and initialize the base classes.
        @param numberOfSources  The number of sources.
         */
        PolarLines(ulong numberOfSources) noexcept :
        m_number_of_sources(numberOfSources)
        {
            m_values_old    = Signal<T>::alloc(m_number_of_sources * 3);
            m_values_new    = Signal<T>::alloc(m_number_of_sources * 3);
            m_values_step   = Signal<T>::alloc(m_number_of_sources * 3);
        }

        //! The destructor.
        /** The destructor free the memory.
         */
        ~PolarLines()
        {
            Signal<T>::free(m_values_old);
            Signal<T>::free(m_values_new);
            Signal<T>::free(m_values_step);
        }

        //! Get the number of sources.
        /** Get the number of sources.
        @return The number of sources.
         */
        inline ulong getNumberOfSources() const noexcept
        {
            return m_number_of_sources;
        }

        //! Get the ramp value.
        /** Get the ramp value.
        @return The ramp value.
         */
        inline ulong getRamp() const noexcept
        {
            return m_ramp;
        }

        //! Get the radius of a source.
        /** Get the radius of a source.
        @param index    The index of the source.
        @return The radius of a source.
         */
        inline T getRadius(const ulong index) const noexcept
        {
            return m_values_new[index];
        }

        //! Get the azimuth of a source.
        /** Get the azimuth of a source.
        @param index    The index of the source.
        @return The azimuth of a source.
         */
        inline T getAzimuth(const ulong index) const noexcept
        {
            return m_values_new[m_number_of_sources + index];
        }

        //! Get the elevation of a source.
        /** Get the elevation of a source.
        @param index    The index of the source.
        @return The elevation of a source.
         */
        inline T getElevation(const ulong index) const noexcept
        {
            return m_values_new[m_number_of_sources * 2 + index];
        }

        //! Set the ramp value.
        /** Set the ramp value.
        @param ramp    The new value of the ramp.
         */
        inline void setRamp(const ulong ramp) noexcept
        {
            m_ramp = max(ramp, (ulong)1);
        }

        //! Set, linearly, the radius of a source.
        /** Set, linearly, the radius of a source.
        @param index    The index of the source.
        @param radi     The new value of the radius.
         */
        inline void setRadius(const ulong index, const T radi) noexcept
        {
            m_values_new[index]  = radi;
            m_values_step[index] = (m_values_new[index] - m_values_old[index]) / (T)m_ramp;
            m_counter = 0;
        }

        //! Set, linearly, the azimuth of a source.
        /** Set, linearly, the azimuth of a source.
        @param index    The index of the source.
        @param azim     The new value of the azimuth.
         */
        inline void setAzimuth(const ulong index, const T azim) noexcept
        {
            m_values_new[index + m_number_of_sources] = Math<T>::wrap_twopi(azim);
            m_values_old[index + m_number_of_sources] = Math<T>::wrap_twopi(m_values_old[index + m_number_of_sources]);

            T distance;
            if(m_values_old[index + m_number_of_sources] > m_values_new[index + m_number_of_sources])
                distance = (m_values_old[index + m_number_of_sources] - m_values_new[index + m_number_of_sources]);
                else
                    distance = (m_values_new[index + m_number_of_sources] - m_values_old[index + m_number_of_sources]);

                    if(distance <= HOA_PI)
                    {
                        m_values_step[index + m_number_of_sources] = (m_values_new[index + m_number_of_sources] - m_values_old[index + m_number_of_sources]) / (T)m_ramp;
                    }
                    else
                    {
                        if(m_values_new[index + m_number_of_sources] > m_values_old[index + m_number_of_sources])
                        {
                            m_values_step[index + m_number_of_sources] = ((m_values_new[index + m_number_of_sources] - HOA_2PI) - m_values_old[index + m_number_of_sources]) / (T)m_ramp;
                        }
                        else
                        {
                            m_values_step[index + m_number_of_sources] = ((m_values_new[index + m_number_of_sources] + HOA_2PI) - m_values_old[index + m_number_of_sources]) / (T)m_ramp;
                        }
                    }
            m_counter = 0;
        }

        //! Set, linearly, the elevation of a source.
        /** Set, linearly, the elevation of a source.
        @param index    The index of the source.
        @param elev     The new value of the elevation.
         */
        inline void setElevation(const ulong index, const T elev) noexcept
        {
            m_values_new[index + m_number_of_sources * 2] = Math<T>::wrap_pi(elev);
            m_values_old[index + m_number_of_sources * 2] = Math<T>::wrap_pi(m_values_old[index + m_number_of_sources * 2]);

            T distance;
            if(m_values_old[index + m_number_of_sources * 2] > m_values_new[index + m_number_of_sources * 2])
                distance = (m_values_old[index + m_number_of_sources * 2] - m_values_new[index + m_number_of_sources * 2]);
            else
                distance = (m_values_new[index + m_number_of_sources * 2] - m_values_old[index + m_number_of_sources * 2]);

            if(distance <= HOA_PI)
            {
                m_values_step[index + m_number_of_sources * 2] = (m_values_new[index + m_number_of_sources * 2] - m_values_old[index + m_number_of_sources * 2]) / (T)m_ramp;
            }
            else
            {
                if(m_values_new[index + m_number_of_sources * 2] > m_values_old[index + m_number_of_sources * 2])
                {
                    m_values_step[index + m_number_of_sources * 2] = ((m_values_new[index + m_number_of_sources * 2] - HOA_2PI) - m_values_old[index + m_number_of_sources * 2]) / (T)m_ramp;
                }
                else
                {
                    m_values_step[index + m_number_of_sources * 2] = ((m_values_new[index + m_number_of_sources * 2] + HOA_2PI) - m_values_old[index + m_number_of_sources * 2]) / (T)m_ramp;
                }
            }
            m_counter = 0;
        }

        //! Set, directly, the radius of a source.
        /** Set, directly, the radius of a source.
        @param index    The index of the source.
        @param radi    The new value of the radius.
         */
        inline void setRadiusDirect(const ulong index, const T radi) noexcept
        {
            m_values_old[index] = m_values_new[index] = radi;
            m_values_step[index] = 0.;
            m_counter = 0;
        }

        //! Set, directly, the azimuth of a source.
        /** Set, directly, the azimuth of a source.
        @param index    The index of the source.
        @param azim     The new value of the azimuth.
         */
        inline void setAzimuthDirect(const ulong index, const T azim) noexcept
        {
            m_values_old[index + m_number_of_sources] = m_values_new[index + m_number_of_sources] = azim;
            m_values_step[index + m_number_of_sources] = 0.;
            m_counter = 0;
        }

        //! Set, directly, the elevation of a source.
        /** Set, directly, the elevation of a source.
        @param index    The index of the source.
        @param elev     The new value of the elevation.
         */
        inline void setElevationDirect(const ulong index, const T elev) noexcept
        {
            m_values_old[index + m_number_of_sources * 2] = m_values_new[index + m_number_of_sources * 2] = elev;
            m_values_step[index + m_number_of_sources * 2] = 0.;
            m_counter = 0;
        }

        //! This method performs the count of the virtual points of the line.
        /** This method performs the count of the virtual points of the line.
         */
        void process(T* vector) noexcept
        {
            Signal<T>::add(m_number_of_sources * 3, m_values_step, m_values_old);
            if(m_counter++ >= m_ramp)
            {
                Signal<T>::copy(m_number_of_sources * 3, m_values_new, m_values_old);
                Signal<T>::clear(m_number_of_sources * 3, m_values_step);
                m_counter    = 0;
            }
            Signal<T>::copy(m_number_of_sources * 3, m_values_old, vector);
        }
    };
}

#endif


