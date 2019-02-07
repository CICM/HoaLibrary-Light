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

#include "Hoa_Signal.hpp"

namespace hoa
{
    // ================================================================================ //
    // LINE //
    // ================================================================================ //
    
    template <typename T>
    class Line
    {
    public:
        
        //! @brief Constructor.
        Line() = default;

        //! @brief The Destructor.
        ~Line() = default;

        //! @brief Get the ramp value.
        //! @return The ramp value.
        inline size_t getRamp() const noexcept
        {
            return m_ramp;
        }

        //! @brief Get the current value.
        //! @return The current value.
        inline T getValue() const noexcept
        {
            return m_value_new;
        }

        //! @brief Set the ramp value.
        //! @param ramp The new value of the ramp.
        inline void setRamp(const size_t ramp) noexcept
        {
            m_ramp = std::max(ramp, (size_t)1);
        }

        //! @brief Set, linearly, the current value.
        //! @param value The new value of the current value.
        inline void setValue(const T value) noexcept
        {
            m_value_new = value;
            m_value_step = (m_value_new - m_value_old) / (T)m_ramp;
            m_counter = 0;
        }

        //! @brief Set, directly, the current value.
        //! @param value The new value of the current value.
        inline void setValueDirect(const T value) noexcept
        {
            m_value_old = m_value_new = value;
            m_value_step = 0.;
            m_counter = 0;
        }

        //! @brief This method performs the count of the virtual points of the line.
        //! @return The old value of the counter.
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
        
    private:
        
        T m_value_old = 0.;
        T m_value_new = 0.;
        T m_value_step = 0.;
        size_t m_counter = 0ul;
        size_t m_ramp = 0ul;
    };
    
    // ================================================================================ //
    // POLAR LINE //
    // ================================================================================ //

    template <Dimension D, typename T>
    class PolarLines;
    
    // ================================================================================ //
    // POLAR LINE 2D //
    // ================================================================================ //

    template <typename T>
    class PolarLines<Hoa2d, T>
    {
    public:
        
        //! @brief Constructor.
        //! @param sources The number of sources.
        PolarLines(size_t sources)
        : m_number_of_sources(sources)
        {
            const auto size = m_number_of_sources * 2;
            m_values_old = Signal<T>::alloc(size);
            m_values_new = Signal<T>::alloc(size);
            m_values_step = Signal<T>::alloc(size);
        }

        //! @brief Destructor.
        ~PolarLines()
        {
            Signal<T>::free(m_values_old);
            Signal<T>::free(m_values_new);
            Signal<T>::free(m_values_step);
        }

        //! @brief Get the number of sources.
        //! @return The number of sources.
        inline size_t getNumberOfSources() const noexcept
        {
            return m_number_of_sources;
        }

        //! @brief Get the ramp value.
        //! @return The ramp value.
        inline size_t getRamp() const noexcept
        {
            return m_ramp;
        }

        //! @brief Get the radius of a source.
        //! @param index    The index of the source.
        //! @return The radius of a source.
        inline T getRadius(const size_t index) const
        {
            return m_values_new[index];
        }

        //! @brief Get the azimuth of a source.
        //! @param index The index of the source.
        //! @return The azimuth of a source.
        inline T getAzimuth(const size_t index) const
        {
            return m_values_new[m_number_of_sources + index];
        }

        //! @brief Set the ramp value.
        //! @param ramp The new value of the ramp.
        inline void setRamp(const size_t ramp) noexcept
        {
            m_ramp = std::max(ramp, (size_t)1);
        }

        //! @brief Set, linearly, the radius of a source.
        //! @param index The index of the source.
        //! @param radius The new value of the radius.
        inline void setRadius(const size_t index, const T radius) noexcept
        {
            m_values_new[index]  = radius;
            m_values_step[index] = (m_values_new[index] - m_values_old[index]) / (T)m_ramp;
            m_counter = 0;
        }

        //! @brief Set, linearly, the azimuth of a source.
        //! @param index The index of the source.
        //! @param azimuth The new value of the azimuth.
        inline void setAzimuth(const size_t index, const T azimuth) noexcept
        {
            const auto idx = index + m_number_of_sources;
            m_values_new[idx] = math<T>::wrap_two_pi(azimuth);
            m_values_old[idx] = math<T>::wrap_two_pi(m_values_old[idx]);

            T distance;
            if(m_values_old[idx] > m_values_new[idx])
                distance = (m_values_old[idx] - m_values_new[idx]);
            else
                distance = (m_values_new[idx] - m_values_old[idx]);

            if(distance <= HOA_PI)
            {
                m_values_step[idx] = (m_values_new[idx] - m_values_old[idx]) / (T)m_ramp;
            }
            else
            {
                if(m_values_new[idx] > m_values_old[idx])
                {
                    m_values_step[idx] = ((m_values_new[idx] - HOA_2PI) - m_values_old[idx]) / (T)m_ramp;
                }
                else
                {
                    m_values_step[idx] = ((m_values_new[idx] + HOA_2PI) - m_values_old[idx]) / (T)m_ramp;
                }
            }
            
            m_counter = 0;
        }

        //! @brief Set, directly, the radius of a source.
        //! @param index The index of the source.
        //! @param radius The new value of the radius.
        inline void setRadiusDirect(const size_t index, const T radius) noexcept
        {
            m_values_old[index] = m_values_new[index] = radius;
            m_values_step[index] = 0.;
            m_counter = 0;
        }

        //! @brief Set, directly, the azimuth of a source.
        //! @param index The index of the source.
        //! @param azimuth The new value of the azimuth.
        inline void setAzimuthDirect(size_t index, const T azimuth) noexcept
        {
            const auto idx = index + m_number_of_sources;
            m_values_old[idx] = m_values_new[idx] = azimuth;
            m_values_step[idx] = 0.;
            m_counter = 0;
        }

        //! @brief This method performs the count of the virtual points of the line.
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
         
    private:
        
         const size_t m_number_of_sources;
         T* m_values_old = nullptr;
         T* m_values_new = nullptr;
         T* m_values_step = nullptr;
         size_t m_counter = 0ul;
         size_t m_ramp = 0ul;
    };
    
    // ================================================================================ //
    // POLAR LINE 3D //
    // ================================================================================ //

    template <typename T>
    class PolarLines<Hoa3d, T>
    {
    public:
        
        //! @brief Constructor.
        //! @param sources The number of sources.
        PolarLines(size_t sources)
        : m_number_of_sources(sources)
        {
            const auto size = m_number_of_sources * 3;
            m_values_old    = Signal<T>::alloc(size);
            m_values_new    = Signal<T>::alloc(size);
            m_values_step   = Signal<T>::alloc(size);
        }

        //! @brief Destructor.
        ~PolarLines()
        {
            Signal<T>::free(m_values_old);
            Signal<T>::free(m_values_new);
            Signal<T>::free(m_values_step);
        }

        //! @brief Get the number of sources.
        //! @return The number of sources.
        inline size_t getNumberOfSources() const noexcept
        {
            return m_number_of_sources;
        }

        //! @brief Get the ramp value.
        //! @return The ramp value.
        inline size_t getRamp() const noexcept
        {
            return m_ramp;
        }

        //! @brief Get the radius of a source.
        //! @param index The index of the source.
        //! @return The radius of a source.
        inline T getRadius(const size_t index) const
        {
            return m_values_new[index];
        }

        //! @brief Get the azimuth of a source.
        //! @param index The index of the source.
        //! @return The azimuth of a source.
        inline T getAzimuth(const size_t index) const
        {
            return m_values_new[m_number_of_sources + index];
        }

        //! @brief Get the elevation of a source.
        //! @param index The index of the source.
        //! @return The elevation of a source.
        inline T getElevation(const size_t index) const
        {
            return m_values_new[m_number_of_sources * 2 + index];
        }

        //! @brief Set the ramp value.
        //! @param ramp The new value of the ramp.
        inline void setRamp(const size_t ramp) noexcept
        {
            m_ramp = std::max(ramp, (size_t)1);
        }

        //! @brief Set, linearly, the radius of a source.
        //! @param index The index of the source.
        //! @param radius The new value of the radius.
        inline void setRadius(const size_t index, const T radius)
        {
            m_values_new[index]  = radius;
            m_values_step[index] = (m_values_new[index] - m_values_old[index]) / (T)m_ramp;
            m_counter = 0;
        }

        //! @brief Set, linearly, the azimuth of a source.
        //! @param index The index of the source.
        //! @param azimuth The new value of the azimuth.
        inline void setAzimuth(const size_t index, const T azimuth)
        {
            const auto idx = index + m_number_of_sources;
            m_values_new[idx] = math<T>::wrap_two_pi(azimuth);
            m_values_old[idx] = math<T>::wrap_two_pi(m_values_old[idx]);

            T distance;
            if(m_values_old[idx] > m_values_new[idx])
            {
                distance = (m_values_old[idx] - m_values_new[idx]);
            }
            else
            {
                distance = (m_values_new[idx] - m_values_old[idx]);
            }
            
            if(distance <= HOA_PI)
            {
                m_values_step[idx] = (m_values_new[idx] - m_values_old[idx]) / (T)m_ramp;
            }
            else
            {
                if(m_values_new[idx] > m_values_old[idx])
                {
                    m_values_step[idx] = ((m_values_new[idx] - HOA_2PI) - m_values_old[idx]) / (T)m_ramp;
                }
                else
                {
                    m_values_step[idx] = ((m_values_new[idx] + HOA_2PI) - m_values_old[idx]) / (T)m_ramp;
                }
            }
            
            m_counter = 0;
        }

        //! @brief Set, linearly, the elevation of a source.
        //! @param index The index of the source.
        //! @param elevation The new value of the elevation.
        inline void setElevation(const size_t index, const T elevation)
        {
            const auto idx = index + m_number_of_sources * 2;
            m_values_new[idx] = wrap_pi(elevation);
            m_values_old[idx] = wrap_pi(m_values_old[idx]);

            T distance;
            if(m_values_old[idx] > m_values_new[idx])
                distance = (m_values_old[idx] - m_values_new[idx]);
            else
                distance = (m_values_new[idx] - m_values_old[idx]);

            if(distance <= HOA_PI)
            {
                m_values_step[idx] = (m_values_new[idx] - m_values_old[idx]) / (T)m_ramp;
            }
            else
            {
                if(m_values_new[idx] > m_values_old[idx])
                {
                    m_values_step[idx] = ((m_values_new[idx] - HOA_2PI) - m_values_old[idx]) / (T)m_ramp;
                }
                else
                {
                    m_values_step[idx] = ((m_values_new[idx] + HOA_2PI) - m_values_old[idx]) / (T)m_ramp;
                }
            }
            m_counter = 0;
        }

        //! @brief Set, directly, the radius of a source.
        //! @param index The index of the source.
        //! @param radius The new value of the radius.
        inline void setRadiusDirect(const size_t index, const T radius) noexcept
        {
            m_values_old[index] = m_values_new[index] = radius;
            m_values_step[index] = 0.;
            m_counter = 0;
        }

        //! @brief Set, directly, the azimuth of a source.
        //! @param index    The index of the source.
        //! @param azimuth     The new value of the azimuth.
        inline void setAzimuthDirect(const size_t index, const T azimuth) noexcept
        {
            const size_t idx = index + m_number_of_sources;
            m_values_old[idx] = m_values_new[idx] = azimuth;
            m_values_step[idx] = 0.;
            m_counter = 0;
        }

        //! @brief Set, directly, the elevation of a source.
        //! @param index The index of the source.
        //! @param elevation The new value of the elevation.
        inline void setElevationDirect(const size_t index, const T elevation) noexcept
        {
            const size_t idx = index + m_number_of_sources * 2;
            m_values_old[idx] = m_values_new[idx] = elevation;
            m_values_step[idx] = 0.;
            m_counter = 0;
        }

        //! @brief This method performs the count of the virtual points of the line.
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
        
    private:
        
        const size_t m_number_of_sources;
        T* m_values_old = nullptr;
        T* m_values_new = nullptr;
        T* m_values_step = nullptr;
        size_t m_counter = 0ul;
        size_t m_ramp = 0ul;
    };
}
