/*
// Copyright (c) 2012-2015 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
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
        Line() noexcept
        {
            ;
        }
        
        ~Line()
        {
            ;
        }
        
        inline ulong getRamp() const noexcept
        {
            return m_ramp;
        }
        
        inline T getValue() const noexcept
        {
            return m_value_new;
        }
        
        inline void setRamp(const ulong ramp) noexcept
        {
            m_ramp = max(ramp, (ulong)1);
        }
        
        inline void setValue(const T value) noexcept
        {
            m_value_new = value;
            m_value_step = (m_value_new - m_value_old) / (T)m_ramp;
            m_counter = 0;
        }
        
        inline void setValueDirect(const T value) noexcept
        {
            m_value_old = m_value_new = value;
            m_value_step = 0.;
            m_counter = 0;
        }
        
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
        
        PolarLines(ulong numberOfSources) noexcept :
        m_number_of_sources(numberOfSources)
        {
            m_values_old    = new T[m_number_of_sources * 2];
            m_values_new    = new T[m_number_of_sources * 2];
            m_values_step   = new T[m_number_of_sources * 2];
        }
        
        ~PolarLines()
        {
            delete [] m_values_old;
            delete [] m_values_new;
            delete [] m_values_step;
        }
        
        inline ulong getNumberOfSources() const noexcept
        {
            return m_number_of_sources;
        }
        
        inline ulong getRamp() const noexcept
        {
            return m_ramp;
        }
        
        inline T getRadius(const ulong index) const noexcept
        {
            return m_values_new[index];
        }
        
        inline T getAzimuth(const ulong index) const noexcept
        {
            return m_values_new[m_number_of_sources +index];
        }
        
        inline void setRamp(const ulong ramp) noexcept
        {
            m_ramp = max(ramp, (ulong)1);
        }
        
        inline void setRadius(const ulong index, const T radius) noexcept
        {
            m_values_new[index]  = radius;
            m_values_step[index] = (m_values_new[index] - m_values_old[index]) / (T)m_ramp;
            m_counter = 0;
        }
        
        inline void setAzimuth(const ulong index, const T azimuth) noexcept
        {
            m_values_new[index + m_number_of_sources] = Math<T>::wrap_twopi(azimuth);
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
        
        inline void setRadiusDirect(const ulong index, const T radius) noexcept
        {
            m_values_old[index] = m_values_new[index] = radius;
            m_values_step[index] = 0.;
            m_counter = 0;
        }
        
        inline void setAzimuthDirect(ulong index, const T azimuth) noexcept
        {
            m_values_old[index + m_number_of_sources] = m_values_new[index + m_number_of_sources] = azimuth;
            m_values_step[index + m_number_of_sources] = 0.;
            m_counter = 0;
        }
        
        void process(T* vector) noexcept
        {
            Signal<T>::vector_add(m_number_of_sources * 2, m_values_step, m_values_old);
            if(m_counter++ >= m_ramp)
            {
                Signal<T>::vector_copy(m_number_of_sources * 2, m_values_new, m_values_old);
                Signal<T>::vector_clear(m_number_of_sources * 2, m_values_step);
                m_counter    = 0;
            }
            Signal<T>::vector_copy(m_number_of_sources * 2, m_values_old, vector);
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
        
        PolarLines(ulong numberOfSources) noexcept :
        m_number_of_sources(numberOfSources)
        {
            m_values_old    = new T[m_number_of_sources * 3];
            m_values_new    = new T[m_number_of_sources * 3];
            m_values_step   = new T[m_number_of_sources * 3];
        }
        
        ~PolarLines()
        {
            delete [] m_values_old;
            delete [] m_values_new;
            delete [] m_values_step;
        }
        
        inline ulong getNumberOfSources() const noexcept
        {
            return m_number_of_sources;
        }
        
        inline ulong getRamp() const noexcept
        {
            return m_ramp;
        }
        
        inline T getRadius(const ulong index) const noexcept
        {
            return m_values_new[index];
        }
        
        inline T getAzimuth(const ulong index) const noexcept
        {
            return m_values_new[m_number_of_sources + index];
        }
        
        inline T getElevation(const ulong index) const noexcept
        {
            return m_values_new[m_number_of_sources * 2 + index];
        }
        
        inline void setRamp(const ulong ramp) noexcept
        {
            m_ramp = max(ramp, (ulong)1);
        }
        
        inline void setRadius(const ulong index, const T radius) noexcept
        {
            m_values_new[index]  = radius;
            m_values_step[index] = (m_values_new[index] - m_values_old[index]) / (T)m_ramp;
            m_counter = 0;
        }
        
        inline void setAzimuth(const ulong index, const T azimuth) noexcept
        {
            m_values_new[index + m_number_of_sources] = Math<T>::wrap_twopi(azimuth);
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
        
        inline void setElevation(const ulong index, const T elevation) noexcept
        {
            m_values_new[index + m_number_of_sources * 2] = Math<T>::wrap_pi(elevation);
            m_values_old[index + m_number_of_sources * 2] = Math<T>::wrap_pi(m_values_old[index + m_number_of_sources * 2]);
            
            double distance;
            if(m_values_old[index + m_number_of_sources * 2] > m_values_new[index + m_number_of_sources * 2])
                distance = (m_values_old[index + m_number_of_sources * 2] - m_values_new[index + m_number_of_sources * 2]);
            else
                distance = (m_values_new[index + m_number_of_sources * 2] - m_values_old[index + m_number_of_sources * 2]);
            
            if(distance <= HOA_PI)
            {
                m_values_step[index + m_number_of_sources * 2] = (m_values_new[index + m_number_of_sources * 2] - m_values_old[index + m_number_of_sources * 2]) / (double)m_ramp;
            }
            else
            {
                if(m_values_new[index + m_number_of_sources * 2] > m_values_old[index + m_number_of_sources * 2])
                {
                    m_values_step[index + m_number_of_sources * 2] = ((m_values_new[index + m_number_of_sources * 2] - HOA_2PI) - m_values_old[index + m_number_of_sources * 2]) / (double)m_ramp;
                }
                else
                {
                    m_values_step[index + m_number_of_sources * 2] = ((m_values_new[index + m_number_of_sources * 2] + HOA_2PI) - m_values_old[index + m_number_of_sources * 2]) / (double)m_ramp;
                }
            }
            m_counter = 0;
        }
        
        inline void setRadiusDirect(const ulong index, const T radius) noexcept
        {
            m_values_old[index] = m_values_new[index] = radius;
            m_values_step[index] = 0.;
            m_counter = 0;
        }
        
        inline void setAzimuthDirect(const ulong index, const T azimuth) noexcept
        {
            m_values_old[index + m_number_of_sources] = m_values_new[index + m_number_of_sources] = azimuth;
            m_values_step[index + m_number_of_sources] = 0.;
            m_counter = 0;
        }
        
        inline void setElevationDirect(const ulong index, const T elevation) noexcept
        {
            m_values_old[index + m_number_of_sources * 2] = m_values_new[index + m_number_of_sources * 2] = elevation;
            m_values_step[index + m_number_of_sources * 2] = 0.;
            m_counter = 0;
        }
        
        void process(T* vector) noexcept
        {
            Signal<T>::vector_add(m_number_of_sources * 3, m_values_step, m_values_old);
            if(m_counter++ >= m_ramp)
            {
                Signal<T>::vector_copy(m_number_of_sources * 3, m_values_new, m_values_old);
                Signal<T>::vector_clear(m_number_of_sources * 3, m_values_step);
                m_counter    = 0;
            }
            Signal<T>::vector_copy(m_number_of_sources * 3, m_values_old, vector);
        }
    };
}

#endif


