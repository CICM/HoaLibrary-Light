/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_TOOLS_LIGHT
#define DEF_HOA_TOOLS_LIGHT

#include "HoaMath.hpp"

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
        
        inline T getValue(const ulong index) const noexcept
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
        
        inline void setRadiusDirect(const T value) noexcept
        {
            m_value_old = m_value_new = value;
            m_value_step = 0.;
            m_counter = 0;
        }
        
        inline T process() noexcept
        {
            if(m_counter < m_ramp)
            {
                m_value_old += m_value_step;
            }
            else if(m_counter++ == m_ramp)
            {
                m_value_old = m_value_new;
                m_value_step = 0.;
                m_counter    = 0;
            }
            return m_value_old;
        }
    };
    
    template <typename T> class PolarLines
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
            m_values_new[index + m_number_of_sources] = wrap_twopi(azimuth);
            m_values_old[index + m_number_of_sources] = wrap_twopi(m_values_old[index + m_number_of_sources]);
            
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
        
        void process(float* vector) noexcept
        {
            cblas_saxpy(m_number_of_sources * 2, 1., m_values_step, 1, m_values_old, 1);
            if(m_counter++ >= m_ramp)
            {
                cblas_scopy(m_number_of_sources * 2, m_values_new, 1, m_values_old, 1);
                memset(m_values_step, 0, sizeof(float) * m_number_of_sources * 2);
                m_counter    = 0;
            }
            cblas_scopy(m_number_of_sources * 2, m_values_old, 1, vector, 1);
        }
        
        void process(double* vector) noexcept
        {
            cblas_daxpy(m_number_of_sources * 2, 1., m_values_step, 1, m_values_old, 1);
            if(m_counter++ >= m_ramp)
            {
                cblas_dcopy(m_number_of_sources * 2, m_values_new, 1, m_values_old, 1);
                memset(m_values_step, 0, sizeof(double) * m_number_of_sources * 2);
                m_counter    = 0;
            }
            cblas_dcopy(m_number_of_sources * 2, m_values_old, 1, vector, 1);
        }
    };
}

#endif


