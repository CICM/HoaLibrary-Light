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
#include "Hoa_Math.hpp"

namespace hoa
{
    // ================================================================================ //
    // PLANEWAVE //
    // ================================================================================ //
    
    template <Dimension D, typename T>
    class ProcessorPlanewaves;
    
    //! @brief The class owns basic plane waves informations.
    //! @details The class allows to retrieves several informations about the plane waves: the
    //! the polar coordinates, azimuth and the elevation and cartesian coodinates,  abscissa,
    //! the ordinate and the height.
    template<Dimension D, typename T>
    class Planewave
    {
    public:

        //! @brief The plane wave polar constructor.
        //! @param index       The index must be at least 1.
        //! @param azimuth     The azimuth \f$\theta\f$.
        //! @param elevation   The elevation \f$\varphi\f$.
        Planewave(size_t index, T azimuth, T elevation) noexcept
        : m_index(index)
        , m_azimuth(math<T>::wrap_two_pi(azimuth))
        , m_elevation(math<T>::wrap_pi(elevation))
        {}

        //! @brief The plane wave cartesian constructor.
        //! @param index       The index must be at least 1.
        //! @param abscissa    The abscissa \f$x\f$.
        //! @param ordinate    The ordinate \f$y\f$.
        //! @param height      The height \f$z\f$.
        Planewave(size_t index, T abscissa, T ordinate, T height) noexcept
        : Planewave(index,
                    xyz2azimuth(abscissa, ordinate, height),
                    xyz2elevation(abscissa, ordinate, height))
        {}

        //! @brief Destructor.
		virtual ~Planewave() noexcept {}

        //! @brief Returns the index of the plane wave.
        inline size_t getIndex() const noexcept { return m_index; }

        //! @brief Returns the azimuth of the planewave.
        inline T getAzimuth() const noexcept { return m_azimuth; }
        
        //! @brief Returns the elevation of the planewave.
        inline T getElevation() const noexcept { return m_elevation; }
        
        //! @brief Returns the abscissa of the plane wave.
        inline T getAbscissa() const noexcept { return ae2abscissa(m_azimuth, m_elevation); }
        
        //! @brief Returns the ordinate of the plane wave.
        inline T getOrdinate() const noexcept { return ae2ordinate(m_azimuth, m_elevation); }

        //! @brief Returns the height of the plane wave.
        inline T getHeight() const noexcept { return ae2height(m_azimuth, m_elevation); }

        //! @brief Sets the azimuth of the plane wave.
        //! @param azimuth The azimuth \f$\theta\f$.
        inline void setAzimuth(const T azimuth) noexcept { m_azimuth = math<T>::wrap_two_pi(azimuth); }

        //! @brief Set the elevation of the plane wave.
        //! @param elevation The elevation \f$\varphi\f$.
        inline void setElevation(const T elevation) noexcept { m_elevation = math<T>::wrap_pi(elevation); }

        //! @brief Returns the name of the plane wave.
        virtual std::string getName() const noexcept
        {
            std::ostringstream ostr;
            ostr <<  "Planewave " << getIndex() << " " << (m_azimuth / HOA_2PI * 360.) << "°";
            if(D == Hoa3d)
            {
                ostr << " " << (m_elevation / HOA_2PI * 360.) << "°";
            }
            
            return ostr.str();
        }
        
        //! @brief Compares the azimuth of two plane waves.
        static bool compare_azimuth(Planewave const& i, Planewave const& j) noexcept
        {
            return i.m_azimuth < j.m_azimuth;
        }
        
    private:
        
        friend class ProcessorPlanewaves<D, T>;
        
        static inline T xyz2azimuth(const T x, const T y, const T z = 0.)
        {
            hoa_unused(z);
            return (x == T(0) && y == T(0)) ? T(0) : math<T>::wrap_two_pi(atan2(y, x) - HOA_PI2);
        }
        
        static inline T xyz2elevation(const T x, const T y, const T z = 0.)
        {
            return (z == T(0)) ? T(0) : math<T>::wrap_pi(asin(z / sqrt(x*x + y*y + z*z)));
        }
        
        static inline T ae2abscissa(const T a, const T e)
        {
            return std::cos(a + T(HOA_PI2)) * std::cos(e);
        }
        
        static inline T ae2ordinate(const T a, const T e)
        {
            return std::sin(a + T(HOA_PI2)) * std::cos(e);
        }
        
        static inline T ae2height(const T a, const T e)
        {
            hoa_unused(a);
            return std::sin(e);
        }
        
        static inline void rotate(Planewave const& p, const T xr, const T yr, const T zr, T& xo, T& yo, T& zo) noexcept
        {
            const T xi = p.getAbscissa();
            const T yi = p.getOrdinate();
            const T zi = p.getHeight();
            // Around x axis
            const T cosxr = std::cos(xr);
            const T sinxr = std::sin(xr);
            // xx
            const T yx = cosxr * yi + -sinxr * zi;
            const T zx = sinxr * yi + cosxr  * zi;
            
            // Around y axis
            const T cosyr = std::cos(yr);
            // yy
            const T sinyr = std::sin(yr);
            const T xy = cosyr  * xi + sinxr * zx;
            zo = -sinyr * xi + cosyr * zx;
            
            // Around z axis
            const T coszr = std::cos(zr);
            const T sinzr = std::sin(zr);
            xo = coszr * xy + -sinzr * yx;
            yo = sinzr * xy + coszr * yx;
            // zz
        }
        
    private:
        
        size_t  m_index = 0ul;
        T       m_azimuth = 0.;
        T       m_elevation = 0.;
    };
}
