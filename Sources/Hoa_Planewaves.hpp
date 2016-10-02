/*
// Copyright (c) 2012-2016 CICM - Universite Paris 8 - Labex Arts H2H.
// Authors :
// 2012: Pierre Guillot, Eliott Paris & Julien Colafrancesco.
// 2012-2015: Pierre Guillot & Eliott Paris.
// 2015: Pierre Guillot & Eliott Paris & Thomas Le Meur (Light version)
// 2016: Pierre Guillot & Eliott Paris.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_PLANEWAVES_LIGHT
#define DEF_HOA_PLANEWAVES_LIGHT

#include "Hoa_Math.hpp"
#include "Hoa_Signal.hpp"

namespace hoa
{
    //! @brief The class owns basic plane waves informations.
    //! @details The class allows to retrieves several informations about the plane waves: the
    //! the polar coordinates, azimuth and the elevation and cartesian coodinates,  abscissa,
    //! the ordinate and the height.
    template<Dimension D, typename T> class Planewave
    {
    public:

        //! @brief The plane wave polar constructor.
        //! @param index       The index must be at least 1.
        //! @param azimuth     The azimuth \f$\theta\f$.
        //! @param elevation   The elevation \f$\varphi\f$.
        Planewave(size_t index, T azimuth, T elevation) hoa_noexcept :
        m_index(index), m_azimuth(wrap_twopi(azimuth)), m_elevation(wrap_pi(elevation)) {}

        //! @brief The plane wave cartesian constructor.
        //! @param index       The index must be at least 1.
        //! @param abscissa    The abscissa \f$x\f$.
        //! @param ordinate    The ordinate \f$y\f$.
        //! @param height      The height \f$z\f$.
        Planewave(size_t index, T abscissa, T ordinate, T height) hoa_noexcept :
        m_index(index), m_azimuth(xyz2azimuth(abscissa, ordinate, height)), m_elevation(xyz2elevation(abscissa, ordinate, height)) {}

        //! @brief The plane wave destructor.
		~Planewave() hoa_noexcept {}

        //! Returns the index of the plane wave.
        inline size_t getIndex() const hoa_noexcept { return m_index; }

        //! @brief Returns the azimuth of the planewave.
        inline T getAzimuth() const hoa_noexcept { return m_azimuth; }
        
        //! @brief Returns the elevation of the planewave.
        inline T getElevation() const hoa_noexcept { return m_elevation; }
        
        //! @brief Returns the abscissa of the plane wave.
        inline T getAbscissa() const hoa_noexcept { return ae2abscissa(m_azimuth, m_elevation); }
        
        //! @brief Returns the ordinate of the plane wave.
        inline T getOrdinate() const hoa_noexcept { return ae2ordinate(m_azimuth, m_elevation); }

        //! @brief Returns the height of the plane wave.
        inline T getHeight() const hoa_noexcept { return ae2height(m_azimuth, m_elevation); }

        //! @brief Sets the azimuth of the plane wave.
        //! @param azimuth The azimuth \f$\theta\f$.
        inline void setAzimuth(const T azimuth) hoa_noexcept { m_azimuth = wrap_twopi(azimuth); }

        //! @brief Set the elevation of the plane wave.
        //! @param elevation The elevation \f$\varphi\f$.
        inline void setElevation(const T elevation) hoa_noexcept { m_elevation = wrap_pi(elevation); }

        //! @brief Returns the name of the plane wave.
        virtual std::string getName() const hoa_noexcept {
            std::ostringstream ostr;
            ostr <<  "Planewave " << getIndex() << " " << (m_azimuth / HOA_2PI * 360.) << "°";
            if(D == Hoa3d) {
                ostr << " " << (m_elevation / HOA_2PI * 360.) << "°"; }
            return ostr.str();
        }
        
        //! @brief Compares the azimuth of two plane waves.
        static bool compare_azimuth(Planewave const& i, Planewave const& j) hoa_noexcept { return i.m_azimuth < j.m_azimuth; }
        
    private:
        
        static inline T xyz2azimuth(const T x, const T y, const T z = 0.) {
            return (x == T(0) && y == T(0)) ? T(0) : wrap_twopi(atan2(y, x) - HOA_PI2); }
        
        static inline T xyz2elevation(const T x, const T y, const T z = 0.) {
            return (z == T(0)) ? T(0) : wrap_pi(asin(z / sqrt(x*x + y*y + z*z))); }
        
        static inline T ae2abscissa(const T a, const T e) {
            return std::cos(a + T(HOA_PI2)) * std::cos(e); }
        
        static inline T ae2ordinate(const T a, const T e) {
            return std::sin(a + T(HOA_PI2)) * std::cos(e); }
        
        static inline T ae2height(const T a, const T e) {
            return std::sin(e); }
        
        static inline T wrap_pi(T const value) {
            if(value < T(-HOA_PI))
                return value + T(HOA_2PI);
            else if(value >= T(HOA_PI))
                return value - T(HOA_2PI);
            return value;
        }
        
        static inline T wrap_twopi(T const value) {
            if(value < 0.)
                return value + T(HOA_2PI);
            else if(value >= T(HOA_2PI))
                return value - T(HOA_2PI);
            return value;
        }
        
        size_t  m_index;
        T       m_azimuth;
        T       m_elevation;
    };
}

#endif


