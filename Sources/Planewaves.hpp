/*
// Copyright (c) 2012-2015 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_PLANEWAVES_LIGHT
#define DEF_HOA_PLANEWAVES_LIGHT

#include "Math.hpp"
#include "Signal.hpp"

namespace hoa
{
    //! The planewave class owns basic position informations.
    /** The planewave allows to retrieves informations about its position, the azimuth and the elevation in polar or the abscissa, the ordinate and the height cartesian.
     */
    template<Dimension D, typename T> class Planewave
    {
    public:
        
        //! The planewave constructor.
        /** The planewave constructor allocates and initializes the generale member values depending on an index and a polar coordinate with an azimuth \f$\theta\f$ and an elevation \f$\varphi\f$ in radian assuming that the radius \f$\rho\f$ is always equal to \f$1\f$.
         @param index       The index must be at least 1.
         @param azimuth     The azimuth \f$\theta\f$.
         @param elevation   The elevation \f$\varphi\f$.
         */
        Planewave(const ulong index, const T azimuth, const T elevation) noexcept = 0;
        
        //! The planewave constructor.
        /** The planewave constructor allocates and initializes the generale member values depending on an index and a cartesian coordinate with an abscissa \f$x\f$, an ordinate \f$y\f$ and an height \f$z\f$ assuming that the point is normalized over the unit circle or unit the sphere.
         @param index       The index must be at least 1.
         @param abscissa    The abscissa \f$x\f$.
         @param ordinate    The ordinate \f$y\f$.
         @param height      The height \f$z\f$.
         */
        Planewave(const ulong index, const T abscissa, const T ordinate, const T height) noexcept = 0;
        
        //! The planewave destructor.
        /** The planewave destructor free the memory.
         */
        virtual ~Planewave() noexcept;
        
        //! Get the index of the planewave.
        /** The method returns the index \f$i\f$ of the planewave.
         @return     The index.
         */
        virtual ulong getIndex() const noexcept;
        
        //! Get the azimuth of the planewave.
        /** The method returns the azimuth \f$\theta\f$ between \f$0\f$ and \f$2\pi\f$ of the planewave. The result will consider a full rotation around the axes x, y, and z. The order the rotation is x, z then y.
         @param x_axe The rotation around the x axe in radian.
         @param y_axe The rotation around the y axe in radian.
         @param z_axe The rotation around the z axe in radian.
         @return     The azimuth.
         */
        virtual T getAzimuth(const T x_axe, const T y_axe, const T z_axe) const noexcept;
        
        //! Set the azimuth of the planewave.
        /** The method sets the azimuth \f$\theta\f$ of the planewave.
         @param azimuth The azimuth \f$\theta\f$.
         @return     The azimuth.
         */
        virtual void setAzimuth(const T azimuth) noexcept;
        
        //! Get the elevation of the planewave.
        /** The method returns the elevation \f$\varphi\f$ between \f$-\pi\f$ and \f$\pi\f$ of the planewave. The result will consider a full rotation around the axes x, y, and z. The order the rotation is x, z then y.
         @param x_axe The rotation around the x axe in radian.
         @param y_axe The rotation around the y axe in radian.
         @param z_axe The rotation around the z axe in radian.
         @return     The elevation.
         */
        virtual T getElevation(const T x_axe, const T y_axe, const T z_axe) const noexcept;
        
        //! Set the elevation of the planewave.
        /** The method sets the elevation \f$\varphi\f$ of the planewave.
         @param elevation The elevation \f$\varphi\f$.
         @return     The elevation.
         */
        virtual void setElevation(const T elevation) noexcept;
        
        //! Get the abscissa of the planewave.
        /** The method returns the abscissa \f$x\f$ between \f$-1\f$ and \f$1\f$ of the planewave. The result will consider a full rotation around the axes x, y, and z. The order the rotation is x, z then y.
         @param x_axe The rotation around the x axe in radian.
         @param y_axe The rotation around the y axe in radian.
         @param z_axe The rotation around the z axe in radian.
         @return     The abscissa.
         */
        virtual T getAbscissa(const T x_axe, const T y_axe, const T z_axe) const noexcept;
        
        //! Get the ordinate of the planewave.
        /** The method returns the ordinate \f$y\f$ between \f$-1\f$ and \f$1\f$ of the planewave. The result will consider a full rotation around the axes x, y, and z. The order the rotation is x, z then y.
         @param x_axe The rotation around the x axe in radian.
         @param y_axe The rotation around the y axe in radian.
         @param z_axe The rotation around the z axe in radian.
         @return     The ordinate.
         */
        virtual T getOrdinate(const T x_axe, const T y_axe, const T z_axe) const noexcept;
        
        //! Get the height of the planewave.
        /** The method returns the height \f$z\f$ between \f$-1\f$ and \f$1\f$ of the planewave. The result will consider a full rotation around the axes x, y, and z. The order the rotation is x, z then y.
         @param x_axe The rotation around the x axe in radian.
         @param y_axe The rotation around the y axe in radian.
         @param z_axe The rotation around the z axe in radian.
         @return     The height.
         */
        virtual T getHeight(const T x_axe, const T y_axe, const T z_axe) const noexcept;
        
        //! Get the name of the planewave.
        /** The method returns the name \f$planewave_{index}\f$ of the planewave.
         @return     The name.
         */
        virtual string getName() const noexcept;
    };
    
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    
    template<typename T> class Planewave<Hoa2d, T>
    {
    private:
        ulong m_index;
        T     m_azimuth;
    public:
        
        //! The planewave constructor.
        /** The planewave constructor allocates and initializes the generale member values depending on an index and a polar coordinate with an azimuth \f$\theta\f$ in radian assuming that the elevation \f$\varphi\f$ is always equal to \f$0\f$ and the radius \f$\rho\f$ is always equal to \f$1\f$.
         @param index       The index must be at least 1.
         @param azimuth     The azimuth \f$\theta\f$.
         @param elevation   The elevation \f$\varphi\f$ (ignored).
         */
        Planewave(const ulong _index, const T _azimuth, const T _elevation) noexcept :
        m_index(_index),
        m_azimuth(_azimuth)
        {
            ;
        }
        
        //! The planewave constructor.
        /** The planewave constructor allocates and initializes the generale member values depending on an index and a cartesian coordinate with an abscissa \f$x\f$, an ordinate \f$y\f$ and an height \f$z\f$ assuming that the point is normalized over the unit circle or unit the sphere and the height \f$z\f$ is always equal to \f$0\f$.
         @param index       The index must be at least 1.
         @param abscissa    The abscissa \f$x\f$.
         @param ordinate    The ordinate \f$y\f$.
         @param height      The height \f$z\f$ (ignored).
         */
        Planewave(const ulong _index, const T _abscissa, const T _ordinate, const T _height) noexcept :
        m_index(_index),
        m_azimuth(Math<T>::azimuth(_abscissa, _ordinate, 0.))
        {
            ;
        }
        
        //! The planewave destructor.
        /** The planewave destructor free the memory.
         */
        ~Planewave() noexcept
        {
            ;
        }
        
        //! Get the index of the planewave.
        /** The method returns the index \f$i\f$ of the planewave.
         @return     The index.
         */
        inline ulong getIndex() const noexcept
        {
            return m_index;
        }
        
        //! Get the azimuth of the planewave.
        /** The method returns the azimuth \f$\theta\f$ between \f$0\f$ and \f$2\pi\f$ of the planewave. The result will consider a rotation around the axes x.
         @param x_axe The rotation around the x axe in radian (ignored).
         @param y_axe The rotation around the y axe in radian (ignored).
         @param z_axe The rotation around the z axe in radian.
         @return     The azimuth.
         */
        inline T getAzimuth(const T x_axe, const T y_axe, const T z_axe) const noexcept
        {
            return Math<T>::wrap_twopi(z_axe + m_azimuth);
        }
        
        //! Set the azimuth of the planewave.
        /** The method sets the azimuth \f$\theta\f$ of the planewave.
         @param azimuth The azimuth \f$\theta\f$.
         @return     The azimuth.
         */
        inline void setAzimuth(const T azimuth) noexcept
        {
            m_azimuth = Math<T>::wrap_twopi(azimuth);
        }
        
        //! Get the elevation of the planewave.
        /** The method returns the elevation \f$\varphi\f$ that is always equal to \f$0\f$.
         @param x_axe The rotation around the x axe in radian (ignored).
         @param y_axe The rotation around the y axe in radian (ignored).
         @param z_axe The rotation around the z axe in radian (ignored).
         @return     The 0.
         */
        inline T getElevation(const T x_axe, const T y_axe, const T z_axe) const noexcept
        {
            return 0.;
        }
        
        //! Set the elevation of the planewave.
        /** The method sets the elevation \f$\varphi\f$ of the planewave.
         @param elevation The elevation \f$\varphi\f$ (ignored).
         @return     The elevation.
         */
        inline void setElevation(const T elevation) const noexcept
        {
            ;
        }
        
        //! Get the abscissa of the planewave.
        /** The method returns the abscissa \f$x\f$ between \f$-1\f$ and \f$1\f$ of the planewave. The result will consider a rotation around the axes x.
         @param x_axe The rotation around the x axe in radian (ignored).
         @param y_axe The rotation around the y axe in radian (ignored).
         @param z_axe The rotation around the z axe in radian.
         @return     The abscissa.
         */
        inline T getAbscissa(const T x_axe, const T y_axe, const T z_axe) const noexcept
        {
            return cos(z_axe + m_azimuth + HOA_PI2);
        }
        
        //! Get the ordinate of the planewave.
        /** The method returns the ordinate \f$y\f$ between \f$-1\f$ and \f$1\f$ of the planewave. The result will consider a rotation around the axes x.
         @param x_axe The rotation around the x axe in radian (ignored).
         @param y_axe The rotation around the y axe in radian (ignored).
         @param z_axe The rotation around the z axe in radian.
         @return     The ordinate.
         */
        inline T getOrdinate(const T x_axe, const T y_axe, const T z_axe) const noexcept
        {
            return sin(z_axe + m_azimuth + HOA_PI2);
        }
        
        //! Get the height of the planewave.
        /** The method returns the height \f$z\f$ that is always equal to \f$0\f$.
         @param x_axe The rotation around the x axe in radian (ignored).
         @param y_axe The rotation around the y axe in radian (ignored).
         @param z_axe The rotation around the z axe in radian (ignored).
         @return     The 0.
         */
        inline T getHeight(const T x_axe, const T y_axe, const T z_axe) const noexcept
        {
            return 0.;
        }
        
        //! Get the name of the planewave.
        /** The method returns the name \f$planewave_{index}\f$ of the planewave.
         @return     The name.
         */
        inline string getName() const noexcept
        {
            return "Planewave " + to_string(getIndex()) + " " + to_string(getAzimuth(0., 0., 0.) / HOA_2PI * 360.) + "°";
        }
        
        //! Compare the azimuth of two planewaves.
        /** The method returns the true if the azimtuh of the planewave i is inferior to the azimtuh of the planewave j, otherwise false.
         @return     The true if the azimtuh of the planewave i is inferior to the azimtuh of the planewave j, otherwise false.
         */
        static bool sort_azimuth(Planewave const& i, Planewave const& j) noexcept
        {
            return i.m_azimuth < j.m_azimuth;
        }
    };
            
    template<typename T> class Planewave<Hoa3d, T>
    {
    private:
        ulong m_index;
        T     m_azimuth;
        T     m_elevation;
    public:
        
        //! The planewave constructor.
        /** The planewave constructor allocates and initializes the generale member values depending on an index and a polar coordinate with an azimuth \f$\theta\f$ and an elevation \f$\varphi\f$ in radian assuming that the radius \f$\rho\f$ is always equal to \f$1\f$.
         @param index       The index must be at least 1.
         @param azimuth     The azimuth \f$\theta\f$.
         @param elevation   The elevation \f$\varphi\f$.
         */
        Planewave(const ulong _index, const T _azimuth, const T _elevation) noexcept :
        m_index(_index),
        m_azimuth(_azimuth),
        m_elevation(_elevation)
        {
            ;
        }
        
        //! The planewave constructor.
        /** The planewave constructor allocates and initializes the generale member values depending on an index and a cartesian coordinate with an abscissa \f$x\f$, an ordinate \f$y\f$ and an height \f$z\f$ assuming that the point is normalized over the unit circle or unit the sphere.
         @param index       The index must be at least 1.
         @param abscissa    The abscissa \f$x\f$.
         @param ordinate    The ordinate \f$y\f$.
         @param height      The height \f$z\f$.
         */
        Planewave(const ulong _index, const T _abscissa, const T _ordinate, const T _height) noexcept :
        m_index(_index),
        m_azimuth(Math<T>::azimuth(_abscissa, _ordinate, _height)),
        m_elevation(Math<T>::elevation(_abscissa, _ordinate, _height))
        {
            ;
        }
    
        //! The planewave destructor.
        /** The planewave destructor free the memory.
         */
        ~Planewave() noexcept
        {
            ;
        }
        
        //! Get the index of the planewave.
        /** The method returns the index \f$i\f$ of the planewave.
         @return     The index.
         */
        inline ulong getIndex() const noexcept
        {
            return m_index;
        }
        
        //! Get the azimuth of the planewave.
        /** The method returns the azimuth \f$\theta\f$ between \f$0\f$ and \f$2\pi\f$ of the planewave. The result will consider a full rotation around the axes x, y, and z. The order the rotation is x, z then y.
         @param x_axe The rotation around the x axe in radian.
         @param y_axe The rotation around the y axe in radian.
         @param z_axe The rotation around the z axe in radian.
         @return     The azimuth.
         */
        inline T getAzimuth(const T x_axe, const T y_axe, const T z_axe) const noexcept
        {
            const T x = getAbscissa(x_axe, y_axe, z_axe);
            const T y = getOrdinate(x_axe, y_axe, z_axe);
            if(x == 0. && y == 0.)
            {
                return 0.;
            }
            else
            {
                return Math<T>::wrap_twopi(atan2(y, x) - HOA_PI2);
            }
        }
        
        //! Set the azimuth of the planewave.
        /** The method sets the azimuth \f$\theta\f$ of the planewave.
         @param azimuth The azimuth \f$\theta\f$.
         @return     The azimuth.
         */
        inline void setAzimuth(const T azimuth) noexcept
        {
            m_azimuth = azimuth;
        }
        
        //! Get the elevation of the planewave.
        /** The method returns the elevation \f$\varphi\f$ between \f$-\pi\f$ and \f$\pi\f$ of the planewave. The result will consider a full rotation around the axes x, y, and z. The order the rotation is x, z then y.
         @param x_axe The rotation around the x axe in radian.
         @param y_axe The rotation around the y axe in radian.
         @param z_axe The rotation around the z axe in radian.
         @return     The elevation.
         */
        inline T getElevation(const T x_axe, const T y_axe, const T z_axe) const noexcept
        {
            const T z = getHeight(x_axe, y_axe, z_axe);
            if(z == 0.)
            {
                return 0.;
            }
            else
            {
                const T x = getAbscissa(x_axe, y_axe, z_axe);
                const T y = getOrdinate(x_axe, y_axe, z_axe);
                return Math<T>::wrap_pi(asin(z / sqrt(x*x + y*y + z*z)));
            }
        }
        
        //! Set the elevation of the planewave.
        /** The method sets the elevation \f$\varphi\f$ of the planewave.
         @param elevation The elevation \f$\varphi\f$.
         @return     The elevation.
         */
        inline void setElevation(const T elevation) noexcept
        {
            m_elevation = elevation;
        }
        
        //! Get the abscissa of the planewave.
        /** The method returns the abscissa \f$x\f$ between \f$-1\f$ and \f$1\f$ of the planewave. The result will consider a full rotation around the axes x, y, and z. The order the rotation is x, z then y.
         @param x_axe The rotation around the x axe in radian.
         @param y_axe The rotation around the y axe in radian.
         @param z_axe The rotation around the z axe in radian.
         @return     The abscissa.
         */
        inline T getAbscissa(const T x_axe, const T y_axe, const T z_axe) const noexcept
        {
            T x = cos(m_azimuth + HOA_PI2) * cos(m_elevation);
            T y = sin(m_azimuth + HOA_PI2) * cos(m_elevation);
            T z = sin(m_elevation);
            
            const T cosAngle = cos(x_axe);
            const T sinAngle = sin(x_axe);
            T ry = y * cosAngle - z * sinAngle;
            z = y * sinAngle + z * cosAngle;
            y = ry;
            
            x = x * cos(z_axe) - y * sin(z_axe);
            
            return x * cos(y_axe) - z * sin(y_axe);
        }
        
        //! Get the ordinate of the planewave.
        /** The method returns the ordinate \f$y\f$ between \f$-1\f$ and \f$1\f$ of the planewave. The result will consider a full rotation around the axes x, y, and z. The order the rotation is x, z then y.
         @param x_axe The rotation around the x axe in radian.
         @param y_axe The rotation around the y axe in radian.
         @param z_axe The rotation around the z axe in radian.
         @return     The ordinate.
         */
        inline T getOrdinate(const T x_axe, const T y_axe, const T z_axe) const noexcept
        {
            return cos(m_azimuth + HOA_PI2) * cos(m_elevation) * sin(z_axe) + ((sin(m_azimuth + HOA_PI2) * cos(m_elevation)) * cos(x_axe) - sin(m_elevation) * sin(x_axe)) * cos(z_axe);
        }
        
        //! Get the height of the planewave.
        /** The method returns the height \f$z\f$ between \f$-1\f$ and \f$1\f$ of the planewave. The result will consider a full rotation around the axes x, y, and z. The order the rotation is x, z then y.
         @param x_axe The rotation around the x axe in radian.
         @param y_axe The rotation around the y axe in radian.
         @param z_axe The rotation around the z axe in radian.
         @return     The height.
         */
        inline T getHeight(const T x_axe, const T y_axe, const T z_axe) const noexcept
        {
            T x = cos(m_azimuth + HOA_PI2) * cos(m_elevation);
            T y = sin(m_azimuth + HOA_PI2) * cos(m_elevation);
            T z = sin(m_elevation);
            
            T cosAngle = cos(x_axe);
            T sinAngle = sin(x_axe);
            T ry = y * cosAngle - z * sinAngle;
            T rz = y * sinAngle + z * cosAngle;
            y = ry;
            z = rz;
            
            cosAngle = cos(z_axe);
            sinAngle = sin(z_axe);
            x = x * cosAngle - y * sinAngle;
            
            return x * sin(y_axe) + z * cos(y_axe);
        }
        
        //! Get the name of the planewave.
        /** The method returns the name \f$planewave_{index}\f$ of the planewave.
         @return     The name.
         */
        inline string getName() const noexcept
        {
            return "Planewave " + to_string(getIndex()) + " " + to_string(getAzimuth(0., 0., 0.) / HOA_2PI * 360.) + "° " + to_string(getElevation(0., 0., 0.) / HOA_2PI * 360.) + "°";
        }
        
        //! Compare the azimuth of two planewaves.
        /** The method returns the true if the azimtuh of the planewave i is inferior to the azimtuh of the planewave j, otherwise false.
         @return     The true if the azimtuh of the planewave i is inferior to the azimtuh of the planewave j, otherwise false.
         */
        static bool sort_azimuth(Planewave const& i, Planewave const& j) noexcept
        {
            return i.m_azimuth < j.m_azimuth;
        }
    };

#endif
}

#endif


