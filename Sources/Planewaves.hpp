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
        
        //! The planewaves processor.
        /**
         The planewaves processor classes, that process on a set of planewaves inherit from this class. It store basic informations like the number of planewaves, the coordinates and the names of the planewaves.
         */
        class Processor
        {
            ;
        };
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
        Planewave(const ulong _index, const T _azimuth, const T _elevation = 0.) noexcept :
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
        Planewave(const ulong _index, const T _abscissa, const T _ordinate, const T _height = 0.) noexcept :
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
        
        //! The planewaves processor.
        /**
         The planewaves processor classes, that process on a set of planewaves inherit from this class. It store basic informations like the number of planewaves, the coordinates and the names of the planewaves.
         */
        class Processor
        {
        private:
            const ulong                 m_number_of_planewaves;
            vector<Planewave<Hoa2d, T> >m_planewaves;
            T                           m_rotation_z;
        public:
            
            //! The planewaves constructor.
            /** The planewaves constructor allocates and initializes the general member values depending on a number of planewaves. The number of planewaves must be a least 1.
             @param     numberOfPlanewaves	The number of planewaves.
             */
            Processor(const ulong numberOfPlanewaves) noexcept :
            m_number_of_planewaves(numberOfPlanewaves),
            m_rotation_z(0.)
            {
                for(ulong i = 0; i < m_number_of_planewaves; i++)
                {
                    m_planewaves.push_back(Planewave<Hoa2d, T>(i+1, (T)i / (m_number_of_planewaves) * HOA_2PI));
                }
            }
            
            //! The ambisonic destructor.
            /** The ambisonic destructor.
             */
            ~Processor() noexcept
            {
                m_planewaves.clear();
            }
            
            //! Retrieve the order of decomposition.
            /** Retrieve the order of decomposition.
             @return The order.
             */
            inline ulong getNumberOfPlanewaves() const noexcept
            {
                return m_number_of_planewaves;
            }
            
            //! Set the offset of the planewaves.
            /**	Set the azimuth offset of the planewaves in radian.
             @param     offset		An azimuth value.
             */
            inline void setPlanewavesRotation(const T z_axe) noexcept
            {
                m_rotation_z = Math<T>::wrap_twopi(z_axe);
            }
            
            //! Get the offset of the planewaves.
            /**	Retreive the azimuth offset of the planewaves in radian.
             @return    The offset of the planewaves.
             */
            inline T getPlanewavesRotation() const noexcept
            {
                return m_rotation_z;
            }
            
            //! Retrieve the index of a planewaves.
            /** Retrieve the index of a planewaves.
             @param      index   The index of the planewaves.
             @return     The index of the planewaves.
             */
            inline ulong getPlanewaveIndex(const ulong index) noexcept
            {
                return m_planewaves[index].getIndex();
            }
            
            //! Sets the azimuth of a planewaves.
            /** The azimuth of the planewaves is in radian, 0 radian is at the front of the soundfield and π is at the back of the sound field. The maximum index must be the number of planewaves - 1.
             
             @param      index   The index of the planewaves.
             @param      azimuth The new azimuth value.
             */
            inline void setPlanewaveAzimuth(const ulong index, const T azimuth) noexcept
            {
                m_planewaves[index].setAzimuth(Math<T>::wrap_twopi(azimuth));
            }
            
            //! Retrieve the azimuth of a planewaves.
            /** The azimuth of the planewaves is in radian, 0 radian is at the front of the soundfield and π is at the back of the sound field. The maximum index must be the number of planewaves - 1.
             
             @param      index      The index of the planewaves.
             @param      rotation   False if you don't want to consider the rotation, otherwise true (default).
             @return     The azimuth of the planewaves.
             */
            inline T getPlanewaveAzimuth(const ulong index, const bool rotation = true) const noexcept
            {
                return m_planewaves[index].getAzimuth(0., 0., rotation ? m_rotation_z : 0.);
            }
            
            //! Retrieve the abscissa of a planewaves.
            /** Retrieve the abscissa of a planewaves. The abscissa is between -1 and 1, -1 is the left of the soundfield, 0 is the center of the soundfield and 1 is the right of the soundfield. The maximum index must be the number of planewaves - 1.
             @param     index    The index of the planewaves.
             @param     rotation   False if you don't want to consider the rotation, otherwise true (default).
             @return    The abscissa of the planewaves.
             */
            inline T getPlanewaveAbscissa(const ulong index, const bool rotation = true) const noexcept
            {
                return m_planewaves[index].getAbscissa(0., 0., rotation ? m_rotation_z : 0.);
            }
            
            //! Retrieve the ordinate of a planewaves.
            /** Retrieve the ordinate of a planewaves. The ordinate is between -1 and 1, -1 is the back of the soundfield, 0 is the center of the soundfield and 1 is the front of the soundfield. The maximum index must be the number of planewaves - 1.
             @param     index	The index of the planewaves.
             @param     rotation   False if you don't want to consider the rotation, otherwise true (default).
             @return    The ordinate of the planewaves.
             */
            inline T getPlanewaveOrdinate(const ulong index, const bool rotation = true) const noexcept
            {
                return m_planewaves[index].getOrdinate(0., 0., rotation ? m_rotation_z : 0.);
            }
            
            //! Retrieve a name for a planewaves.
            /** Retrieve a name for a planewaves in a string format that will be "Planewave index azimuth (in degrees)".
             @param     index	The index of a planewaves.
             @return    The method returns a name for the planewaves.
             */
            inline string getPlanewaveName(const ulong index) const noexcept
            {
                return m_planewaves[index].getName();
            }
        };
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
            return "Planewave " + to_string(getIndex()) + " " + to_string(getAzimuth(0., 0., 0.) / HOA_2PI * 360.) + "°" " " + to_string(getElevation(0., 0., 0.) / HOA_2PI * 360.) + "°";
        }
        
        //! Compare the azimuth of two planewaves.
        /** The method returns the true if the azimtuh of the planewave i is inferior to the azimtuh of the planewave j, otherwise false.
         @return     The true if the azimtuh of the planewave i is inferior to the azimtuh of the planewave j, otherwise false.
         */
        static bool sort_azimuth(Planewave const& i, Planewave const& j) noexcept
        {
            return i.m_azimuth < j.m_azimuth;
        }
        
        
        //! The planewaves class.
        /**
         The planewaves classes, that process on a set of planewaves inherit from this class. It store basic informations like the number of planewaves, the coordinates and the names of the planewaves.
         */
        class Processor
        {
        private:
            const ulong                 m_number_of_planewaves;
            vector<Planewave<Hoa3d, T> >m_planewaves;
            T                           m_rotation_z;
            T                           m_rotation_y;
            T                           m_rotation_x;
        public:
            
            //! The planewaves constructor.
            /** The planewaves constructor allocates and initializes the general member values depending on a number of planewaves. The number of planewaves must be a least 1.
             @param     numberOfPlanewaves	The number of planewaves.
             */
            Processor(const ulong numberOfPlanewaves) noexcept :
            m_number_of_planewaves(numberOfPlanewaves),
            m_rotation_z(0.),
            m_rotation_y(0.),
            m_rotation_x(0.)
            {
                if(m_number_of_planewaves == 4)
                {
                    const double oh = -(sqrt(2. / 3.) / sqrt(3. / 8.) - 1.);
                    const double hc = sqrt(1. - oh * oh);
                    const double el = asin(oh / sqrt(hc*hc + oh*oh));
                    m_planewaves.push_back(Planewave<Hoa3d, T>(1,  0., HOA_PI2));
                    m_planewaves.push_back(Planewave<Hoa3d, T>(2, 0., el));
                    m_planewaves.push_back(Planewave<Hoa3d, T>(3, HOA_2PI / 3., el));
                    m_planewaves.push_back(Planewave<Hoa3d, T>(4, 2. * HOA_2PI / 3., el));
                }
                else if(m_number_of_planewaves == 6)
                {
                    m_planewaves.push_back(Planewave<Hoa3d, T>(1, 0., HOA_PI2));
                    m_planewaves.push_back(Planewave<Hoa3d, T>(2, 0., 0.));
                    m_planewaves.push_back(Planewave<Hoa3d, T>(3, HOA_PI2, 0.));
                    m_planewaves.push_back(Planewave<Hoa3d, T>(4, 2. * HOA_PI2, 0.));
                    m_planewaves.push_back(Planewave<Hoa3d, T>(5, 3. * HOA_PI2, 0.));
                    m_planewaves.push_back(Planewave<Hoa3d, T>(6, 0., -HOA_PI2));
                }
                else if(m_number_of_planewaves == 8)
                {
                    m_planewaves.push_back(Planewave<Hoa3d, T>(1, 1., 1., 1.));
                    m_planewaves.push_back(Planewave<Hoa3d, T>(2, -1., 1., 1.));
                    m_planewaves.push_back(Planewave<Hoa3d, T>(3, -1., -1., 1.));
                    m_planewaves.push_back(Planewave<Hoa3d, T>(4, 1., -1., 1.));
                    
                    m_planewaves.push_back(Planewave<Hoa3d, T>(5, 1., 1., -1.));
                    m_planewaves.push_back(Planewave<Hoa3d, T>(6, -1., 1., -1.));
                    m_planewaves.push_back(Planewave<Hoa3d, T>(7, -1., -1., -1.));
                    m_planewaves.push_back(Planewave<Hoa3d, T>(8, 1., -1., -1.));
                }
                else if(m_number_of_planewaves == 12)
                {
                    m_planewaves.push_back(Planewave<Hoa3d, T>(1, 0., HOA_PI2));
                    for(ulong i = 1; i < 6; i++)
                    {
                        m_planewaves.push_back(Planewave<Hoa3d, T>(i*2, T(i - 1.) / 5. * HOA_2PI, atan(0.5)));
                        m_planewaves.push_back(Planewave<Hoa3d, T>(i*2+1, T(i - 1.) / 5. * HOA_2PI - HOA_PI / 5., -atan(0.5)));
                    }
                    m_planewaves.push_back(Planewave<Hoa3d, T>(12, 0., -HOA_PI2));
                }
                else if(m_number_of_planewaves == 20)
                {
                    const T phi = (sqrt(5.) - 1.) / 2.;
                    const T R = 1. / sqrt(3.);
                    const T a = R;
                    const T b = R / phi;
                    const T c = R * phi;
                    ulong index = 1;
                    for(long i = -1; i < 2; i += 2)
                    {
                        for(long j = -1; j < 2; j += 2)
                        {
                            m_planewaves.push_back(Planewave<Hoa3d, T>(index++, 0., i * c * R, -j * b * R));
                            m_planewaves.push_back(Planewave<Hoa3d, T>(index++, i * c * R, j * b * R, 0.));
                            m_planewaves.push_back(Planewave<Hoa3d, T>(index++, i * b * R, 0., -j * c * R));
                            for(long k = -1; k < 2; k += 2)
                            {
                                m_planewaves.push_back(Planewave<Hoa3d, T>(index++, i * a * R, j * a * R, -k * a * R));
                            }
                        }
                    }
                }
                else
                {
                    if(m_number_of_planewaves % 2)
                    {
                        m_planewaves.push_back(Planewave<Hoa3d, T>(1, 0., HOA_PI2));
                    }
                    const T     phi     = (sqrt(5.) - 1.) / 4.;
                    const ulong limit   = (m_number_of_planewaves - (m_number_of_planewaves % 2)) / 2;
                    const T     offset  = 1. / T(limit) * HOA_PI;
                    for(ulong i = 0; i < m_number_of_planewaves - (m_number_of_planewaves % 2); i++)
                    {
                        if(i < limit)
                        {
                            m_planewaves.push_back(Planewave<Hoa3d, T>(i+1, T(i) / T(limit) * HOA_2PI, phi));
                        }
                        else
                        {
                            m_planewaves.push_back(Planewave<Hoa3d, T>(i+1, T(i - (limit)) / T(limit) * HOA_2PI + offset, phi * 2));
                        }
                    }
                }
            }
            
            //! The ambisonic destructor.
            /** The ambisonic destructor.
             */
            ~Processor() noexcept
            {
                m_planewaves.clear();
            }
            
            //! Retrieve the order of decomposition.
            /** Retrieve the order of decomposition.
             @return The order.
             */
            inline ulong getNumberOfPlanewaves() const noexcept
            {
                return m_number_of_planewaves;
            }
            
            //! Set the offset of the planewaves.
            /**	Set the azimuth offset of the planewaves in radian.
             @param     offset		An azimuth value.
             */
            inline void setPlanewavesRotation(const T x_axe, const T y_axe, const T z_axe) noexcept
            {
                m_rotation_x = Math<T>::wrap_twopi(x_axe);
                m_rotation_y = Math<T>::wrap_twopi(y_axe);
                m_rotation_z = Math<T>::wrap_twopi(z_axe);
            }
            
            //! Get the offset of the planewaves.
            /**	Retreive the azimuth offset of the planewaves in radian.
             @return    The offset of the planewaves.
             */
            inline T getPlanewavesRotationX() const noexcept
            {
                return m_rotation_x;
            }
            
            //! Get the offset of the planewaves.
            /**	Retreive the azimuth offset of the planewaves in radian.
             @return    The offset of the planewaves.
             */
            inline T getPlanewavesRotationY() const noexcept
            {
                return m_rotation_y;
            }
            
            //! Get the offset of the planewaves.
            /**	Retreive the azimuth offset of the planewaves in radian.
             @return    The offset of the planewaves.
             */
            inline T getPlanewavesRotationZ() const noexcept
            {
                return m_rotation_z;
            }
            
            //! Retrieve the index of a planewaves.
            /** Retrieve the index of a planewaves.
             @param      index   The index of the planewaves.
             @return     The index of the planewaves.
             */
            inline ulong getPlanewaveIndex(const ulong index) noexcept
            {
                return m_planewaves[index].getIndex();
            }
            
            //! Retrieve the azimuth of a planewaves.
            /** Retrieve the azimuth of a planewaves. The azimuth of the planewaves is in radian, 0 radian is at the front of the soundfield and π is at the back of the sound field. The maximum index must be the number of planewaves - 1.
             
             @param      index   The index of the planewaves.
             @return     The azimuth of the planewaves.
             */
            inline void setPlanewaveAzimuth(const ulong index, const T azimuth) noexcept
            {
                m_planewaves[index].setAzimuth(Math<T>::wrap_twopi(azimuth));
            }
            
            //! Retrieve the azimuth of a planewaves.
            /** Retrieve the azimuth of a planewaves. The azimuth of the planewaves is in radian, 0 radian is at the front of the soundfield and π is at the back of the sound field. The maximum index must be the number of planewaves - 1.
             
             @param      index   The index of the planewaves.
             @param      rotation   False if you don't want to consider the rotation, otherwise true (default).
             @return     The azimuth of the planewaves.
             */
            inline T getPlanewaveAzimuth(const ulong index, const bool rotation = true) const noexcept
            {
                return m_planewaves[index].getAzimuth(rotation ? m_rotation_x : T(0.), rotation ? m_rotation_y : T(0.), rotation ? m_rotation_z : T(0.));
            }
            
            //! Retrieve the azimuth of a planewaves.
            /** Retrieve the azimuth of a planewaves. The azimuth of the planewaves is in radian, 0 radian is at the front of the soundfield and π is at the back of the sound field. The maximum index must be the number of planewaves - 1.
             
             @param      index   The index of the planewaves.
             @return     The azimuth of the planewaves.
             */
            inline void setPlanewaveElevation(const ulong index, const T azimuth) noexcept
            {
                m_planewaves[index].setElevation(Math<T>::wrap_pi(azimuth));
            }
            
            //! Retrieve the azimuth of a planewaves.
            /** Retrieve the azimuth of a planewaves. The azimuth of the planewaves is in radian, 0 radian is at the front of the soundfield and π is at the back of the sound field. The maximum index must be the number of planewaves - 1.
             
             @param      index      The index of the planewaves.
             @param      rotation   False if you don't want to consider the rotation, otherwise true (default).
             @return     The azimuth of the planewaves.
             */
            inline T getPlanewaveElevation(const ulong index, const bool rotation = true) const noexcept
            {
                return m_planewaves[index].getElevation(rotation ? m_rotation_x : T(0.), rotation ? m_rotation_y : T(0.), rotation ? m_rotation_z : T(0.));
            }
            
            //! Retrieve the abscissa of a planewaves.
            /** Retrieve the abscissa of a planewaves. The abscissa is between -1 and 1, -1 is the left of the soundfield, 0 is the center of the soundfield and 1 is the right of the soundfield. The maximum index must be the number of planewaves - 1.
             @param     index    The index of the planewaves.
             @param      rotation   False if you don't want to consider the rotation, otherwise true (default).
             @return    The abscissa of the planewaves.
             */
            inline T getPlanewaveAbscissa(const ulong index, const bool rotation = true) const noexcept
            {
                return m_planewaves[index].getAbscissa(rotation ? m_rotation_x : T(0.), rotation ? m_rotation_y : T(0.), rotation ? m_rotation_z : T(0.));
            }
            
            //! Retrieve the ordinate of a planewaves.
            /** Retrieve the ordinate of a planewaves. The ordinate is between -1 and 1, -1 is the back of the soundfield, 0 is the center of the soundfield and 1 is the front of the soundfield. The maximum index must be the number of planewaves - 1.
             @param     index	The index of the planewaves.
             @param      rotation   False if you don't want to consider the rotation, otherwise true (default).
             @return    The ordinate of the planewaves.
             */
            inline T getPlanewaveOrdinate(const ulong index, const bool rotation = true) const noexcept
            {
                return m_planewaves[index].getOrdinate(rotation ? m_rotation_x : T(0.), rotation ? m_rotation_y : T(0.), rotation ? m_rotation_z : T(0.));
            }
            
            //! Retrieve the height of a planewaves.
            /** Retrieve the ordinate of a planewaves. The height is between -1 and 1, -1 is the back of the soundfield, 0 is the center of the soundfield and 1 is the front of the soundfield. The maximum index must be the number of planewaves - 1.
             @param     index	The index of the planewaves.
             @param      rotation   False if you don't want to consider the rotation, otherwise true (default).
             @return    The height of the planewaves.
             */
            inline T getPlanewaveHeight(const ulong index, const bool rotation = true) const noexcept
            {
                return m_planewaves[index].getHeight(rotation ? m_rotation_x : T(0.), rotation ? m_rotation_y : T(0.), rotation ? m_rotation_z : T(0.));
            }
            
            //! Retrieve a name for a planewaves.
            /** Retrieve a name for a planewaves in a string format that will be "Planewave index azimuth (in degrees)".
             @param     index	The index of a planewaves.
             @return    The method returns a name for the planewaves.
             */
            inline string getPlanewaveName(const ulong index) const noexcept
            {
                return m_planewaves[index].getName();
            }
        };
    };

#endif
}

#endif


