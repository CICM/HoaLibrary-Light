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
    template<Dimension D, typename T> class Planewave;
    
    template<typename T> class Planewave<Hoa2d, T>
    {
    private:
        ulong m_index;
        T     m_azimuth;
 
    public:
        
        Planewave(const ulong _index, const T _azimuth) noexcept :
        m_index(_index),
        m_azimuth(_azimuth)
        {
            
        }
        
        ~Planewave()
        {
            ;
        }
        
        inline ulong getIndex() const noexcept
        {
            return m_index;
        }
        
        inline T getAzimuth() const noexcept
        {
            return m_azimuth;
        }
        
        inline void setAzimuth(const T azimuth) noexcept
        {
            m_azimuth = azimuth;
        }
        
        inline T getAbscissa(const T x_axe = 0.) const noexcept
        {
            return cos(m_azimuth + x_axe + HOA_PI2);
        }
        
        inline T getOrdinate(const T x_axe = 0.) const noexcept
        {
            return sin(m_azimuth + x_axe + HOA_PI2);
        }
        
        inline string getName() const noexcept
        {
            return "Planewave " + to_string(getIndex()) + " " + to_string(getAzimuth() / HOA_2PI * 360.) + "°";
        }
        
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
            vector<Planewave<Hoa2d, T> > m_planewaves;
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
            ~Processor()
            {
                m_planewaves.clear();
            }
            
            //! Retrieve the decomposition order.
            /** Retrieve the decomposition order.
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
             
             @param      index   The index of the planewaves.
             @return     The azimuth of the planewaves.
             */
            inline T getPlanewaveAzimuth(const ulong index) const noexcept
            {
                return m_planewaves[index].getAzimuth();
            }
            
            //! Retrieve the azimuth of a planewaves (including rotation).
            /** Retrieve the azimuth of a planewaves. The azimuth of the planewaves is in radian, 0 radian is at the front of the soundfield and π is at the back of the sound field. The maximum index must be the number of planewaves - 1.
             
             @param      index   The index of the planewaves.
             @return     The azimuth of the planewaves.
             */
            inline T getPlanewaveAzimuthRotated(const ulong index) const noexcept
            {
                return Math<T>::wrap_twopi(m_planewaves[index].getAzimuth() + m_rotation_z);
            }
            
            //! Retrieve the abscissa of a planewaves.
            /** Retrieve the abscissa of a planewaves. The abscissa is between -1 and 1, -1 is the left of the soundfield, 0 is the center of the soundfield and 1 is the right of the soundfield. The maximum index must be the number of planewaves - 1.
             @param     index    The index of the planewaves.
             @return    The abscissa of the planewaves.
             */
            inline T getPlanewaveAbscissa(const ulong index) const noexcept
            {
                return m_planewaves[index].getAbscissa(m_rotation_z);
            }
            
            //! Retrieve the ordinate of a planewaves.
            /** Retrieve the ordinate of a planewaves. The ordinate is between -1 and 1, -1 is the back of the soundfield, 0 is the center of the soundfield and 1 is the front of the soundfield. The maximum index must be the number of planewaves - 1.
             @param     index	The index of the planewaves.
             @return    The ordinate of the planewaves.
             */
            inline T getPlanewaveOrdinate(const ulong index) const noexcept
            {
                return m_planewaves[index].getOrdinate(m_rotation_z);
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
        
        Planewave(const ulong _index, const T _azimuth, const T _elevation) noexcept :
        m_index(_index),
        m_azimuth(_azimuth),
        m_elevation(_elevation)
        {
            ;
        }
        
        ~Planewave()
        {
            ;
        }
        
        inline ulong getIndex() const noexcept
        {
            return m_index;
        }
        
        inline T getAzimuth() const noexcept
        {
            return m_azimuth;
        }
        
        inline void setAzimuth(const T azimuth) noexcept
        {
            m_azimuth = azimuth;
        }
        
        inline T getElevation() const noexcept
        {
            return m_elevation;
        }
        
        inline void setElevation(const T elevation) noexcept
        {
            m_elevation = elevation;
        }
        
        inline T getAbscissa(const T x_axe = 0., const T y_axe = 0., const T z_axe = 0.) const noexcept
        {
            T x = cos(m_azimuth + HOA_PI2) * cos(m_elevation);
            T y = sin(m_azimuth + HOA_PI2) * cos(m_elevation);
            T z = sin(m_elevation);
            
            // Rotation around x
            T cosAngle = cos(x_axe);
            T sinAngle = sin(x_axe);
            T ry = y * cosAngle - z * sinAngle;
            z = y * sinAngle + z * cosAngle;
            y = ry;
            
            // Rotation around z
            cosAngle = cos(z_axe);
            sinAngle = sin(z_axe);
            x = x * cosAngle - y * sinAngle;
            
            // Rotation around y
            cosAngle = cos(y_axe);
            sinAngle = sin(y_axe);
            return x * cosAngle - z * sinAngle;
        }
        
        inline T getOrdinate(const T x_axe = 0., const T y_axe = 0., const T z_axe = 0.) const noexcept
        {
            T x = cos(m_azimuth + HOA_PI2) * cos(m_elevation);
            T y = sin(m_azimuth + HOA_PI2) * cos(m_elevation);
            T z = sin(m_elevation);
            
            // Rotation around x
            T cosAngle = cos(x_axe);
            T sinAngle = sin(x_axe);
            y = y * cosAngle - z * sinAngle;
            
            // Rotation around z
            cosAngle = cos(z_axe);
            sinAngle = sin(z_axe);
            return x * sinAngle + y * cosAngle;
        }
        
        inline T getHeight(const T x_axe = 0., const T y_axe = 0., const T z_axe = 0.) const noexcept
        {
            T x = cos(m_azimuth + HOA_PI2) * cos(m_elevation);
            T y = sin(m_azimuth + HOA_PI2) * cos(m_elevation);
            T z = sin(m_elevation);
            
            // Rotation around x
            T cosAngle = cos(x_axe);
            T sinAngle = sin(x_axe);
            T ry = y * cosAngle - z * sinAngle;
            z = y * sinAngle + z * cosAngle;
            y = ry;
            
            // Rotation around z
            cosAngle = cos(z_axe);
            sinAngle = sin(z_axe);
            x = x * cosAngle - y * sinAngle;;
            
            // Rotation around y
            cosAngle = cos(y_axe);
            sinAngle = sin(y_axe);
            return x * sinAngle + z * cosAngle;
        }
        
        inline string getName() const noexcept
        {
            return "Planewave " + to_string(getIndex()) + " " + to_string(getAzimuth() / HOA_2PI * 360.) + "°";
        }
        
        bool operator<(Planewave const& j) const noexcept
        {
            return this->m_azimuth < j.m_azimuth;
        }
        
        //! The planewaves class.
        /**
         The planewaves classes, that process on a set of planewaves inherit from this class. It store basic informations like the number of planewaves, the coordinates and the names of the planewaves.
         */
        class Processor
        {
        private:
            const ulong                 m_number_of_planewaves;
            vector<Planewave<Hoa3d, T> > m_planewaves;
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
                for(ulong i = 0; i < m_number_of_planewaves; i++)
                {
                    m_planewaves.push_back(Planewave<Hoa3d, T>(i+1, (T)i / (m_number_of_planewaves) * HOA_2PI, 0.));
                }
            }
            
            //! The ambisonic destructor.
            /** The ambisonic destructor.
             */
            ~Processor()
            {
                m_planewaves.clear();
            }
            
            //! Retrieve the decomposition order.
            /** Retrieve the decomposition order.
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
             @return     The azimuth of the planewaves.
             */
            inline T getPlanewaveAzimuth(const ulong index) const noexcept
            {
                return m_planewaves[index].getAzimuth();
            }
            
            //! Retrieve the azimuth of a planewaves.
            /** Retrieve the azimuth of a planewaves. The azimuth of the planewaves is in radian, 0 radian is at the front of the soundfield and π is at the back of the sound field. The maximum index must be the number of planewaves - 1.
             
             @param      index   The index of the planewaves.
             @return     The azimuth of the planewaves.
             */
            inline T getPlanewaveAzimuthRotated(const ulong index) const noexcept
            {
                return Math<T>::azimuth(getPlanewaveAbscissa(index), getPlanewaveOrdinate(index), getPlanewaveHeight(index));
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
             
             @param      index   The index of the planewaves.
             @return     The azimuth of the planewaves.
             */
            inline T getPlanewaveElevation(const ulong index) const noexcept
            {
                return m_planewaves[index].getElevation();
            }
            
            //! Retrieve the azimuth of a planewaves.
            /** Retrieve the azimuth of a planewaves. The azimuth of the planewaves is in radian, 0 radian is at the front of the soundfield and π is at the back of the sound field. The maximum index must be the number of planewaves - 1.
             
             @param      index   The index of the planewaves.
             @return     The azimuth of the planewaves.
             */
            inline T getPlanewaveElevationRotated(const ulong index) const noexcept
            {
                return Math<T>::elevation(getPlanewaveAbscissa(index), getPlanewaveOrdinate(index), getPlanewaveHeight(index));
            }
            
            //! Retrieve the abscissa of a planewaves.
            /** Retrieve the abscissa of a planewaves. The abscissa is between -1 and 1, -1 is the left of the soundfield, 0 is the center of the soundfield and 1 is the right of the soundfield. The maximum index must be the number of planewaves - 1.
             @param     index    The index of the planewaves.
             @return    The abscissa of the planewaves.
             */
            inline T getPlanewaveAbscissa(const ulong index) const noexcept
            {
                return m_planewaves[index].getAbscissa(m_rotation_x, m_rotation_y, m_rotation_z);
            }
            
            //! Retrieve the ordinate of a planewaves.
            /** Retrieve the ordinate of a planewaves. The ordinate is between -1 and 1, -1 is the back of the soundfield, 0 is the center of the soundfield and 1 is the front of the soundfield. The maximum index must be the number of planewaves - 1.
             @param     index	The index of the planewaves.
             @return    The ordinate of the planewaves.
             */
            inline T getPlanewaveOrdinate(const ulong index) const noexcept
            {
                return m_planewaves[index].getOrdinate(m_rotation_x, m_rotation_y, m_rotation_z);
            }
            
            //! Retrieve the height of a planewaves.
            /** Retrieve the ordinate of a planewaves. The height is between -1 and 1, -1 is the back of the soundfield, 0 is the center of the soundfield and 1 is the front of the soundfield. The maximum index must be the number of planewaves - 1.
             @param     index	The index of the planewaves.
             @return    The height of the planewaves.
             */
            inline T getPlanewaveHeight(const ulong index) const noexcept
            {
                return m_planewaves[index].getHeight(m_rotation_x, m_rotation_y, m_rotation_z);
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
}

#endif


