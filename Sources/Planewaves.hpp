/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_PLANEWAVES_LIGHT
#define DEF_HOA_PLANEWAVES_LIGHT

#include "HoaMath.hpp"

namespace hoa
{
    template<Dimension D, typename T> class Planewave
    {
        ;
    };
    
    template<typename T> class Planewave<Hoa2d, T>
    {
    private:
        ulong m_index;
        T     m_azimuth;
    public:
        
        Planewave() noexcept :
        m_index(0),
        m_azimuth(0)
        {
            ;
        }
        
        Planewave(const ulong _index, const T _azimuth) noexcept :
        m_index(_index),
        m_azimuth(_azimuth)
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
        
        inline T getAbscissa(const T offset = 0.) const noexcept
        {
            return std::cos(m_azimuth + offset + HOA_PI2);
        }
        
        inline T getOrdinate(const T offset = 0.) const noexcept
        {
            return std::sin(m_azimuth + offset + HOA_PI2);
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
            vector<Planewave<Hoa2d, T>> m_planewaves;
            T                           m_offset;
        public:
            
            //! The planewaves constructor.
            /** The planewaves constructor allocates and initializes the general member values depending on a number of planewaves. The number of planewaves must be a least 1.
             @param     numberOfPlanewaves	The number of planewaves.
             */
            Processor(const ulong numberOfPlanewaves) noexcept :
            m_number_of_planewaves(numberOfPlanewaves),
            m_offset(0.)
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
            inline void setPlanewavesOffset(T offset) noexcept
            {
                m_offset = wrap_twopi(offset);
            }
            
            //! Get the offset of the planewaves.
            /**	Retreive the azimuth offset of the planewaves in radian.
             @return    The offset of the planewaves.
             */
            inline T getPlanewavesOffset() const noexcept
            {
                return m_offset;
            }
            
            //! Retrieve the azimuth of a planewaves.
            /** Retrieve the azimuth of a planewaves. The azimuth of the planewaves is in radian, 0 radian is at the front of the soundfield and π is at the back of the sound field. The maximum index must be the number of planewaves - 1.
             
             @param      index   The index of the planewaves.
             @return     The azimuth of the planewaves.
             */
            inline void setPlanewaveAzimuth(const ulong index, const T azimuth) noexcept
            {
                m_planewaves[index].setAzimuth(azimuth);
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
            
            //! Retrieve the abscissa of a planewaves.
            /** Retrieve the abscissa of a planewaves. The abscissa is between -1 and 1, -1 is the left of the soundfield, 0 is the center of the soundfield and 1 is the right of the soundfield. The maximum index must be the number of planewaves - 1.
             @param     index    The index of the planewaves.
             @return    The abscissa of the planewaves.
             */
            inline T getPlanewaveAbscissa(const ulong index) const noexcept
            {
                return m_planewaves[index].getAbscissa(m_offset);
            }
            
            //! Retrieve the ordinate of a planewaves.
            /** Retrieve the ordinate of a planewaves. The ordinate is between -1 and 1, -1 is the back of the soundfield, 0 is the center of the soundfield and 1 is the front of the soundfield. The maximum index must be the number of planewaves - 1.
             @param     index	The index of the planewaves.
             @return    The ordinate of the planewaves.
             */
            inline T getPlanewaveOrdinate(const ulong index) const noexcept
            {
                return m_planewaves[index].getOrdinate(m_offset);
            }
            
            //! Retrieve a name for a planewaves.
            /** Retrieve a name for a planewaves in a std::string format that will be "Planewave index azimuth (in degrees)".
             @param     index	The index of a planewaves.
             @return    The method returns a name for the planewaves.
             */
            inline std::string getPlanewaveName(const ulong index) const noexcept
            {
                return m_planewaves[index].getName();
            }
        };
    };
}

#endif


