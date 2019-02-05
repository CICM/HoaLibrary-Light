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

#include "Hoa_Harmonics.hpp"
#include "Hoa_Planewaves.hpp"
#include <cassert>

namespace hoa
{
    // ================================================================================ //
    // PROCESSOR //
    // ================================================================================ //
    
    //! @brief The processor.
    //! @details The processor the performs digital signal processing.
    template <Dimension D, typename T>
    class Processor
    {
    public:
        //! @brief Destructor.
        virtual ~Processor() = default;
        
        //! @brief The pure virtual method performs that performs the digital signal processing.
        //! @param inputs  The inputs array.
        //! @param outputs The outputs array.
        virtual void process(const T* inputs, T* outputs) noexcept
        {
            hoa_unused(inputs);
            hoa_unused(outputs);
            
            assert(0 && "This method must be overriden");
        }
    };
    
    // ================================================================================ //
    // PROCESSOR HARMONICS //
    // ================================================================================ //
    
    //! @brief The harmonic processor.
    //! @details The harmonic processor owns a set of harmonics depending on the order of decomposition.
    template <Dimension D, typename T>
    class ProcessorHarmonics
    : virtual public Processor<D, T>
    {
    public:

        //! @brief The harmonics constructor.
        //! @param order The order of decomposition, must be at least 1.
        ProcessorHarmonics(const size_t order) noexcept
        : m_order_of_decomposition(order)
        , m_harmonics(createVector(order))
        {}

        //! @brief The harmonics destructor.
        virtual ~ProcessorHarmonics() = default;

        //! @brief Returns the order of decomposition.
        inline size_t getDecompositionOrder() const noexcept { return m_order_of_decomposition; }

        //! @brief Returns the number of harmonics.
        inline size_t getNumberOfHarmonics() const noexcept { return m_harmonics.size(); }

        //! @brief Returns the degree of an harmonic.
        //! @param index The index of an harmonic.
        inline size_t getHarmonicDegree(const size_t index) const
        {
            return m_harmonics[index].getDegree();
        }

        //! @brief Returns the azimuthal order of an harmonic.
        //! @param index The index of an harmonic.
        inline long getHarmonicOrder(const size_t index) const
        {
            return m_harmonics[index].getOrder();
        }

        //! @brief Returns the index of an harmonic given the degree and the azimuthal order.
        //! @param degree The degree of the harmonic.
        //! @param order  The azimuthal order of the harmonic.
        inline size_t getHarmonicIndex(const size_t degree, const long order) const
        {
            return Harmonic<D, T>::getIndex(degree, order);
        }

        //! @brief Returns the name of an harmonic.
        //! @param index The index of an harmonic.
        inline std::string getHarmonicName(const size_t index) const
        {
            return m_harmonics[index].getName();
        }
        
        //! @brief Returns the normalization of an harmonic.
        //! @param index The index of an harmonic.
        inline T getHarmonicNormalization(const size_t index) const
        {
            return m_harmonics[index].getNormalization();
        }
        
        //! @brief Returns the semi-normalization of an harmonic.
        //! @param index The index of an harmonic.
        inline T getHarmonicSemiNormalization(const size_t index) const
        {
            return m_harmonics[index].getSemiNormalization();
        }
        
    private:
        
        static inline std::vector< Harmonic<D, T> > createVector(const size_t order)
        {
            std::vector<Harmonic<D, T>> harmonics {};
            for(size_t i = 0; i < Harmonic<D, T>::getNumberOfHarmonics(order); ++i)
            {
                harmonics.push_back(Harmonic<D, T>(i));
            }
            return harmonics;
        }
        
    private:
        
        const size_t m_order_of_decomposition = 0ul;
        const std::vector< Harmonic<D, T> > m_harmonics;
    };
    
    // ================================================================================ //
    // PROCESSOR PLANEWAVES //
    // ================================================================================ //

    //! @brief The planewaves processor.
    //! @details The planewaves processor owns a set of planewaves.
    template <Dimension D, typename T>
    class ProcessorPlanewaves
    : virtual public Processor<D, T>
    {
    public:
        
        //! @brief The planewaves constructor.
        //! @param nplws The number of planewaves.
        ProcessorPlanewaves(const size_t nplws)
        {
            if(D == Hoa2d)
            {
                generateCicle(nplws);
            }
            else
            {
                switch(nplws)
                {
                    case 4  : { generateTetrahedron(); break; }
                    case 6  : { generateOctahedron(); break; }
                    case 8  : { generateHexahedron(); break; }
                    case 12 : { generateIcosahedron(); break; }
                    case 20 : { generateDodecahedron(); break; }
                    default : { generatePolyhedron(nplws); break; }
                }
            }
        }

        //! @brief The planewaves destructor.
        virtual ~ProcessorPlanewaves() noexcept { m_planewaves.clear(); }

        //! @brief Returns the number of planewaves.
        inline size_t getNumberOfPlanewaves() const noexcept { return m_planewaves.size(); }

        //! @brief Sets the global rotation of the planewaves.
        //! @param x_axe The angle of rotation around the x axis in radian.
        //! @param y_axe The angle of rotation around the y axis in radian.
        //! @param z_axe The angle of rotation around the z axis in radian.
        void setPlanewavesRotation(const T x_axe, const T y_axe, const T z_axe) noexcept
        {
            m_rotation_x = math<T>::wrap_two_pi(x_axe);
            m_rotation_y = math<T>::wrap_two_pi(y_axe);
            m_rotation_z = math<T>::wrap_two_pi(z_axe);
        }

        //! @brief Returns the angle of rotation around the x axis in radian of the planewaves.
        inline T getPlanewavesRotationX() const noexcept { return m_rotation_x; }

        //! @brief Returns the angle of rotation around the y axis in radian of the planewaves.
        inline T getPlanewavesRotationY() const noexcept { return m_rotation_y; }

        //! @brief Returns the angle of rotation around the z axis in radian of the planewaves.
        inline T getPlanewavesRotationZ() const noexcept { return m_rotation_z; }

        //! @brief Returns the index of a planewave.
        inline size_t getPlanewaveIndex(const size_t index)
        {
            return m_planewaves[index].getIndex();
        }

        //! @brief Sets the azimuth of a planewave in radian.
        //! @details The azimuth of the planewaves is in radian, 0 radian is at the front of
        //! the soundfield and π is at the back of the sound field.
        //! @param index   The index of the planewave.
        //! @param azimuth The azimuth of the planewave.
        inline void setPlanewaveAzimuth(const size_t index, const T azimuth)
        {
            m_planewaves[index].setAzimuth(azimuth);
        }
        
        //! @brief Sets the azimuth of a planewave in radian.
        //! @details The elevation of the planewaves is in radian, 0 radian is at the equator
        //! of the soundfield, π/2 radian is the top of the sound field and -π/2 radian is the
        //! bottom of the sound field.
        //! @param index   The index of the planewave.
        //! @param elevation The elevation of the planewave.
        inline void setPlanewaveElevation(const size_t index, const T elevation)
        {
            m_planewaves[index].setElevation(elevation);
        }

        //! @brief Returns the azimuth of a planewave.
        //! @details The azimuth of the planewaves is in radian, 0 radian is at the front of
        //! the soundfield and π is at the back of the sound field.
        //! @param index    The index of the planewave.
        //! @param rotation false if you don't want to consider the rotation, otherwise true (default).
        inline T getPlanewaveAzimuth(const size_t index, const bool rotation = true) const
        {
            if(rotation)
            {
                T x, y, z;
                Planewave<D, T>::rotate(m_planewaves[index], m_rotation_x, m_rotation_y, m_rotation_z, x, y, z);
                return Planewave<D, T>::xyz2azimuth(x, y, z);
            }
            
            return m_planewaves[index].getAzimuth();
        }

        //! @brief Returns the elevation of a planewave.
        //! @details The elevation of the planewaves is in radian, 0 radian is at the equator
        //! of the soundfield, π/2 radian is the top of the sound field and -π/2 radian is the
        //! bottom of the sound field.
        //! @param index    The index of the planewave.
        //! @param rotation false if you don't want to consider the rotation, otherwise true (default).
        inline T getPlanewaveElevation(const size_t index, const bool rotation = true) const
        {
            if(rotation)
            {
                T x, y, z;
                Planewave<D, T>::rotate(m_planewaves[index], m_rotation_x, m_rotation_y, m_rotation_z, x, y, z);
                return Planewave<D, T>::xyz2elevation(x, y, z);
            }
            
            return m_planewaves[index].getElevation();
        }

        //! @brief Returns the abscissa of a planewave.
        //! @details The abscissa is between -1 and 1, -1 is the left of the soundfield, 0 is
        //! the center of the soundfield and 1 is the right of the soundfield.
        //! @param index    The index of the planewave.
        //! @param rotation false if you don't want to consider the rotation, otherwise true (default).
        inline T getPlanewaveAbscissa(const size_t index, const bool rotation = true) const
        {
            if(rotation)
            {
                T x, y, z;
                Planewave<D, T>::rotate(m_planewaves[index], m_rotation_x, m_rotation_y, m_rotation_z, x, y, z);
                return x;
            }
            
            return m_planewaves[index].getAbscissa();
        }

        //! @brief Returns the ordinate of a planewave.
        //! @details The ordinate is between -1 and 1, -1 is the back of the soundfield, 0 is
        //! the center of the soundfield and 1 is the front of the soundfield.
        //! @param index    The index of the planewave.
        //! @param rotation false if you don't want to consider the rotation, otherwise true (default).
        inline T getPlanewaveOrdinate(const size_t index, const bool rotation = true) const
        {
            if(rotation)
            {
                T x, y, z;
                Planewave<D, T>::rotate(m_planewaves[index], m_rotation_x, m_rotation_y, m_rotation_z, x, y, z);
                return y;
            }
            
            return m_planewaves[index].getOrdinate();
        }

        //! @brief Returns the height of a planewave.
        //! @details The height is between -1 and 1, -1 is the back of the soundfield, 0 is
        //! the center of the soundfield and 1 is the front of the soundfield.
        //! @param index    The index of the planewave.
        //! @param rotation false if you don't want to consider the rotation, otherwise true (default).
        inline T getPlanewaveHeight(const size_t index, const bool rotation = true) const
        {
            if(rotation)
            {
                T x, y, z;
                Planewave<D, T>::rotate(m_planewaves[index], m_rotation_x, m_rotation_y, m_rotation_z, x, y, z);
                return z;
            }
            
            return m_planewaves[index].getHeight();
        }
        
        //! @brief Retruns the name for a planewave.
        //! @details The a name for a planewave is in a std::string format that will be
        //! "Planewave index azimuth (in degrees)".
        //! @param index The index of a planewave.
        inline std::string getPlanewaveName(const size_t index) const noexcept
        {
            return m_planewaves[index].getName();
        }
        
    private:

        void generateCicle(size_t nplws) noexcept
        {
            for(size_t i = 0; i < nplws; i++)
            {
                m_planewaves.push_back(Planewave<D, T>(i+1, static_cast<T>(i) / static_cast<T>(nplws) * static_cast<T>(HOA_2PI), 0.));
            }
        }
        
        void generateTetrahedron() noexcept
        {
            const double oh = -(sqrt(2. / 3.) / sqrt(3. / 8.) - 1.);
            const double hc = sqrt(1. - oh * oh);
            const double el = asin(oh / sqrt(hc*hc + oh*oh));
            m_planewaves.push_back(Planewave<D, T>(1ul,  0., HOA_PI2));
            m_planewaves.push_back(Planewave<D, T>(2ul, 0., el));
            m_planewaves.push_back(Planewave<D, T>(3ul, T(HOA_2PI / 3.), el));
            m_planewaves.push_back(Planewave<D, T>(4ul, T(HOA_2PI / 1.5), el));
        }
        
        void generateHexahedron() noexcept
        {
            m_planewaves.push_back(Planewave<D, T>(1ul, 1., 1., 1.));
            m_planewaves.push_back(Planewave<D, T>(2ul, -1., 1., 1.));
            m_planewaves.push_back(Planewave<D, T>(3ul, -1., -1., 1.));
            m_planewaves.push_back(Planewave<D, T>(4ul, 1., -1., 1.));
            
            m_planewaves.push_back(Planewave<D, T>(5ul, 1., 1., -1.));
            m_planewaves.push_back(Planewave<D, T>(6ul, -1., 1., -1.));
            m_planewaves.push_back(Planewave<D, T>(7ul, -1., -1., -1.));
            m_planewaves.push_back(Planewave<D, T>(8ul, 1., -1., -1.));
        }
        
        void generateOctahedron() noexcept
        {
            m_planewaves.push_back(Planewave<D, T>(1ul, 0., HOA_PI2));
            m_planewaves.push_back(Planewave<D, T>(2ul, 0., 0.));
            m_planewaves.push_back(Planewave<D, T>(3ul, HOA_PI2, 0.));
            m_planewaves.push_back(Planewave<D, T>(4ul, (2. * HOA_PI2), 0.));
            m_planewaves.push_back(Planewave<D, T>(5ul, (3. * HOA_PI2), 0.));
            m_planewaves.push_back(Planewave<D, T>(6ul, 0., -HOA_PI2));
        }
        
        void generateIcosahedron() noexcept
        {
            m_planewaves.push_back(Planewave<D, T>(1, 0., HOA_PI2));
            for(size_t i = 1; i < 6; i++)
            {
                m_planewaves.push_back(Planewave<D, T>(i*2, static_cast<T>(i - 1.) / 5. * HOA_2PI, atan(0.5)));
                m_planewaves.push_back(Planewave<D, T>(i*2+1, T(i - 1.) / 5. * HOA_2PI - HOA_PI / 5., -atan(0.5)));
            }
            m_planewaves.push_back(Planewave<D, T>(12, 0., static_cast<T>(-HOA_PI2)));
        }
        
        void generateDodecahedron() noexcept
        {
            const T phi = (sqrt(5.) - 1.) / 2.;
            const T R = 1. / sqrt(3.);
            const T a = R;
            const T b = R / phi;
            const T c = R * phi;
            size_t index = 1;
            for(long i = -1; i < 2; i += 2)
            {
                for(long j = -1; j < 2; j += 2)
                {
                    m_planewaves.push_back(Planewave<D, T>(index++, 0., i * c * R, -j * b * R));
                    m_planewaves.push_back(Planewave<D, T>(index++, i * c * R, j * b * R, 0.));
                    m_planewaves.push_back(Planewave<D, T>(index++, i * b * R, 0., -j * c * R));
                    for(long k = -1; k < 2; k += 2)
                    {
                        m_planewaves.push_back(Planewave<D, T>(index++, i * a * R, j * a * R, -k * a * R));
                    }
                }
            }
        }
        
        void generatePolyhedron(size_t nplws) noexcept
        {
            if(nplws % 2)
            {
                m_planewaves.push_back(Planewave<D, T>(1, 0., static_cast<T>(HOA_PI2)));
            }
            const T     phi     = (sqrt(5.) - 1.) / 4.;
            const size_t limit   = (nplws - (nplws % 2)) / 2;
            const T     offset  = 1. / T(limit) * HOA_PI;
            for(size_t i = 0; i < nplws - (nplws % 2); i++)
            {
                if(i < limit)
                {
                    m_planewaves.push_back(Planewave<D, T>(i+1, T(i) / T(limit) * HOA_2PI, phi));
                }
                else
                {
                    m_planewaves.push_back(Planewave<D, T>(i+1, T(i - (limit)) / T(limit) * HOA_2PI + offset, phi * 2));
                }
            }
        }
        
    private:
        
        std::vector<Planewave<D, T>> m_planewaves {};
        T m_rotation_z = 0.;
        T m_rotation_y = 0.;
        T m_rotation_x = 0.;
    };
}
