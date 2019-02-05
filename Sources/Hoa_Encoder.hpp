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

#include "Hoa_Processor.hpp"
#include "Hoa_Math.hpp"

namespace hoa
{
    // ================================================================================ //
    // ENCODER //
    // ================================================================================ //
    
    //! @brief The class encodes a signal into the harmonics domain depending on coodinates.
    //! @details The class generates the signals associated to the harmonics \f$Y_{l,m}\f$
    //! according to a radius \f$\rho\f$, an azimuth \f$\theta\f$ and an elevation \f$\varphi\f$.<br>
    //! The coefficients of the harmonics are defined by:<br>
    //! \f[Y_{l,m}(\theta, \varphi) = G_{l,m}(\rho) P_{l, \left|m\right|}(\sin{(\varphi)}) e^{+im\theta} k_{l, m} \f]
    //! with<br>
    //! \f$G_{l,m}(\rho)\f$ the radius part of the equation,<br>
    //! \f$e^{+im\theta}\f$ the azimuth part with \f$i\f$ the imaginary,<br>
    //! \f$P_{l, \left|m\right|}(\cos{(\varphi)})\f$ the elevation part of the equation with
    //! \f$P_{l, \left|m\right|}(x)\f$ the associated Legendre polynomials,<br>
    //! \f$k_{l, m}\f$ the normalization and<br>
    //! \f$N\f$ the order of decomposition, \f$l\f$ the degree, \f$m\f$ the azimuthal order,
    //! \f$\rho\f$ the radius, \f$\theta\f$ the azimuth and \f$\varphi\f$ the elevation in radian.<br><br>
    //! The radius \f$\rho\f$ is included between \f$0\f$ and \f$+\infty\f$. \f$0\f$ is the
    //! center of the soundfield, \f$1\f$ is the radius of the ambisonics circle or sphere,
    //! beyond this limit the gains \f$G\f$ of the harmonics decrease globally and before the
    //! gains \f$G\f$ decrease independently:<br>
    //! if \f$\rho\ge1\f$
    //! \f[G_{l,m}(\rho) = \frac{1}{\rho}\f]
    //! else
    //! \f[G_{l,m}(\rho) = \rho^l((1-\rho)(N-l)+1)\f]
    //! The azimuth \f$\theta\f$ in radian is included between \f$0\f$ and \f$2\pi\f$. The
    //! direction of rotation is counterclockwise. The \f$0\f$ radian is \f$\frac{\pi}{2}\f$
    //! phase shifted relative to a mathematical representation of a circle, then the \f$0\f$
    //! radian is at the "front" of the soundfield. The azimuth part in the imaginary form of
    //! the equation \f$e^{+im\theta}\f$ can be expressed with the real form :<br>
    //! if \f$m\geq0\f$
    //! \f[e^{+im\theta} = \cos{(\left|m\right|\theta)}\f]
    //! else
    //! \f[e^{+im\theta} = \sin{(\left|m\right|\theta)}\f]
    //! The elevation \f$\varphi\f$ is included between \f$-\pi\f$ and \f$\pi\f$. The
    //! direction of rotation is from bottom to the top. The \f$0\f$ radian is centered at the
    //! "front" of the soundfield, then \f$\frac{\pi}{2}\f$ is at the top, \f$-\frac{\pi}{2}\f$
    //! is at the bottom and \f$\pi\f$ is behind. Note that if the angle of elevation is
    //! between \f$\frac{\pi}{2}\f$ and \f$\frac{3\pi}{2}\f$, the azimuth is reversed. The
    //! elevation part \f$P_{l, \left|m\right|}(x)\f$ of the formula can be expressed with
    //! the recursives formulas:<br>
    //! \f[P_{l+1,l+1}(x) = -(2l+1)\sqrt{(1-x^2)}P_{(l,l)}(x) \f]
    //! \f[P_{l+1,l}(x) = x(2l+1)P_{(l,l)}(x) \f]
    //! \f[P_{l+1,m}(x) = \frac{x(2l+1)P_{(l,m)}(x) - (l+m)P_{(l-1,m)}(x)}{l-m+1} \f]
    //! and with \f[P_{0, 0}(x) = 1\f]
    //! The normalization part \f$k_{l, m}\f$ is equivalent to :<br>
    //! if \f$m = 0\f$ then
    //! \f[k_{l, m} = 1\f]
    //! else
    //! \f[k_{l, m} = \sqrt{\frac{(l - \left|m\right|)!}{(l + \left|m\right|)!}}\sqrt{2} \f]
    template <Dimension D, typename T>
    class Encoder
    : public ProcessorHarmonics<D, T>
    {
    public:

        //! @brief Constructor.
        //! @param order The order of decomposition.
        Encoder(const size_t order)
        : ProcessorHarmonics<D, T>(order)
        , m_radius_coeffs(order+1)
        , m_azimuth_coeffs(order*2+1+2)
        , m_elevation_coeffs(((order + 1) * (order + 1)) / 2 + (order + 1)+3)
        , m_normalization_coeffs(ProcessorHarmonics<D, T>::getNumberOfHarmonics())
        {
            setRadius(1.);
            setAzimuth(0.);
            setElevation(0.);
            
            for(size_t i = 0; i < ProcessorHarmonics<D, T>::getNumberOfHarmonics(); ++i)
            {
                m_normalization_coeffs[i] = ProcessorHarmonics<D, T>::getHarmonicSemiNormalization(i) * std::pow(static_cast<T>(-1), static_cast<T>(ProcessorHarmonics<D, T>::getHarmonicOrder(i)));
            }
        }

        //! Destructor.
        ~Encoder() = default;
        
        //! @brief Returns the radius.
        inline T getRadius() const noexcept { return m_radius; }
        
        //! @brief Returns the azimuth.
        inline T getAzimuth() const noexcept { return m_azimuth; }
        
        //! @brief Returns the elevation.
        inline T getElevation()  const noexcept { return m_elevation; }
        
        //! @brief Sets the coordinates.
        //! @param radius The new radius.
        //! @param azimuth The new azimuth.
        //! @param elevation The new elevation.
        inline void setCoordinates(const T radius, const T azimuth, const T elevation) noexcept
        {
            setRadius(radius);
            setAzimuth(azimuth);
            setElevation(elevation);
        }
        
        //! @brief Sets the radius.
        //! @param radius The new radius.
        void setRadius(const T radius) noexcept
        {
            m_radius = std::max(radius, static_cast<T>(0));
            if(m_radius >= 1)
            {
                std::fill(m_radius_coeffs.begin(), m_radius_coeffs.end(), static_cast<T>(1) / m_radius);
            }
            else
            {
                const size_t order  = ProcessorHarmonics<D, T>::getDecompositionOrder();
                const T widening    = m_radius;
                const T temp        = 1. - widening;
                T* coeff            = m_radius_coeffs.data();
                (*coeff++) = static_cast<T>(order) * temp + static_cast<T>(1);
                for(size_t i = 1; i <= order; ++i)
                {
                    (*coeff++) = (std::pow(widening, static_cast<T>(i))
                                  * (temp * static_cast<T>(order - i) + static_cast<T>(1.)));
                }
            }
            
        }
        
        //! @brief Sets the azimuth.
        //! @param azimuth The new azimuth.
        void setAzimuth(const T azimuth) noexcept
        {
            m_azimuth   = math<T>::wrap_two_pi(azimuth);
            T* coeffs   = m_azimuth_coeffs.data();
            const size_t order  = ProcessorHarmonics<D, T>::getDecompositionOrder();
            T const ccos_x  = std::cos(azimuth);
            T const csin_x  = std::sin(azimuth);
            if(D == Hoa2d)
            {
                // Organization 0, -1, 1, -2, 2, -3, 3,...
                T cos_x     = ccos_x;
                T sin_x     = csin_x;
                T tcos_x    = cos_x;
                (*coeffs++) = 1;                               // Hamonic [0]
                (*coeffs++) = sin_x;                           // Hamonic [-1]
                (*coeffs++) = cos_x;                           // Hamonic [1]
                for(size_t i = 2; i <= order; i++)
                {
                    cos_x   = tcos_x * ccos_x - sin_x * csin_x;
                    sin_x   = tcos_x * csin_x + sin_x * ccos_x;
                    tcos_x  = cos_x;
                    (*coeffs++) = sin_x;                        // Hamonic [l,-l]
                    (*coeffs++) = cos_x;                        // Hamonic [l,l]
                }
            }
            else
            {
                // Organization ..., -3, -2, -1, 0, 1, 2, 3, ...
                T cos_x     = ccos_x;
                T sin_x     = csin_x;
                T tcos_x    = cos_x;
                coeffs[order] = 1;                                 // Hamonic [0]
                coeffs[order-1] = sin_x;                           // Hamonic [-1]
                coeffs[order+1] = cos_x;                           // Hamonic [1]
                for(size_t i = 2; i <= order; i++)
                {
                    cos_x   = tcos_x * ccos_x - sin_x * csin_x;
                    sin_x   = tcos_x * csin_x + sin_x * ccos_x;
                    tcos_x  = cos_x;
                    coeffs[order-i] = sin_x;                        // Hamonic [l,-l]
                    coeffs[order+i] = cos_x;                        // Hamonic [l,l]
                }
            }
        }
        
        //! @brief Sets the elevation.
        //! @param elevation The new elevation.
        void setElevation(const T elevation) noexcept
        {
            m_elevation = math<T>::wrap_pi(elevation);
            const size_t order = ProcessorHarmonics<D, T>::getDecompositionOrder();
            const T _x      = std::sin(elevation);
            const T _x_pow  = _x * _x;                    // x^2
            const T _x_powm = static_cast<T>(1) - _x_pow; // 1 - x^2
            const T _x_sqpm = std::sqrt(_x_powm);         // sqrt(1 - x^2)
            
            // Organization [0, 0], [1, 1], [1, 0], [2, 2], [2, 1], [2, 0], [3, 3], [3, 2], [3, 1], [3, 0], [4, 4], ...
            // P(l+1, l+1)(x) = -(2l+1)sqrt(1-x^2)P(l,l)(x)
            // P(l+1,l)(x)    = x(2l+1)P(l,l)(x)
            // P(l+1,m)(x)    = x(2l+1)P(l,m)(x) - (l+m)P(l-1,m)(x) / (l-m+1)
            
            T* coeffs   = m_elevation_coeffs.data();
            *(coeffs++) = 1;
            *(coeffs++) = -_x_sqpm;
            *(coeffs++) = _x;
            
            if(order <= 1) { return; }
            *(coeffs++) = static_cast<T>(3) * _x_powm;
            *(coeffs++) = static_cast<T>(-3) * _x * _x_sqpm;
            *(coeffs++) = static_cast<T>(1.5) * _x_pow - static_cast<T>(0.5);
            
            if(order <= 2) { return; }
            *(coeffs++) = static_cast<T>(-15) * _x_powm * _x_sqpm;
            *(coeffs++) = static_cast<T>(15) * _x_powm * _x;
            *(coeffs++) = (static_cast<T>(-7.5) * _x_pow + static_cast<T>(1.5)) * _x_sqpm;
            *(coeffs++) = (static_cast<T>(2.5) * _x_pow + static_cast<T>(-1.5)) * _x;
            
            if(order <= 3) { return; }
            *(coeffs++) = static_cast<T>(105) * _x_powm * _x_powm;
            *(coeffs++) = static_cast<T>(-105) * _x_powm * _x_sqpm * _x;
            *(coeffs++) = (static_cast<T>(52.5) * _x_pow - static_cast<T>(7.5)) * _x_powm;
            *(coeffs++) = (static_cast<T>(-17.5) * _x_pow + static_cast<T>(7.5)) * _x_sqpm  * _x;
            *(coeffs++) = (static_cast<T>(4.375) * _x_pow - static_cast<T>(3.75)) * _x_pow + static_cast<T>(0.375);
            
            if(order <= 4) { return; }
            *(coeffs++) = static_cast<T>(-945) * _x_powm * _x_powm * _x_sqpm;
            *(coeffs++) = static_cast<T>(945) * _x_powm * _x_powm * _x;
            *(coeffs++) = (static_cast<T>(52.5) - static_cast<T>(472.5) * _x_pow) * _x_powm * _x_sqpm;
            *(coeffs++) = (static_cast<T>(157.5) * _x_pow - static_cast<T>(52.5)) * _x_powm * _x;
            *(coeffs++) = ((static_cast<T>(-39.375) * _x_pow + static_cast<T>(26.25)) * _x_pow - static_cast<T>(1.875)) * _x_sqpm;
            *(coeffs++) = ((static_cast<T>(7.875) * _x_pow - static_cast<T>(8.75)) * _x_pow + static_cast<T>(1.875)) * _x;
            
            for(size_t i = 6; i <= order; ++i)
            {
                const T ratio = static_cast<T>(2 * (i-1) + 1);
                const T previous = *(coeffs-i);
                // P(l+1,l+1)(x) = -(2l+1) sqrt(1-x^2) P(l,l)(x)
                *(coeffs++) = ratio * -_x_sqpm * previous;
                // P(l+1,l)(x)   = x(2l+1) P(l,l)(x)
                *(coeffs++) = ratio * _x * previous;
                for(size_t j = i - 2; j > 0; --j)
                {
                    const T previous1 = *(coeffs-(i+1));
                    const T previous2 = *(coeffs-(2*i+1));
                    // P(l+1,m)(x)   = (x(2l+1)P(l,m)(x) - (l+m)P(l-1, m)(x)) / (l-m+1)
                    *(coeffs++) = ((ratio * _x * previous1) - static_cast<T>(i-1+j) * previous2) / static_cast<T>(i-j);
                }
                const T previous1 = *(coeffs-(i+1));
                const T previous2 = *(coeffs-(2*i+1));
                // P(l+1,0)(x)   = (x(2l+1)P(l,0)(x) - (l)P(l-1, 0)(x)) / (l+1)
                *(coeffs++) = ((ratio * _x * previous1) - static_cast<T>(i-1) * previous2) / static_cast<T>(i);
            }
        
        }
    
        //! @brief The method performs the encoding of the harmonics signal.
        //! @details The input pointer must be the sample to encode and the outputs array
        //! contains the spherical harmonics samples thus the minimum size of the array must
        //! be the number of harmonics.
        //! @param input   The input pointer.
        //! @param outputs The outputs array.
        void process(const T* input, T* outputs) noexcept override
        {
            if(D == Hoa2d)
            {
                const size_t order      = ProcessorHarmonics<D, T>::getDecompositionOrder();
                T const* radius_coeffs  = m_radius_coeffs.data();
                T const* azimuth_coeffs = m_azimuth_coeffs.data();
                (*outputs++) = (*input) * (*radius_coeffs++) * (*azimuth_coeffs++);
                for(size_t i = 1; i <= order; ++i)
                {
                    const T factor = (*radius_coeffs++);
                    (*outputs++) = (*input) * factor * (*azimuth_coeffs++);
                    (*outputs++) = (*input) * factor * (*azimuth_coeffs++);
                }
            }
            else
            {
                const size_t order          = ProcessorHarmonics<D, T>::getDecompositionOrder();
                T* tout = outputs;
                T const* radius_coeffs      = m_radius_coeffs.data();
                T const* azimuth_coeffs     = m_azimuth_coeffs.data()+order;
                T const* elevation_coeffs   = m_elevation_coeffs.data();
                T const* norm_coeffs        = m_normalization_coeffs.data();
                (*outputs++) = (*input) * (*radius_coeffs++) * (*azimuth_coeffs--) * (*elevation_coeffs++) * (*norm_coeffs++);
                for(size_t i = 1; i <= order; ++i)
                {
                    const T factor      = (*radius_coeffs++);
                    T const* az_coeffs2 = azimuth_coeffs--;
                    T const* el_coeffs2 = elevation_coeffs;
                    elevation_coeffs = elevation_coeffs+(i+1);
                    for(size_t j = i; j; --j)
                    {
                        (*outputs++) = (*input) * factor * (*az_coeffs2++) * (*el_coeffs2++) * (*norm_coeffs++);
                    }
                    (*outputs++) = (*input) * factor * (*az_coeffs2++) * (*el_coeffs2--) * (*norm_coeffs++);
                    for(size_t j = i; j; --j)
                    {
                        (*outputs++) = (*input) * factor * (*az_coeffs2++) * (*el_coeffs2--) * (*norm_coeffs++);
                    }
                }
                
                if(!(m_elevation >= static_cast<T>(-HOA_PI2) && m_elevation <= static_cast<T>(HOA_PI2)))
                {
                    for(size_t i = 1; i < ProcessorHarmonics<D, T>::getNumberOfHarmonics(); i+=2)
                    {
                        tout[i] = -tout[i];
                    }
                }
            }
            
        }
        
    private:
        
        T m_radius = 0.;
        T m_azimuth = 0.;
        T m_elevation = 0.;
        std::vector<T> m_radius_coeffs {};
        std::vector<T> m_azimuth_coeffs {};
        std::vector<T> m_elevation_coeffs {};
        std::vector<T> m_normalization_coeffs {};
    };
}
