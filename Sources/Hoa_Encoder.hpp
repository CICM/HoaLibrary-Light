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

#ifndef DEF_HOA_ENCODER_LIGHT
#define DEF_HOA_ENCODER_LIGHT

#include "Hoa_Processor.hpp"

namespace hoa
{

    //! @brief The class encodes a signal into the harmonics domain depending on coodinates.
    //! @details The class generates the signals associated to the harmonics \f$Y_{l,m}\f$
    //! according to a radius \f$\rho\f$, an azimuth \f$\theta\f$ and an elevation \f$\varphi\f$.\n
    //! The coefficients of the harmonics are defined by:\n
    //! \f[Y_{l,m}(\theta, \varphi) = G_{l,m}(\rho) P_{l, \left|m\right|}(\cos{(\varphi)}) e^{+im\theta} k_{l, m} \f]
    //! with\n
    //! \f$G_{l,m}(\rho)\f$ the radius part of the equation,\n
    //! \f$e^{+im\theta}\f$ the azimuth part with \f$i\f$ the imaginary,\n
    //! \f$P_{l, \left|m\right|}(\cos{(\varphi)})\f$ the elevation part of the equation with
    //! \f$P_{l, \left|m\right|}(x)\f$ the associated Legendre polynomials,\n
    //! \f$k_{l, m}\f$ the normalization and\n
    //! \f$N\f$ the order of decomposition, \f$l\f$ the degree, \f$m\f$ the azimuthal order,
    //! \f$\rho\f$ the radius, \f$\theta\f$ the azimuth and \f$\varphi\f$ the elevation in radian.\n\n
    //! The radius \f$\rho\f$ is included between \f$0\f$ and \f$+\infty\f$. \f$0\f$ is the
    //! center of the soundfield, \f$1\f$ is the radius of the ambisonics circle or sphere,
    //! beyond this limit the gains \f$G\f$ of the harmonics decrease globally and before the
    //! gains \f$G\f$ decrease independently:\n
    //! if \f$\rho\ge1\f$
    //! \f[G_{l,m}(\rho) = \frac{1}{\rho}\f]
    //! else
    //! \f[G_{l,m}(\rho) = \rho^l((1-\rho)(N-l)+1)\f]
    //! The azimuth \f$\theta\f$ in radian is included between \f$0\f$ and \f$2\pi\f$. The
    //! direction of rotation is counterclockwise. The \f$0\f$ radian is \f$\frac{\pi}{2}\f$
    //! phase shifted relative to a mathematical representation of a circle, then the \f$0\f$
    //! radian is at the "front" of the soundfield. The azimuth part in the imaginary form of
    //! the equation \f$e^{+im\theta}\f$ can be expressed with the real form :\n
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
    //! the recursives formulas:\n
    //! \f[P_{l+1,l+1}(x) = -(2l+1)\sqrt{(1-x^2)}P_{(l,l)}(x) \f]
    //! \f[P_{l+1,l}(x) = x(2l+1)P_{(l,l)}(x) \f]
    //! \f[P_{l+1,m}(x) = \frac{x(2l+1)P_{(l,m)}(x) - (l+m)P_{(l-1,m)}(x)}{l-m+1} \f]
    //! and with \f[P_{0, 0}(x) = 1\f]
    //! The normalization part \f$k_{l, m}\f$ is equivalent to :\n
    //! if \f$m = 0\f$ then
    //! \f[k_{l, m} = 1\f]
    //! else
    //! \f[k_{l, m} = \sqrt{\frac{(l - \left|m\right|)!}{(l + \left|m\right|)!}}\sqrt{2} \f]
    //! @todo Check the azimuth offset
    template <Dimension D, typename T> class Encoder : public ProcessorHarmonics<D, T>
    {
    public:

        //! @brief The constructor.
        //! @param order The order of decomposition.
        Encoder(size_t order) : ProcessorHarmonics<D, T>(order)
        {
            order = order >= 1 ? order : 1;
            size_t const nharmo = Harmonic<Hoa3d, T>::getNumberOfHarmonics(order);
            m_radius_coeffs         = new T[order + 1];
            //m_azimuth_coeffs        = new T[order * 2 + 1];
            //m_elevation_coeffs      = new T[nharmo / 2 + (order + 1)];
            m_normalization_coeffs  = new T[nharmo];
            setRadius(1.);
            //setAzimuth(0.);
            //setElevation(0.);
            
            for(size_t i = 0; i < nharmo; ++i)
            {
                m_normalization_coeffs[i] = ProcessorHarmonics<D, T>::getHarmonicSemiNormalization(i);
            }
        }

        //! The destructor.
        ~Encoder() hoa_noexcept
        {
            delete [] m_radius_coeffs;
            delete [] m_azimuth_coeffs;
            delete [] m_elevation_coeffs;
            delete [] m_normalization_coeffs;
        }
        
        
        //! @brief Returns the radius.
        inline T getRadius() const hoa_noexcept { return m_radius; }
        
        //! @brief Returns the azimuth.
        inline T getAzimuth() const hoa_noexcept { return m_azimuth; }
        
        //! @brief Returns the elevation.
        inline T getElevation()  const hoa_noexcept { return m_elevation; }
        
        //! @brief Sets the radius.
        //! @param radius The new radius.
        //! @param azimuth The new azimuth.
        //! @param elevation The new elevation.
        inline void setCoordinates(T radius, T azimuth, T elevation) hoa_noexcept {
            setRadius(radius);
            elevation = Math<T>::wrap_pi(elevation);
            if(elevation >= -HOA_PI2 && elevation <= HOA_PI2)
            {
                setAzimuth(azimuth);
                setElevation(elevation);
            }
            else
            {
                setAzimuth(azimuth + HOA_PI);
                setElevation(elevation);
            }
        }
        
        //! @brief Sets the radius.
        //! @param radius The new radius.
        void setRadius(T radius) hoa_noexcept {
            m_radius = std::max(radius, T(0.));
            
            if(m_radius >= 1)
            {
                const size_t order  = ProcessorHarmonics<D, T>::getDecompositionOrder();
                T* coeff     = m_radius_coeffs;
                const T gain = T(1) / m_radius;
                (*coeff++) = gain;
                for(size_t i = 1; i <= order; ++i)
                {
                    (*coeff++) = gain;
                }
            }
            else
            {
                const size_t order  = ProcessorHarmonics<D, T>::getDecompositionOrder();
                T* coeff            = m_radius_coeffs;
                const T widening    = m_radius;
                const T temp        = 1 - widening;
                (*coeff++) = (T(order) * temp + T(1.));
                for(size_t i = 1; i <= order; ++i)
                {
                    (*coeff++) = std::pow(widening, T(i)) * (temp * T(order - i) + T(1.));
                }
            }
            
        }
        
        //! @brief Sets the azimuth.
        //! @param azimuth The new azimuth.
        void setAzimuth(T azimuth) hoa_noexcept {
            m_azimuth   = azimuth;
            T* coeffs   = m_azimuth_coeffs;
            const size_t order  = ProcessorHarmonics<D, T>::getDecompositionOrder();
            T const ccos_x  = std::cos(azimuth);
            T const csin_x  = std::sin(azimuth);
            if(D == Hoa2d)
            {
                // Organization 0, -1, 1, -2, 2, -3, 3,...
                T cos_x     = ccos_x;
                T sin_x     = csin_x;
                T tcos_x    = cos_x;
                (*coeffs++) = 1;                                // Hamonic [0]
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
        void setElevation(T elevation) hoa_noexcept {
            m_elevation = elevation;
            T* coeffs   = m_elevation_coeffs;
            const size_t order = ProcessorHarmonics<D, T>::getDecompositionOrder();
            const T cos_theta = T(std::sin(elevation)); // In fact we use sinus to offset the elvation
            const T sqr_theta = -T(std::sqrt(T(1.) - cos_theta * cos_theta));
            // Organization [0, 0], [1, 1], [1, 0], [2, 2], [2, 1], [2, 0], [3, 3],
            // [3, 2], [3, 1], [3, 0], [4, 4], ...
            
            *(coeffs++) = 1;            // P(0, 0)(x) = 1
            *(coeffs++) = sqr_theta;    // P(1, 1)(x) = -(2*0+1)sqrt(1-x^2)P(0,0)(x) = -sqrt(1-x^2)
            *(coeffs++) = cos_theta;    // P(1, 0)(x) = x(2*0+1)P(0,0)(x) = x
            for(size_t i = 2; i <= order; ++i)
            {
                const T ratio = 2 * T(i-1) + 1;
                const T previous = *(coeffs-i);
                // P(l+1,l+1)(x) = -(2l+1) sqrt(1-x^2) P(l,l)(x)
                *(coeffs++) = ratio * sqr_theta * previous;
                // P(l+1,l)(x)   = x(2l+1) P(l,l)(x)
                *(coeffs++) = ratio * cos_theta * previous;
                for(size_t j = i - 2; j > 0; --j)
                {
                    const T previous1 = *(coeffs-(i+1));
                    const T previous2 = *(coeffs-(2*i+1));
                    // P(l+1,m)(x)   = (x(2l+1)P(l,m)(x) - (l+m)P(l-1, m)(x)) / (l-m+1)
                    *(coeffs++) = ((ratio * cos_theta * previous1) - (T(i-1) + T(j)) * previous2) / (T(i-1) - T(j) + 1);
                }
                const T previous1 = *(coeffs-(i+1));
                const T previous2 = *(coeffs-(2*i+1));
                // P(l+1,0)(x)   = (x(2l+1)P(l,0)(x) - (l)P(l-1, 0)(x)) / (l+1)
                *(coeffs++) = ((ratio * cos_theta * previous1) - T(i-1) * previous2) / (T(i-1) + 1);
            }
        }
    
        //! @brief The method performs the encoding of the harmonics signal.
        //! @details The input pointer must be the sample to encode and the outputs array
        //! contains the spherical harmonics samples thus the minimum size of the array must
        //! be the number of harmonics.
        //! @param input   The input pointer.
        //! @param outputs The outputs array.
        void process(const T* input, T* outputs) hoa_noexcept hoa_override
        {
            if(D == Hoa2d)
            {
                const size_t order      = ProcessorHarmonics<D, T>::getDecompositionOrder();
                T const* radius_coeffs  = m_radius_coeffs;
                T const* azimuth_coeffs = m_azimuth_coeffs;
                (*outputs++) = (*input) * (*radius_coeffs++) * (*azimuth_coeffs++);
                for(size_t i = 1; i <= order; ++i)
                {
                    const T      factor = (*radius_coeffs++);
                    (*outputs++) = (*input) * factor * (*azimuth_coeffs++);
                    (*outputs++) = (*input) * factor * (*azimuth_coeffs++);
                }
            }
            else
            {
                const size_t order          = ProcessorHarmonics<D, T>::getDecompositionOrder();
                T const* radius_coeffs      = m_radius_coeffs;
                T const* azimuth_coeffs     = m_azimuth_coeffs+order;
                T const* elevation_coeffs   = m_elevation_coeffs;
                T const* norm_coeffs        = m_normalization_coeffs;
                (*outputs++) = (*input) * (*radius_coeffs++) * (*azimuth_coeffs--) * (*elevation_coeffs++) * (*norm_coeffs++);
                for(size_t i = 1; i <= order; ++i)
                {
                    const T factor      = (*radius_coeffs++);
                    T const* az_coeffs2 = azimuth_coeffs--;
                    T const* el_coeffs2 = elevation_coeffs;
                    elevation_coeffs = elevation_coeffs+(i+1);
                    for(size_t j = i; j; --j)
                    {
                        //(*outputs++) = (*el_coeffs2++);
                        (*outputs++) = (*input) * factor * (*az_coeffs2++) * (*el_coeffs2++) * (*norm_coeffs++);
                    }
                    //(*outputs++) = (*el_coeffs2--);
                    (*outputs++) = (*input) * factor * (*az_coeffs2++) * (*el_coeffs2--) * (*norm_coeffs++);
                    for(size_t j = i; j; --j)
                    {
                        //(*outputs++) = (*el_coeffs2--);
                        (*outputs++) = (*input) * factor * (*az_coeffs2++) * (*el_coeffs2--) * (*norm_coeffs++);
                    }
                }
            }
            
        }
        
    private:

        T   m_radius;
        T   m_azimuth;
        T   m_elevation;
        T*  m_radius_coeffs;
        T*  m_azimuth_coeffs;
        T*  m_elevation_coeffs;
        T*  m_normalization_coeffs;
    };
}

#endif
