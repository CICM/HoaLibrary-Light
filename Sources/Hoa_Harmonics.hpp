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

namespace hoa
{
    // ================================================================================ //
    // HARMONIC //
    // ================================================================================ //
    
    //! @brief The class owns basic harmonics informations.
    //! @details The class allows to retrieves several informations about the harmonics: the
    //! numbering (ACN), the normalization (SN3D/N3D), the degree, the aziumthal order, etc.
    template<Dimension D, typename T>
    class Harmonic
    {
    public:

        //! @brief The harmonic constructor.
        //! @param index The index must be at least 1.
        inline Harmonic(const size_t index) noexcept
        : m_index(index)
        , m_degree(getDegree(index))
        , m_order(getOrder(index))
        {}

        //! @brief The harmonic destructor.
        inline ~Harmonic() noexcept {}

        //! @brief Returns the ACN index of the harmonic.
        inline size_t getIndex() const noexcept { return m_index; }

        //! @brief Returns the precomputed degree of the harmonic.
        inline size_t getDegree() const noexcept { return m_degree; }

        //! @brief Returns the precomputed azimuthal order of the harmonic.
        inline long getOrder() const noexcept { return m_order; }
        
        //! @brief Returns the name of the harmonic.
        std::string getName() const
        {
            std::ostringstream ostr;
            ostr <<  "Harmonic (" << m_index << ") : " << m_degree << " " << m_order;
            return ostr.str();
        }

        //! @brief Returns the normalization of the harmonic.
        inline T getNormalization() const noexcept
        {
            return getNormalization(m_degree, m_order);
        }

        //! @brief Returns the semi-normalization of the harmonic.
        inline T getSemiNormalization() const noexcept
        {
            return getSemiNormalization(m_degree, m_order);
        }

        //! @brief Returns the index of an harmonic for a specific degree and azimuthal order.
        //! @details The computation is \f$i=2|m|-(m<0)\f$ in 2D and \f$i=(l(l+1))+m\f$ in 3D.
        //! with \f$l\f$ the degree, \f$m\f$ the azimuthal order and \f$i\f$ the index of the
        //! harmonic.
        //! @param degree  The degree of the harmonic.
        //! @param order   The order of the harmonic.
        static inline long getIndex(const size_t degree, const long order) noexcept
        {
            return (D == Hoa2d) ? (std::abs(order) *  2 - long(order < 0)) : (size_t(long(degree * (degree + 1)) + order));
        }
        
        //! @brief Returns the degree of an harmonic for an index.
        //! @details The computation is \f$l=\frac{i+i\bmod{2}}{2}\f$ in 2D and \f$l=\sqrt{i}\f$ in 3D.
        //! with \f$l\f$ the degree and \f$i\f$ the index of the harmonic.
        //! @param index  The index of the harmonic.
        static inline size_t getDegree(const size_t index) noexcept
        {
            return (D == Hoa2d) ? ((index + index % 2) / size_t(2)) : size_t(sqrt(double(index)));
        }

        //! @brief Returns the azimuthal order of an harmonic for an index.
        //! @details The computation is \f$m=l(1-2(i\bmod{2}))\f$ in 2D and \f$m=i-(l(l+1))\f$ in 3D.
        //! with \f$l\f$ the degree, \f$m\f$ the azimuthal order and \f$i\f$ the index of the
        //! harmonic.
        //! @param index  The index of the harmonic.
        static inline long getOrder(const size_t index) noexcept
        {
            const long l = long(getDegree(index));
            return (D == Hoa2d) ? (l * (1l - static_cast<long>(index % 2) * 2l)) : (static_cast<long>(index) - (l * (l + 1)));
        }

        //! @brief Returns the number of harmonics for an order of decomposition.
        //! @details The computation is \f$2N+1\f$ in 2D and \f$(N+1)^{2}\f$ in 3D.
        //! with \f$N\f$ the order of decomposition.
        //! @param order   The order of decomposition.
        static inline constexpr size_t getNumberOfHarmonics(const size_t order) noexcept
        {
            return (D == Hoa2d) ? (order * 2 + 1) : (order + 1) * (order + 1);
        }
        
        //! @brief Returns the number of harmonics for a degree.
        //! @details The computation is \f$(degree\neq0)+1\f$ in 2D and \f$2l+1\f$ in 3D.
        //! with \f$l\f$ the degree.
        //! @param degree  The degree.
        static inline constexpr size_t getNumberOfHarmonicsInDegree(const size_t degree) noexcept
        {
            return (D == Hoa2d) ? ((degree != 0) + 1) : (degree * 2 + 1);
        }

        //! @brief Returns the normalization N3D or N2D of an harmonic.
        //! @details The semi-normalization \f$k^{n2d}_{l, m}\f$ and \f$k^{n3d}_{l, m}\f$
        //! are defined by:
        //! \f[k^{n2d}_{l, m} = 1 \f]
        //! \f[k^{n3d}_{l, m} = k^{sn3d}_{l, m}\sqrt{2l+1} \f]
        //! with\f$l\f$ the degree, \f$m\f$ the azimuthal order,
        //! @param degree  The degree of the harmonic.
        //! @param order   The order of the harmonic.
        static inline T getNormalization(const size_t degree, const long order) noexcept
        {
            return (D == Hoa2d) ? 1. : (getSemiNormalization(degree, order) * T(sqrt(T(2.) * T(degree) + T(1.))));
        }

        //! @brief Returns the semi-normalization of an harmonic.
        //! @param degree  The degree of the harmonic.
        //! @param order   The order of the harmonic.
        static T getSemiNormalization(const size_t degree, const long order) noexcept
        {
            if(D == Hoa2d)
                return 1.;
            
            const long double fac1 = hfactorial(long(degree) - long(std::abs(order)));
            const long double fac2 = hfactorial(long(degree) + long(std::abs(order)));
            
            return std::sqrt(static_cast<T>(bool(order == 0) ? 1. : 2.) * static_cast<T>(fac1 / fac2));
        }
        
    private:
        
        static inline long double hfactorial(long n)
        {
            long double result = n;
            if(n == 0)
                return 1;
            while(--n > 0)
                result *= n;
            
            return result;
        }
        
    private:
        
        size_t m_index = 0ul;
        size_t m_degree = 0ul;
        long   m_order = 0ul;
    };    
}
