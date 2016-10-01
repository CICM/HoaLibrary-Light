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

#ifndef DEF_HOA_HARMONICS_LIGHT
#define DEF_HOA_HARMONICS_LIGHT

#include "Hoa_Math.hpp"
#include "Hoa_Signal.hpp"

namespace hoa
{
    //! @brief The class owns basic harmonics informations.
    //! @details The class allows to retrieves several informations about the harmonics: the
    //! numbering (ACN), the normalization (SN3D/N3D), the degree, the aziumthal order, etc.
    template <Dimension D, typename T> class Harmonic
    {
    public:

        //! @brief The harmonic constructor.
        //! @param index    The index must be at least 1.
        Harmonic(const size_t index) hoa_noexcept;

        //! @brief The harmonic destructor.
        /** The harmonic destructor free the memory.
         */
		~Harmonic() hoa_noexcept = 0;

        //! @brief Returns the index of the harmonic.
        size_t getIndex() const hoa_noexcept;

        //! @brief Returns the degree of the harmonic.
        size_t getDegree() const hoa_noexcept;

        //! @brief Returns the azimuthal order of the harmonic.
        long getOrder() const hoa_noexcept;

        //! @brief Returns the name of the harmonic.
        std::string getName() const hoa_noexcept;

        //! @brief Returns the normalization of the harmonic.
        T getNormalization() const hoa_noexcept;

        //! @brief Returns the semi-normalization of the harmonic.
        T getSemiNormalization() const hoa_noexcept;

        //! @brief Returns the index of an harmonic for a specific degree and azimuthal order.
        //! @param degree  The degree of the harmonic.
        //! @param order   The order of the harmonic.
        static size_t getHarmonicIndex(const size_t degree, const long order) hoa_noexcept;

        //! @brief Returns the degree of an harmonic for an index.
        //! @param index  The index of the harmonic.
        static size_t getHarmonicDegree(const size_t index) hoa_noexcept;

        //! @brief Returns the azimuthal order of an harmonic for an index.
        //! @param index  The index of the harmonic.
        static size_t getHarmonicOrder(const size_t index) hoa_noexcept;

        //! @brief Returns the number of harmonics for an order of decomposition.
        //! @param order   The order of decomposition.
        static size_t getNumberOfHarmonics(const size_t order) hoa_noexcept;
        
        //! @brief Returns the number of harmonics for a degree.
        //! @param degree  The degree.
        static inline size_t getNumberOfHarmonicsInDegree(const size_t degree) hoa_noexcept;

        //! @brief Returns the normalization of an harmonic.
        //! @param degree  The degree of the harmonic.
        //! @param order   The order of the harmonic.
        static T getNormalization(const size_t degree, const long order) hoa_noexcept;

        //! @brief Returns the semi-normalization of an harmonic.
        //! @param degree  The degree of the harmonic.
        //! @param order   The order of the harmonic.
        static T getSemiNormalization(const size_t degree, const long order) hoa_noexcept;
    };
    
#ifndef DOXYGEN_SHOULD_SKIP_THIS

    //! @brief The class owns basic harmonics informations.
    //! @details The class allows to retrieves several informations about the harmonics: the
    //! numbering (ACN), the normalization (SN3D/N3D), the degree, the aziumthal order, etc.
    template<typename T> class Harmonic<Hoa2d, T>
    {
    private:
        size_t m_index;
        size_t m_degree;
        long   m_order;
    public:

        //! @brief The harmonic constructor.
        //! @param index    The index must be at least 1.
        inline Harmonic(const size_t index) hoa_noexcept :
        m_index(index), m_degree(getDegree(index)), m_order(getOrder(index)) {}

        //! @brief The harmonic destructor.
        inline ~Harmonic() hoa_noexcept hoa_default_f

        //! @brief Returns the index of the harmonic.
        inline size_t getIndex() const hoa_noexcept { return m_index; }

        //! @brief Returns the degree of the harmonic.
        inline size_t getDegree() const hoa_noexcept { return m_degree; }

        //! @brief Returns the azimuthal order of the harmonic.
        inline long getOrder() const hoa_noexcept { return m_order; }
        
        //! @brief Returns the name of the harmonic.
        std::string getName() const hoa_noexcept {
            std::ostringstream ostr;
            ostr <<  "Harmonic " << m_degree << " " << m_order;
            return ostr.str();
        }

        //! @brief Returns the normalization of the harmonic.
        inline T getNormalization() const hoa_noexcept { return getNormalization(m_degree, m_order); }

        //! @brief Returns the semi-normalization of the harmonic.
        inline T getSemiNormalization() const hoa_noexcept { return getSemiNormalization(m_degree, m_order); }

        //! @brief Returns the index of an harmonic for a specific degree and azimuthal order.
        //! @param degree  The degree of the harmonic.
        //! @param order   The order of the harmonic.
        static inline long getIndex(const size_t degree, const long order) hoa_noexcept {
            return std::abs(order) *  2 - long(order < 0); }
        
        //! @brief Returns the degree of an harmonic for an index.
        //! @param index  The index of the harmonic.
        static inline size_t getDegree(const size_t index) hoa_noexcept {
            return (index + index % 2) / 2ul; }

        //! @brief Returns the azimuthal order of an harmonic for an index.
        //! @param index  The index of the harmonic.
        static inline long getOrder(const size_t index) hoa_noexcept {
            return long(long(index + index % 2l) / 2l) * (1l - (long)(index % 2) * 2l); }

        //! @brief Returns the number of harmonics for an order of decomposition.
        //! @param order   The order of decomposition.
        static inline size_t getNumberOfHarmonics(const size_t order) hoa_noexcept { return order * 2 + 1; }
        
        //! @brief Returns the number of harmonics for a degree.
        //! @param degree  The degree.
        static inline size_t getNumberOfHarmonicsInDegree(const size_t degree) hoa_noexcept { return (degree != 0) + 1; }

        //! @brief Returns the normalization of an harmonic.
        //! @param degree  The degree of the harmonic.
        //! @param order   The order of the harmonic.
        static inline T getNormalization(const size_t degree, const long order) hoa_noexcept { return 1.; }

        //! @brief Returns the semi-normalization of an harmonic.
        //! @param degree  The degree of the harmonic.
        //! @param order   The order of the harmonic.
        static inline T getSemiNormalization(const size_t degree, const long order) hoa_noexcept { return 1.; }
    };

    //! @brief The class owns basic harmonics informations.
    //! @details The class allows to retrieves several informations about the harmonics: the
    //! numbering (ACN), the normalization (SN3D/N3D), the degree, the aziumthal order, etc.
    template <typename T> class Harmonic<Hoa3d, T>
    {
    private:
        size_t m_index;
        size_t m_degree;
        long   m_order;
    public:

        //! @brief The harmonic constructor.
        //! @param index    The index must be at least 1.
        Harmonic(const size_t index) hoa_noexcept :
        m_index(index), m_degree(getDegree(index)), m_order(getOrder(index)) {}

        //! @brief The harmonic destructor.
        ~Harmonic() hoa_noexcept hoa_default_f

        //! @brief Returns the index of the harmonic.
        inline size_t getIndex() const hoa_noexcept { return m_index; }

        //! @brief Returns the degree of the harmonic.
        inline size_t getDegree() const hoa_noexcept {return m_degree; }

        //! @brief Returns the azimuthal order of the harmonic.
        inline long getOrder() const hoa_noexcept { return m_order; }

        //! @brief Returns the name of the harmonic.
        std::string getName() const hoa_noexcept {
            std::ostringstream ostr;
            ostr <<  "Harmonic " << m_degree << " " << m_order;
            return ostr.str();
        }

        //! @brief Returns the normalization of the harmonic.
        inline T getNormalization() const hoa_noexcept { return getNormalization(m_degree, m_order); }

        //! @brief Returns the semi-normalization of the harmonic.
        inline T getSemiNormalization() const hoa_noexcept { return getSemiNormalization(m_degree, m_order); }
        
        //! @brief Returns the degree of an harmonic for an index.
        //! @param index  The index of the harmonic.
        static inline size_t getDegree(const size_t index) hoa_noexcept {
            return size_t(sqrt(double(index))); }

        //! @brief Returns the azimuthal order of an harmonic for an index.
        //! @param index  The index of the harmonic.
        static inline long getOrder(const size_t index) hoa_noexcept {
            return long(index) - (long(sqrt(double(index))) * (long(sqrt(double(index))) + 1)); }

        //! @brief Returns the index of an harmonic for a specific degree and azimuthal order.
        //! @param degree  The degree of the harmonic.
        //! @param order   The order of the harmonic.
        static inline size_t getIndex(const size_t degree, const long order) hoa_noexcept {
            return size_t(long(degree * (degree + 1)) + order); }

        //! @brief Returns the number of harmonics for an order of decomposition.
        //! @param order   The order of decomposition.
        static inline size_t getNumberOfHarmonics(const size_t order) hoa_noexcept  {
            return (order + 1) * (order + 1); }

        //! @brief Returns the number of harmonics for a degree.
        //! @param degree  The degree.
        static inline size_t getNumberOfHarmonicsInDegree(const size_t degree) hoa_noexcept {
            return degree * 2 + 1; }
        
        //! @brief Returns the normalization of an harmonic.
        //! @param degree  The degree of the harmonic.
        //! @param order   The order of the harmonic.
        static inline T getNormalization(const size_t degree, const long order) hoa_noexcept
        {
            return getSemiNormalization(degree, order) * T(sqrt(T(2.) * T(degree) + T(1.)));
        }

        //! @brief Returns the semi-normalization of an harmonic.
        //! @param degree  The degree of the harmonic.
        //! @param order   The order of the harmonic.
        static inline T getSemiNormalization(const size_t degree, const long order) hoa_noexcept
        {
            const long double fac1 = hfactorial(long(degree) - long(std::abs(order)));
            const long double fac2 = hfactorial(long(degree) + long(std::abs(order)));
            //return T(sqrt(fac1 / fac2)) * (bool(order == 0) ? T(0.2820947918) : T(0.3989422804));
            //return T(sqrt(fac1 / fac2)) * (bool(order == 0) ? T(1.) : T(1.41421356237309504880168872420969808));
            return T(sqrt((bool(order == 0) ? T(1.) : T(2.)) * T(fac1 / fac2)));
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
    };
#endif
    
}

#endif


