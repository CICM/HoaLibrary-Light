/*
// Copyright (c) 2012-2015 Eliott Paris, Julien Colafrancesco, Thomas Le Meur & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_HARMONICS_LIGHT
#define DEF_HOA_HARMONICS_LIGHT

#include "Math.hpp"
#include "Signal.hpp"

namespace hoa
{
    //! The harmonic class owns basic harmonics ordering informations.
    /** The harmonic allows to retrieves informations about its ACN ordering, the degree and the order.
     */
    template <Dimension D, typename T> class Harmonic
    {
    public:

        //! The harmonic constructor.
        /** The harmonic constructor allocates and initializes the general member values depending on an index.
         @param index    The index must be at least 1.
         */
        Harmonic(const ulong index) noexcept;

        //! The harmonic destructor.
        /** The harmonic destructor free the memory.
         */
		~Harmonic() noexcept = 0;

        //! Get the index of the harmonic.
        /** The method returns the index \f$i\f$ of the harmonic.
         @return     The index.
         */
        ulong getIndex() const noexcept;

        //! Get the degree of the harmonic.
        /** The method returns the degree \f$l\f$ of the harmonic.
         @return     The degree.
         */
        ulong getDegree() const noexcept;

        //! Get the order of the harmonic.
        /** The method returns the order \f$m\f$ of the harmonic.
         @return     The order.
         */
        long getOrder() const noexcept;

        //! Get the name of the harmonic.
        /** The method returns the name \f$harmonic_{l,m}\f$ of the harmonic.
         @return     The name.
         */
        string getName() const noexcept;

        //! Get the normalization of the harmonic.
        /** The method returns the normalization of the harmonics.
         @return        The normalization of the harmonics.
         */
        T getNormalization() const noexcept;

        //! Get the semi-normalization of the harmonic.
        /** The method returns the semi-normalization of the harmonics.
         @return        The semi-normalization of the harmonics.
         */
        T getSemiNormalization() const noexcept;

        //! Get the index of an harmonic with its degree and its order.
        /** The method returns the index of the harmonic.
         @param degree  The degree of the harmonic.
         @param order   The order of the harmonic.
         @return        The index.
         */
        static ulong getHarmonicIndex(const ulong degree, const long order) noexcept;

        //! Get the degree of an harmonic with an index.
        /** The method returns the degree of the harmonic.
         @param index  The index of the harmonic.
         @return        The degree.
         */
        static ulong getHarmonicDegree(const ulong index) noexcept;

        //! Get the order of an harmonic with an index.
        /** The method returns the order of the harmonic.
         @param index  The index of the harmonic.
         @return        The degree.
         */
        static ulong getHarmonicOrder(const ulong index) noexcept;

        //! Get the number of harmonics for an order of decomposition.
        /** The method returns the number of harmonics for a order of decomposition \f$N\f$.
         @param order   The order of decomposition.
         @return        The number of harmonics.
         */
        static ulong getNumberOfHarmonics(const ulong order) noexcept;

        //! Get the normalization of an harmonic.
        /** The method returns the normalization of an harmonics.
         @param degree  The degree of the harmonic.
         @param order   The order of the harmonic.
         @return        The normalization of the harmonics.
         */
        static T getNormalization(const ulong degree, const long order) noexcept;

        //! Get the semi-normalization of an harmonic.
        /** The method returns the semi-normalization of an harmonics.
         @param degree  The degree of the harmonic.
         @param order   The order of the harmonic.
         @return        The semi-normalization of the harmonics.
         */
        static T getSemiNormalization(const ulong degree, const long order) noexcept;
    };

#ifndef DOXYGEN_SHOULD_SKIP_THIS

    template<typename T> class Harmonic<Hoa2d, T>
    {
    private:
        ulong m_index;
    public:

        //! The harmonic constructor.
        /** The harmonic constructor allocates and initializes the general member values depending on an index.
         @param index    The index must be at least 1.
         */
        Harmonic(const ulong index) noexcept :
        m_index(index)
        {
            ;
        }

        //! The harmonic destructor.
        /** The harmonic destructor free the memory.
         */
        ~Harmonic() noexcept
        {
            ;
        }

        //! Get the index of the harmonic.
        /** The method returns the index of the harmonic.
         @return     The index.
         */
        inline ulong getIndex() const noexcept
        {
            return m_index;
        }

        //! Get the degree of the harmonic.
        /** The method returns the degree of the harmonic.
         @return     The degree.
         */
        inline ulong getDegree() const noexcept
        {
            return (m_index + m_index % 2) * 0.5;
        }

        //! Get the order of the harmonic.
        /** The method returns the order of the harmonic.
         @return     The order.
         */
        inline long getOrder() const noexcept
        {
            return getDegree() * (1 - (m_index % 2) * 2);
        }

        //! Get the name of the harmonic.
        /** The method returns the name of the harmonic.
         @return     The name.
         */
        inline string getName() const noexcept
        {
            return "Harmonic " + to_string(getDegree()) + " " + to_string(getOrder());
        }

        //! Get the normalization of the harmonic.
        /** The method returns the normalization of the harmonics.
         @return        The normalization of the harmonics.
         */
        inline T getNormalization() const noexcept
        {
            return 1.;
        }

        //! Get the semi-normalization of the harmonic.
        /** The method returns the semi-normalization of the harmonics.
         @return        The semi-normalization of the harmonics.
         */
        inline T getSemiNormalization() const noexcept
        {
            return 1.;
        }

        //! Get the degree of an harmonic with an index.
        /** The method returns the degree of the harmonic.
         @param index  The index of the harmonic.
         @return        The degree.
         */
        static inline ulong getDegree(const ulong index) noexcept
        {
            return (index + index % 2) / 2ul;
        }

        //! Get the order of an harmonic with an index.
        /** The method returns the order of the harmonic.
         @param index  The index of the harmonic.
         @return        The order.
         */
        static inline long getOrder(const ulong index) noexcept
        {
            return long(long(index + index % 2l) / 2l) * (1l - (long)(index % 2) * 2l);
        }

        //! Get the index of an harmonic with its degree and its order.
        /** The method returns the index of the harmonic.
         @param degree  The degree of the harmonic.
         @param order   The order of the harmonic.
         @return        The index.
         */
        static inline long getIndex(const ulong degree, const long order) noexcept
        {
            return abs(order) *  2 - long(order < 0);
        }

        //! Get the number of harmonics for an order of decomposition.
        /** The method returns the number of harmonics for a order of decomposition.
         @param order   The order of decomposition.
         @return        The number of harmonics.
         */
        static inline ulong getNumberOfHarmonics(const ulong order) noexcept
        {
            return order * 2 + 1;
        }

        //! Get the normalization of an harmonic.
        /** The method returns the normalization of an harmonics.
         @param degree  The degree of the harmonic.
         @param order   The order of the harmonic.
         @return        The normalization of the harmonics.
         */
        static inline T getNormalization(const ulong degree, const long order) noexcept
        {
            return 1.;
        }

        //! Get the semi-normalization of an harmonic.
        /** The method returns the semi-normalization of an harmonics.
         @param degree  The degree of the harmonic.
         @param order   The order of the harmonic.
         @return        The semi-normalization of the harmonics.
         */
        static inline T getSemiNormalization(const ulong degree, const long order) noexcept
        {
            return 1.;
        }
    };

    template <typename T> class Harmonic<Hoa3d, T>
    {
    private:
        ulong m_index;
    public:

        //! The harmonic constructor.
        /** The harmonic constructor allocates and initializes the general member values depending on an index.
         @param index    The index must be at least 1.
         */
        Harmonic(const ulong index) noexcept :
        m_index(index)
        {

        }

        //! The harmonic destructor.
        /** The harmonic destructor free the memory.
         */
        ~Harmonic() noexcept
        {
            ;
        }

        //! Get the index of the harmonic.
        /** The method returns the index of the harmonic.
         @return     The index.
         */
        inline ulong getIndex() const noexcept
        {
            return m_index;
        }

        //! Get the degree of the harmonic.
        /** The method returns the degree of the harmonic.
         @return     The degree.
         */
        inline ulong getDegree() const noexcept
        {
            return ulong(sqrt(m_index));
        }

        //! Get the order of the harmonic.
        /** The method returns the order of the harmonic.
         @return     The order.
         */
        inline long getOrder() const noexcept
        {
            return long(m_index) - long(getDegree() * (getDegree() + 1));
        }

        //! Get the name of the harmonic.
        /** The method returns the name of the harmonic.
         @return     The name.
         */
        inline string getName() const noexcept
        {
            return "Harmonic " + to_string(getDegree()) + " " + to_string(getOrder());
        }

        //! Get the normalization of the harmonic.
        /** The method returns the normalization of the harmonics.
         @return        The normalization of the harmonics.
         */
        inline T getNormalization() const noexcept
        {
            return getNormalization(getDegree(), getOrder());
        }

        //! Get the semi-normalization of the harmonic.
        /** The method returns the semi-normalization of the harmonics.
         @return        The semi-normalization of the harmonics.
         */
        inline T getSemiNormalization() const noexcept
        {
            return getSemiNormalization(getDegree(), getOrder());
        }

        //! Get the degree of an harmonic with an index.
        /** The method returns the degree of the harmonic.
         @param index  The index of the harmonic.
         @return        The degree.
         */
        static inline ulong getDegree(const ulong index) noexcept
        {
            return ulong(sqrt(index));
        }

        //! Get the order of an harmonic with an index.
        /** The method returns the order of the harmonic.
         @param index  The index of the harmonic.
         @return        The order.
         */
        static inline long getOrder(const ulong index) noexcept
        {
            return long(index) - (long(sqrt(index)) * (long(sqrt(index)) + 1));
        }

        //! Get the index of an harmonic with its degree and its order.
        /** The method returns the index of the harmonic.
         @param degree  The degree of the harmonic.
         @param order   The order of the harmonic.
         @return        The index.
         */
        static inline ulong getIndex(const ulong degree, const long order) noexcept
        {
            return ulong(long(degree * (degree + 1)) + order);
        }

        //! Get the number of harmonics for an order of decomposition.
        /** The method returns the number of harmonics for a order of decomposition.
         @param order   The order of decomposition.
         @return        The number of harmonics.
         */
        static inline ulong getNumberOfHarmonics(const ulong order) noexcept
        {
            return (order + 1) * (order + 1);
        }

        //! Get the normalization of an harmonic.
        /** The method returns the normalization of an harmonics.
         @param degree  The degree of the harmonic.
         @param order   The order of the harmonic.
         @return        The normalization of the harmonics.
         */
        static inline T getNormalization(const ulong degree, const long order) noexcept
        {
            return getSemiNormalization(degree, order) * sqrt(2. * double(degree) + 1.);
        }

        //! Get the semi-normalization of an harmonic.
        /** The method returns the semi-normalization of an harmonics.
         @param degree  The degree of the harmonic.
         @param order   The order of the harmonic.
         @return        The semi-normalization of the harmonics.
         */
        static inline T getSemiNormalization(const ulong degree, const long order) noexcept
        {
            if(order == 0)
            {
                return T(sqrt(Math<T>::factorial(long(degree) - long(abs(order))) / Math<T>::factorial(long(order) + long(abs(order))))) * T(sqrt(T(1.) / T(4. * HOA_PI)));
            }
            else
            {
                return T(sqrt(Math<T>::factorial(long(degree) - abs(order)) / Math<T>::factorial(long(degree) + abs(order))) * sqrt(T(2.) / T(4. * HOA_PI)));
            }
        }
    };

#endif

}

#endif


