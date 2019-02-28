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

namespace hoa
{
    // ================================================================================ //
    // WIDER //
    // ================================================================================ //
    
    //! @brief The class widens the propagation of the sounds in a sound field.
    //! @details The class simulates fractional orders of decomposition to reduce the
    //! precision of the sound field. When the factor of widening is \f$0\f$ sound field, only
    //! the first hamonic \f$Y_{0,0}\f$ remains and the sound field is omni directional. By
    //! increasing the factor of widening toward \f$1\f$, the other harmonics appears in a
    //! logarithmic way, increasing the precision of the sound field that becomes more and
    //! more directional until all the harmonics appeared. The weight of the harmonics are
    //! defined by:
    //! \f[W_{l,m}(x) = x^l((1-x)(N-l)+1)\f]<br>
    //! with \f$N\f$ the order of decomposition, \f$l\f$ the degree, \f$m\f$ the
    //! azimuthal order and \f$x\f$ the factor of widening.
    template <Dimension D, typename T>
    class Wider
    : public ProcessorHarmonics<D, T>
    {
    public:
        
        //! @brief Constructor.
        //! @param order The order of decomposition.
        Wider(const size_t order)
        : ProcessorHarmonics<D, T>(order)
        , m_coeffs(order + 1)
        {
            setWidening(1.);
        }

        //! @brief Destructor.
        ~Wider() = default;

        //! @brief This method set factor of widening.
        //! @param value The factor of widening.
        inline void setWidening(const T value) noexcept
        {
            m_widening = std::max(std::min(value, T(1.)), T(0.));
            const size_t order  = ProcessorHarmonics<D, T>::getDecompositionOrder();
            const T      temp   = T(1) - m_widening;
            T* coeff            = m_coeffs.data();
            
            (*coeff++) = (T(order) * temp + T(1.));
            for(size_t i = 1; i <= order; ++i)
            {
                (*coeff++) = std::pow(m_widening, T(i)) * (temp * T(order - i) + T(1.));
            }
        }

        //! @brief Returns the the widening value.
		inline T getWidening() const noexcept { return m_widening; }

        //! @brief The method performs the widening on the harmonics signal.
        //! @details The method can be used for in-place or not-in-place processing and sample
        //! by sample. The inputs array and outputs array contains the spherical harmonics
        //! samples thus the minimum size of the array must be the number of harmonics.
        //! @param inputs  The inputs array.
        //! @param outputs The outputs array.
		void process(const T* inputs, T* outputs) noexcept final
        {
            const size_t order  = ProcessorHarmonics<D, T>::getDecompositionOrder();
            T* coeff     = m_coeffs.data();
            (*outputs++) = (*inputs++) * (*coeff++);
            for(size_t i = 1; i <= order; ++i)
            {
                const size_t nharmo = Harmonic<D, T>::getNumberOfHarmonicsInDegree(i);
                const T      factor = (*coeff++);
                for(size_t j = 0; j < nharmo; ++j)
                {
                    (*outputs++) = (*inputs++) * factor;
                }
            }
        }
        
    private:
        
        T m_widening = 0.;
        std::vector<T> m_coeffs {};
    };
}
