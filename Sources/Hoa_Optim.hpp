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

#ifndef DEF_HOA_OPTIM_LIGHT
#define DEF_HOA_OPTIM_LIGHT

#include "Hoa_Processor.hpp"

namespace hoa
{
    //! @brief The class optimizes the ambisonic sound field for several restitution systems.
    //! @details The class should be used to optimize the ambisonic sound field. There are 3
    //! optimizations, basic (no optimization), max-re (energy vector optimization) and
    //! in-phase (energy and velocity vector optimization).\n
    //! The basic optimization has no effect, it should be used (or not) with a perfect
    //! ambisonic channels arrangement where all the channels are to equal distance on a
    //! circle or a sphere, and for a listener placed at the perfect center of the circle of
    //! the sphere.\n
    //! The max-re should be used should be used for an auditory confined to the center of the
    //! circle or the sphere.\n
    //! The in-phase optimization should be used when the auditory covers the entire channels
    //! area and/or when the channels arrangement is not a perfect circle or a perfect sphere
    //! (when the channels are not to equal distance for example).\n
    //! Note that the optimizations decrease the precision of the sound field restitution thus
    //! it can be compared to particular cases of the fractional orders.
    template <Dimension D, typename T> class Optim : public Processor<D, T>::Harmonics
    {
    public:
        
        //! @brief The different optimization mode.
        //! @see getMode(), setMode()
        enum Mode
        {
            Basic = 0,  //!< The basic optimization
            MaxRe = 1,  //!< The max-re optimization
            InPhase = 2 //!< The in-phase optimization
        };

        //! @brief The constructor.
        //! @param order The order of decomposition.
        Optim(size_t order) hoa_noexcept : Processor<D, T>::Harmonics(order),
        m_weights(Signal<T>::alloc(Processor<D, T>::Harmonics::getNumberOfHarmonics()))
        { setMode(InPhase); }

        //! @brief The destructor.
		~Optim() hoa_noexcept { Signal<T>::free(m_weights); }
        
        //! @brief Returns the current optimization mode.
        inline Mode getMode() const hoa_noexcept { return m_mode; }
        
        //! @brief Set the optimization mode.
        //! @details The basic optimization generates a set of weights defined by:\n
        //!  \f[Y^{optimized}_{l,m} = Y_{l,m}\f]\n
        //! The max-re optimization generates a set of weights defined by:\n
        //! \f[Y^{optimized}_{l,m} = \cos{(l \times \frac{\pi}{2N + 2})} Y_{l,m} \f]\n
        //! The in-phase optimization generates a set of weights defined by:\n
        //! \f[Y^{optimized}_{l,m} = \frac{N!^2}{(N + l)!(N -l)!} Y_{l,m} \f]\n
        //! with \f$N\f$ the order of decomposition, \f$l\f$ the degree and \f$m\f$ the
        //! azimuthal order.
        //! @param mode The mode of optimization.
        void setMode(Mode mode) hoa_noexcept
        {
            const size_t size   = Processor<D, T>::Harmonics::getNumberOfHarmonics();
            const size_t order  = Processor<D, T>::Harmonics::getDecompositionOrder();
            if(mode == Basic)
            {
                for(size_t i = 0; i < size;  ++i) {
                    m_weights[i] = 1.;
                }
            }
            else if(mode == MaxRe)
            {
                m_weights[0] = T(1.);
                for(size_t i = 1; i < size; i++) {
                    const size_t degree = Processor<D, T>::Harmonics::getHarmonicDegree(i);
                    m_weights[i] = cos(T(degree) *  T(HOA_PI) / T(2. * order + 2.));
                }
            }
            else
            {
                m_weights[0] = T(1.);
                const T facn = Math<T>::factorial(long(order));
                for(size_t i = 1; i < size; i++) {
                    const size_t degree = Processor<D, T>::Harmonics::getHarmonicDegree(i);
                    m_weights[i] = facn / Math<T>::factorial(long(order - degree)) * facn / Math<T>::factorial(long(order + degree));
                }
                
            }
            m_mode = mode;
        }

        //! @brief The method performs the optimization on the harmonics signal.
        //! @details The method can be used for in-place or not-in-place processing and sample
        //! by sample. The inputs array and outputs array contains the spherical harmonics
        //! samples thus the minimum size of the array must be the number of harmonics.
        //! @param inputs  The inputs array.
        //! @param outputs The outputs array.
        void process(T const* inputs, T* outputs) hoa_noexcept hoa_final
        {
            const  size_t size = Processor<D, T>::Harmonics::getNumberOfHarmonics();
            const T* weights = m_weights;
            for(size_t i = size>>3; i; --i, inputs += 8, weights += 8, outputs += 8)
            {
                T const f0 = inputs[0], f1 = inputs[1], f2 = inputs[2], f3 = inputs[3];
                T const f4 = inputs[4], f5 = inputs[5], f6 = inputs[6], f7 = inputs[7];
                
                T const g0 = weights[0], g1 = weights[1], g2 = weights[2], g3 = weights[3];
                T const g4 = weights[4], g5 = weights[5], g6 = weights[6], g7 = weights[7];
                
                outputs[0] = f0 * g0; outputs[1] = f1 * g1; outputs[2] = f2 * g2; outputs[3] = f3 * g3;
                outputs[4] = f4 * g4; outputs[5] = f5 * g5; outputs[6] = f6 * g6; outputs[7] = f7 * g7;
            }
            for(size_t i = size&7; i; --i, inputs++, weights++, outputs++)
            {
                outputs[0] = inputs[0] * weights[0];
            }
        }
    protected:
        Mode m_mode;
        T*   m_weights;
    };
}

#endif
