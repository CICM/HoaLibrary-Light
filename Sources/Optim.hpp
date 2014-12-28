/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_OPTIM
#define DEF_HOA_OPTIM

#include "Harmonics.hpp"

namespace hoa
{
    //! The ambisonic optimization.
    /** The optimization should be used to optimize the ambisonic sound field. There are 3 optimization modes, Basic (no optimization), MaxRe (energy vector optimization) and InPhase (energy and velocity vector optimization). Basic has no effect, it should be used with a perfect ambisonic channels, arrengement where all the channels are to equal distance on a circle, and for a listener placed at the perfect center of the circle. MaxRe should be used for auditory confined to the center of the circle. InPhase should be used when the auditory covers the entire channels area and when the channels arragement is not a perfect circle or when the channels are not to equal distance. Note that the optimizations decrease the precision of the sound field restitution thus it can be compared to particular cases of the fractional orders.
     */
    template <Dimension D, typename T> class Optim : public Harmonic<D>::Processor
    {
    private:
        Optimization    m_mode;
        T*              m_gains;
    public:
        
        //! The optimization constructor.
        /**	The optimization constructor allocates and initialize the member values to computes spherical harmonics weighted coefficients depending of a decomposition order. The order must be at least 1.
         @param     order	The order.
         @param     optim	The optimization.
         */
        Optim(const ulong order, const Optimization optim = InPhase) noexcept :
        Harmonic<D>::Processor(order)
        {
            m_gains = new T[Harmonic<D>::Processor::getDecompositionOrder()];
            setOptimization(optim);
        }
        
        //! The optimization destructor.
        /**	The optimization destructor free the memory.
         */
        ~Optim()
        {
            delete [] m_gains;
        }
        
        //! Retrieve the optimization mode.
        /** Retrieve the optimization mode : Basic, MaxRe or InPhase.
         @return The method returns the optimization mode.
         */
        inline Optimization getOptimization() const noexcept
        {
            return m_mode;
        }
        
        //! This method set the optimization mode.
        /**	The mode should be one of the 3 optimization modes, Basic, MaxRe or InPhase.
         @param     mode The optimization mode.
         */
        inline void setOptimization(const Optimization mode) noexcept
        {
            m_mode = mode;
            if(m_mode == Optimization::Basic)
            {
                for(ulong i = 0; i < Harmonic<D>::Processor::getDecompositionOrder(); i++)
                {
                    m_gains[i] = 1.;
                }
            }
            else if(m_mode == Optimization::MaxRe)
            {
                for(ulong i = 0; i < Harmonic<D>::Processor::getDecompositionOrder(); i++)
                {
                    m_gains[i] = cos((T)i * HOA_PI / (T)(2. * Harmonic<D>::Processor::getDecompositionOrder() + 2.));
                }
            }
            else
            {
                for(ulong i = 0; i < Harmonic<D>::Processor::getDecompositionOrder(); i++)
                {
                    const long double temp1 = factorial(Harmonic<D>::Processor::getDecompositionOrder()) / factorial(Harmonic<D>::Processor::getDecompositionOrder() + i + 1.);
                    const long double temp2 = factorial(Harmonic<D>::Processor::getDecompositionOrder() + 1.) / factorial(Harmonic<D>::Processor::getDecompositionOrder() - i);
                    m_gains[i] = temp1 * temp2;
                }
            }
        }
        
        //! This method performs the optimization.
        /**	You should use this method for in-place or not-in-place processing and performs the optimization sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        inline void process(T const* inputs, T* outputs) const noexcept
        {
            T const* factor = m_gains;
            *(outputs++) = (*inputs++) * (*factor++);
            *(outputs++) = (*inputs++) * (*factor);
            *(outputs++) = (*inputs++) * (*factor++);
            for(ulong i = 2; i < Harmonic<D>::Processor::getDecompositionOrder(); i++)
            {
                *(outputs++) = (*inputs++) * (*factor);
                *(outputs++) = (*inputs++) * (*factor++);
            }
        }
    };
}

#endif



