/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_OPTIM_LIGHT
#define DEF_HOA_OPTIM_LIGHT

#include "Harmonics.hpp"

namespace hoa
{
    enum Optimization
    {
        Basic = 0,
        MaxRe = 1,
        InPhase=2
    };
    
    //! The ambisonic optimization.
    /** The optimization should be used to optimize the ambisonic sound field. There are 3 optimization modes, Basic (no optimization), MaxRe (energy vector optimization) and InPhase (energy and velocity vector optimization). Basic has no effect, it should be used with a perfect ambisonic channels, arrengement where all the channels are to equal distance on a circle, and for a listener placed at the perfect center of the circle. MaxRe should be used for auditory confined to the center of the circle. InPhase should be used when the auditory covers the entire channels area and when the channels arragement is not a perfect circle or when the channels are not to equal distance. Note that the optimizations decrease the precision of the sound field restitution thus it can be compared to particular cases of the fractional orders.
     */
    template <Dimension D, typename T, Optimization O> class Optim;
    
    //! The ambisonic optimization.
    /** The basic.
     */
    template <typename T> class Optim<Hoa2d, T, Basic> : public Harmonic<Hoa2d, T>::Processor
    {
    public:
        
        //! The optimization constructor.
        /**	The optimization constructor allocates and initialize the member values to computes spherical harmonics weighted coefficients depending of a decomposition order. The order must be at least 1.
         @param     order	The order.
         @param     optim	The optimization.
         */
        Optim(const ulong order) noexcept :
        Harmonic<Hoa2d, T>::Processor(order)
        {
            ;
        }
        
        //! The optimization destructor.
        /**	The optimization destructor free the memory.
         */
        ~Optim()
        {
            ;
        }
        
        //! This method performs the basic optimization.
        /**	You should use this method for in-place or not-in-place processing and performs the basic optimization sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        inline void processInphase(T const* inputs, T* outputs) const noexcept
        {
            (*outputs++)  = (*inputs++);
            (*outputs++)  = (*inputs++);
            (*outputs++)  = (*inputs++);
            for(ulong i = 2; i <= Harmonic<Hoa2d, T>::Processor::getDecompositionOrder(); i++)
            {
                (*outputs++)  = (*inputs++);
                (*outputs++)  = (*inputs++);
            }
        }
    };
    
    //! The ambisonic optimization.
    /** The optimization should be used to optimize the ambisonic sound field. There are 3 optimization modes, Basic (no optimization), MaxRe (energy vector optimization) and InPhase (energy and velocity vector optimization). Basic has no effect, it should be used with a perfect ambisonic channels, arrengement where all the channels are to equal distance on a circle, and for a listener placed at the perfect center of the circle. MaxRe should be used for auditory confined to the center of the circle. InPhase should be used when the auditory covers the entire channels area and when the channels arragement is not a perfect circle or when the channels are not to equal distance. Note that the optimizations decrease the precision of the sound field restitution thus it can be compared to particular cases of the fractional orders.
     */
    template <typename T> class Optim<Hoa2d, T, MaxRe> : public Harmonic<Hoa2d, T>::Processor
    {
    private:
        const T   m_cosmaxRe;
    public:
        
        //! The optimization constructor.
        /**	The optimization constructor allocates and initialize the member values to computes spherical harmonics weighted coefficients depending of a decomposition order. The order must be at least 1.
         @param     order	The order.
         @param     optim	The optimization.
         */
        Optim(const ulong order) noexcept :
        Harmonic<Hoa2d, T>::Processor(order),
        m_cosmaxRe(cos(HOA_PI / (T)(2. * order + 2.)))
        {
            ;
        }
        
        //! The optimization destructor.
        /**	The optimization destructor free the memory.
         */
        ~Optim()
        {
            ;
        }
        
        //! This method performs the max-re optimization.
        /**	You should use this method for in-place or not-in-place processing and performs the max-re optimization sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        inline void process(T const* inputs, T* outputs) const noexcept
        {
            T cos_re = m_cosmaxRe;
            (*outputs++)  = (*inputs++);
            (*outputs++)  = (*inputs++) * cos_re;
            (*outputs++)  = (*inputs++) * cos_re;
            for(ulong i = 2; i <= Harmonic<Hoa2d, T>::Processor::getDecompositionOrder(); i++)
            {
                cos_re = 2. * cos_re * cos_re - 1.;
                (*outputs++)  = (*inputs++) * cos_re;
                (*outputs++)  = (*inputs++) * cos_re;
            }
        }
    };
    
    //! The ambisonic optimization.
    /** The basic.
     */
    template <typename T> class Optim<Hoa2d, T, InPhase> : public Harmonic<Hoa2d, T>::Processor
    {
    private:
        const T   m_facorder;
        const T   m_facinphase;
        const T   m_facorder1;
        const T   m_facorder2;
    public:
        
        //! The optimization constructor.
        /**	The optimization constructor allocates and initialize the member values to computes spherical harmonics weighted coefficients depending of a decomposition order. The order must be at least 1.
         @param     order	The order.
         @param     optim	The optimization.
         */
        Optim(const ulong order) noexcept :
        Harmonic<Hoa2d, T>::Processor(order),
        m_facorder(factorial(order)),
        m_facinphase(m_facorder * m_facorder * order),
        m_facorder1(m_facorder * (order + 2) * order),
        m_facorder2(m_facorder / order)
        {
            ;
        }
        
        //! The optimization destructor.
        /**	The optimization destructor free the memory.
         */
        ~Optim()
        {
            ;
        }
        
        //! This method performs the in-phase optimization.
        /**	You should use this method for in-place or not-in-place processing and performs the in-phase optimization sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        inline void process(T const* inputs, T* outputs) const noexcept
        {
            T order1 = Harmonic<Hoa2d, T>::Processor::getDecompositionOrder() + 3;
            T order2 = Harmonic<Hoa2d, T>::Processor::getDecompositionOrder() - 1;
            T factor1 = m_facorder1;
            T factor2 = m_facorder2;
            T factor  = m_facinphase / (factor1 * factor2);
            
            (*outputs++)  = (*inputs++);
            (*outputs++)  = (*inputs++) * factor;
            (*outputs++)  = (*inputs++) * factor;
            for(ulong i = 2; i <= Harmonic<Hoa2d, T>::Processor::getDecompositionOrder(); i++)
            {
                factor1 *= order1++;
                factor2 /= order2--;
                factor = m_facinphase / (factor1 * factor2);
                
                (*outputs++)  = (*inputs++) * factor;
                (*outputs++)  = (*inputs++) * factor;
            }
        }
    };
}

#endif



