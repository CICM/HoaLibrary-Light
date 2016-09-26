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
    //! The optim class optimizes the ambisonic sound field for several restitution systems.
    /** The optim should be used to optimize the ambisonic sound field. There are 3 optimizations, Basic (no optimization), MaxRe (energy vector optimization) and InPhase (energy and velocity vector optimization). Basic has no effect, it should be used (or not) with a perfect ambisonic channels arrangement where all the channels are to equal distance on a circle or a sphere, and for a listener placed at the perfect center of the circle of the sphere. MaxRe should be used should be used for an auditory confined to the center of the circle of the sphere. InPhase should be used when the auditory covers the entire channels area and when the channels arrangement is not a perfect circle or a perfect sphere or when the channels are not to equal distance. Note that the optimizations decrease the precision of the sound field restitution thus it can be compared to particular cases of the fractional orders.
     */
    template <Dimension D, typename T> class Optim : public Processor<D, T>::Harmonics
    {
    protected:
        T*  m_weights;
    public:

        //! The optim constructor.
        /**	The optim constructor allocates and initialize the member values to computes spherical harmonics weighted coefficients depending on a order of decomposition. The order must be at least 1.
         @param     order	The order.
         */
        Optim(const size_t order) hoa_noexcept : Processor<D, T>::Harmonics(order),
        m_weights(Signal<T>::alloc(Processor<D, T>::Harmonics::getNumberOfHarmonics())) {}

        //! The optim destructor.
        /**	The optim destructor free the memory.
         */
		virtual ~Optim() hoa_noexcept hoa_default_f

        //! This method performs the optimization.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        virtual void process(T const* inputs, T* outputs) hoa_noexcept
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
    };
    
    //! The basic optim.
    /** The basic optim has no effect, it should be used (or not) with a perfect ambisonic channels arrangement where all the channels are to equal distance on a circle or a sphere, and for a listener placed at the perfect center of the circle of the sphere. This method performs the basic optimization.
        \f[Y^{optimized}_{l,m} = Y_{l,m}\f]
        with \f$l\f$ the degree and \f$m\f$ the order.
     */
    template <Dimension D, typename T> class OptimBasic : public Optim<D, T>
    {
    public:
        
        //! The optim constructor.
        /**	The optim constructor allocates and initialize the member values to computes spherical harmonics weighted coefficients depending on a order of decomposition. The order must be at least 1.
         @param     order	The order.
         */
        OptimBasic(const size_t order) hoa_noexcept : Optim<D, T>(order)
        {
            const  size_t size = Processor<D, T>::Harmonics::getNumberOfHarmonics();
            T* weights = Optim<D, T>::m_weights;
            for(size_t i = 0; i < size;  ++i)
            {
                weights[i] = 1.;
            }
        }
        
        //! The optim destructor.
        /**	The optim destructor free the memory.
         */
        ~OptimBasic() hoa_noexcept hoa_default_f
    };
    
    //! The maxre optim.
    /** The maxre optim should be used for an auditory confined to the center of the circle of the sphere. \f[Y^{optimized}_{l,m} = \cos{(l \times \frac{\pi}{2N + 2})} Y_{l,m} \f]
     with \f$N\f$ the order of decomposition, \f$l\f$ the degree and \f$m\f$ the order.
     */
    template <Dimension D, typename T> class OptimMaxRe : public Optim<D, T>
    {
    public:
        //! The optim constructor.
        /**	The optim constructor allocates and initialize the member values to computes spherical harmonics weighted coefficients depending on a order of decomposition. The order must be at least 1. 
         @param     order	The order.
         */
        OptimMaxRe(const size_t order) hoa_noexcept : Optim<D, T>(order)
        {
            const  size_t size = Processor<D, T>::Harmonics::getNumberOfHarmonics();
            T* weights = Optim<D, T>::m_weights;
            weights[0] = T(1.);
            for(size_t i = 1; i < size; i++)
            {
                const size_t degree = Processor<D, T>::Harmonics::getHarmonicDegree(i);
                weights[i] = cos(T(degree) *  T(HOA_PI) / T(2. * order + 2.));
            }
        }
        
        //! The optim destructor.
        /**	The optim destructor free the memory.
         */
        ~OptimMaxRe() hoa_noexcept hoa_default_f
    };
    
    //! The inphase optim.
    /** The inphase optim should be used when the auditory covers the entire channels area and when the channels arrangement is not a perfect circle or a perfect sphere or when the channels are not to equal distance. \f[Y^{optimized}_{l,m} = \frac{N!^2}{(N + l)!(N -l)!} Y_{l,m} \f]
     with \f$N\f$ the order of decomposition, \f$l\f$ the degree and \f$m\f$ the order.
     */
    template <Dimension D, typename T> class OptimInPhase : public Optim<D, T>
    {
    public:
        //! The optim constructor.
        /**	The optim constructor allocates and initialize the member values to computes spherical harmonics weighted coefficients depending on a order of decomposition. The order must be at least 1.
         @param     order	The order.
         */
        OptimInPhase(const size_t order) hoa_noexcept : Optim<D, T>(order)
        {
            const  size_t size = Processor<D, T>::Harmonics::getNumberOfHarmonics();
            T* weights = Optim<D, T>::m_weights;
            weights[0] = T(1.);
            const T facn = Math<T>::factorial(long(order));
            for(size_t i = 1; i < size; i++)
            {
                const size_t degree = Processor<D, T>::Harmonics::getHarmonicDegree(i);
                weights[i] = facn / Math<T>::factorial(long(order - degree)) * facn / Math<T>::factorial(long(order + degree));
            }
        }
        
        //! The optim destructor.
        /**	The optim destructor free the memory.
         */
        ~OptimInPhase() hoa_noexcept hoa_default_f;
    };

}

#endif
