/*
// Copyright (c) 2012-2015 Pierre Guillot, Eliott Paris & Thomas Le Meur CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_OPTIM_LIGHT
#define DEF_HOA_OPTIM_LIGHT

#include "Processor.hpp"

namespace hoa
{
    //! The optim class optimizes the ambisonic sound field for several restitution systems.
    /** The optim should be used to optimize the ambisonic sound field. There are 3 optimizations, Basic (no optimization), MaxRe (energy vector optimization) and InPhase (energy and velocity vector optimization). Basic has no effect, it should be used (or not) with a perfect ambisonic channels arrangement where all the channels are to equal distance on a circle or a sphere, and for a listener placed at the perfect center of the circle of the sphere. MaxRe should be used should be used for an auditory confined to the center of the circle of the sphere. InPhase should be used when the auditory covers the entire channels area and when the channels arrangement is not a perfect circle or a perfect sphere or when the channels are not to equal distance. Note that the optimizations decrease the precision of the sound field restitution thus it can be compared to particular cases of the fractional orders.
     */
    template <Dimension D, typename T> class Optim : public Processor<D, T>::Harmonics
    {
    public:

        //! The optim constructor.
        /**	The optim constructor allocates and initialize the member values to computes spherical harmonics weighted coefficients depending on a order of decomposition. The order must be at least 1.
         @param     order	The order.
         */
        Optim(const size_t order) hoa_noexcept;

        //! The optim destructor.
        /**	The optim destructor free the memory.
         */
		virtual ~Optim() hoa_noexcept = 0;

        //! This method performs the optimization.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        virtual void process(T const* inputs, T* outputs) hoa_noexcept;
    };
    
    //! The basic optim.
    /** The basic optim has no effect, it should be used (or not) with a perfect ambisonic channels arrangement where all the channels are to equal distance on a circle or a sphere, and for a listener placed at the perfect center of the circle of the sphere.
     */
    template <Dimension D, typename T> class OptimBasic : public Optim<D, T>
    {
    public:
        
        //! The optim constructor.
        /**	The optim constructor allocates and initialize the member values to computes spherical harmonics weighted coefficients depending on a order of decomposition. The order must be at least 1.
         @param     order	The order.
         */
        OptimBasic(const size_t order) hoa_noexcept;
        
        //! The optim destructor.
        /**	The optim destructor free the memory.
         */
        virtual ~OptimBasic() hoa_noexcept = 0;
        
        //! This method performs the basic optimization.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         \f[Y^{optimized}_{l,m} = Y_{l,m}\f]
         with \f$l\f$ the degree and \f$m\f$ the order.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        virtual void process(T const* inputs, T* outputs) hoa_noexcept;
    };
    
    //! The maxre optim.
    /** The maxre optim should be used for an auditory confined to the center of the circle of the sphere.
     */
    template <Dimension D, typename T> class OptimMaxRe : public Optim<D, T>
    {
    public:
        //! The optim constructor.
        /**	The optim constructor allocates and initialize the member values to computes spherical harmonics weighted coefficients depending on a order of decomposition. The order must be at least 1.
         @param     order	The order.
         */
        OptimMaxRe(const size_t order) hoa_noexcept;
        
        //! The optim destructor.
        /**	The optim destructor free the memory.
         */
        virtual ~OptimMaxRe() hoa_noexcept = 0;
        
        //! This method performs the max-re optimization.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         \f[Y^{optimized}_{l,m} = \cos{(l \times \frac{\pi}{2N + 2})} Y_{l,m} \f]
         with \f$N\f$ the order of decomposition, \f$l\f$ the degree and \f$m\f$ the order.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        virtual void process(T const* inputs, T* outputs) hoa_noexcept;
    };
    
    //! The inphase optim.
    /** The inphase optim should be used when the auditory covers the entire channels area and when the channels arrangement is not a perfect circle or a perfect sphere or when the channels are not to equal distance.
     */
    template <Dimension D, typename T> class OptimInPhase : public Optim<D, T>
    {
    public:
        //! The optim constructor.
        /**	The optim constructor allocates and initialize the member values to computes spherical harmonics weighted coefficients depending on a order of decomposition. The order must be at least 1.
         @param     order	The order.
         */
        OptimInPhase(const size_t order) hoa_noexcept;
        
        //! The optim destructor.
        /**	The optim destructor free the memory.
         */
        virtual ~OptimInPhase() hoa_noexcept = 0;
        
        //! This method performs the in-phase optimization.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         \f[Y^{optimized}_{l,m} = \frac{N!^2}{(N + l)!(N -l)!} Y_{l,m} \f]
         with \f$N\f$ the order of decomposition, \f$l\f$ the degree and \f$m\f$ the order.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        virtual void process(T const* inputs, T* outputs) hoa_noexcept;
    };

    
    
    
    
    
    
    
    
    
    
    
    
    
    //! The optim class optimizes the ambisonic sound field for several restitution systems.
    /** The optim should be used to optimize the ambisonic sound field. There are 3 optimizations, Basic (no optimization), MaxRe (energy vector optimization) and InPhase (energy and velocity vector optimization). Basic has no effect, it should be used (or not) with a perfect ambisonic channels arrangement where all the channels are to equal distance on a circle or a sphere, and for a listener placed at the perfect center of the circle of the sphere. MaxRe should be used should be used for an auditory confined to the center of the circle of the sphere. InPhase should be used when the auditory covers the entire channels area and when the channels arrangement is not a perfect circle or a perfect sphere or when the channels are not to equal distance. Note that the optimizations decrease the precision of the sound field restitution thus it can be compared to particular cases of the fractional orders.
     */
    template <typename T> class Optim<Hoa2d, T> : public Processor<Hoa2d, T>::Harmonics
    {
    public:

        //! The optim constructor.
        /**	The optim constructor allocates and initialize the member values to computes spherical harmonics weighted coefficients depending on a order of decomposition. The order must be at least 1.
         @param     order	The order.
         */
        Optim(const size_t order) hoa_noexcept : Processor<Hoa2d, T>::Harmonics(order) {}

        //! The optim destructor.
        /**	The optim destructor free the memory.
         */
        virtual ~Optim() hoa_noexcept hoa_default_f

        //! This method performs the optimization.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        virtual void process(T const* inputs, T* outputs) hoa_noexcept = 0;
    };
    
    //! The basic optim.
    /** The basic optim has no effect, it should be used (or not) with a perfect ambisonic channels arrangement where all the channels are to equal distance on a circle or a sphere, and for a listener placed at the perfect center of the circle of the sphere.
     */
    template <typename T> class OptimBasic<Hoa2d, T> : public Optim<Hoa2d, T>
    {
    public:

        //! The optimization constructor.
        /**	The optimization constructor allocates and initialize the member values to computes spherical harmonics weighted coefficients depending on a order of decomposition. The order must be at least 1.
         @param     order	The order.
         */
        OptimBasic(const size_t order) hoa_noexcept :  Optim<Hoa2d, T>(order) {}

        //! The optimization destructor.
        /**	The optimization destructor free the memory.
         */
        ~OptimBasic() hoa_noexcept hoa_default_f

        //! This method performs the basic optimization.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        inline void process(T const* inputs, T* outputs) hoa_noexcept hoa_override
        {
            (*outputs++)  = (*inputs++);
            (*outputs++)  = (*inputs++);
            (*outputs++)  = (*inputs++);
            for(size_t i = 2; i <= Processor<Hoa2d, T>::Harmonics::getDecompositionOrder(); i++)
            {
                (*outputs++)  = (*inputs++);
                (*outputs++)  = (*inputs++);
            }
        }
    };

    //! The maxre optim.
    /** The maxre optim should be used for an auditory confined to the center of the circle of the sphere.
     */
    template <typename T> class OptimMaxRe<Hoa2d, T> : public Optim<Hoa2d, T>
    {
    private:
        static T* generate(const size_t order) hoa_noexcept
        {
            T* vector = Signal<T>::alloc(order);
            for(size_t i = 1; i <= order; i++)
            {
                vector[i-1] = cos(T(i) *  T(HOA_PI) / (T)(2. * order + 2.));
            }
            return vector;
        }
        const T*  m_weights;
    public:

        //! The optimization constructor.
        /**	The optimization constructor allocates and initialize the member values to computes spherical harmonics weighted coefficients depending on a order of decomposition. The order must be at least 1.
         @param     order	The order.
         */
        OptimMaxRe(const size_t order) hoa_noexcept : Optim<Hoa2d, T>(order),
        m_weights(generate(order)) {}

        //! The optimization destructor.
        /**	The optimization destructor free the memory.
         */
        ~OptimMaxRe() hoa_noexcept
        {
             Signal<T>::free(const_cast<T*>(m_weights));
        }

        //! This method performs the max-re optimization.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        inline void process(T const* inputs, T* outputs) hoa_noexcept hoa_override
        {
            const T* weights = m_weights;
            *outputs    = *inputs;
            *(++outputs)  = *(++inputs) * *weights;
            *(++outputs)  = *(++inputs) * *weights;
            for(size_t i = 2; i <= Processor<Hoa2d, T>::Harmonics::getDecompositionOrder(); i++)
            {
                const T weight = *(++weights);
                *(++outputs) = *(++inputs) * weight;
                *(++outputs) = *(++inputs) * weight;
            }
        }
    };

    //! The inphase optim.
    /** The inphase optim should be used when the auditory covers the entire channels area and when the channels arrangement is not a perfect circle or a perfect sphere or when the channels are not to equal distance.
     */
    template <typename T> class OptimInPhase<Hoa2d, T> : public  Optim<Hoa2d, T>
    {
    private:
        static T* generate(const size_t order)
        {
            T* vector = Signal<T>::alloc(order);
            const T facn = Math<T>::factorial(long(order));
            for(size_t i = 1; i <= order; i++)
            {
                vector[i-1] = facn / Math<T>::factorial(long(order - i)) * facn / Math<T>::factorial(long(order + i));
            }
            return vector;
        }
        const T*  m_weights;
    public:

        //! The optimization constructor.
        /**	The optimization constructor allocates and initialize the member values to computes spherical harmonics weighted coefficients depending on a order of decomposition. The order must be at least 1.
         @param     order	The order.
         */
        OptimInPhase(const size_t order) hoa_noexcept :  Optim<Hoa2d, T>(order),
        m_weights(generate(order)) {}

        //! The optimization destructor.
        /**	The optimization destructor free the memory.
         */
        ~OptimInPhase() hoa_noexcept
        {
             Signal<T>::free(const_cast<T*>(m_weights));
        }

        //! This method performs the in-phase optimization.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        inline void process(T const* inputs, T* outputs) hoa_noexcept hoa_override
        {
            const T* weights = m_weights;
            *outputs    = *inputs;
            *(++outputs)  = *(++inputs) * *weights;
            *(++outputs)  = *(++inputs) * *weights;
            for(size_t i = 2; i <= Processor<Hoa2d, T>::Harmonics::getDecompositionOrder(); i++)
            {
                const T weight = *(++weights);
                *(++outputs) = *(++inputs) * weight;
                *(++outputs) = *(++inputs) * weight;
            }
        }
    };

    
    
    
    
    
    template <typename T> class Optim<Hoa3d, T> : public Processor<Hoa3d, T>::Harmonics
    {
    public:

        //! The optim constructor.
        /**	The optim constructor allocates and initialize the member values to computes spherical harmonics weighted coefficients depending on a order of decomposition. The order must be at least 1.
         @param     order	The order.
         */
        Optim(const size_t order) hoa_noexcept : Processor<Hoa3d, T>::Harmonics(order) {}

        //! The optim destructor.
        /**	The optim destructor free the memory.
         */
        virtual ~Optim() hoa_noexcept {}

        //! This method performs the optimization.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        virtual void process(T const* inputs, T* outputs) hoa_noexcept = 0;
    };

    //! The basic optim.
    /** The basic optim has no effect, it should be used (or not) with a perfect ambisonic channels arrangement where all the channels are to equal distance on a circle or a sphere, and for a listener placed at the perfect center of the circle of the sphere.
     */
    template <typename T> class OptimBasic<Hoa3d, T> : public Optim<Hoa3d, T>
    {
    public:

        //! The optimization constructor.
        /**	The optimization constructor allocates and initialize the member values to computes spherical harmonics weighted coefficients depending on a order of decomposition. The order must be at least 1.
         @param     order	The order.
         */
        OptimBasic(const size_t order) hoa_noexcept : Optim<Hoa3d, T>(order) {}

        //! The optimization destructor.
        /**	The optimization destructor free the memory.
         */
        ~OptimBasic() hoa_noexcept hoa_default_f

        //! This method performs the basic optimization.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        inline void process(T const* inputs, T* outputs) hoa_noexcept hoa_override
        {
            (*outputs++)  = (*inputs++);
            (*outputs++)  = (*inputs++);
            (*outputs++)  = (*inputs++);
            (*outputs++)  = (*inputs++);
            for(size_t i = 2; i <= Processor<Hoa3d, T>::Harmonics::getDecompositionOrder(); i++)
            {
                for(size_t j = 0; j < 2 * i + 1; j++)
                {
                    (*outputs++)    = (*inputs++);    // Hamonic [i, ~j]
                }
            }
        }
    };

    //! The maxre optim.
    /** The maxre optim should be used for an auditory confined to the center of the circle of the sphere.
     */
    template <typename T> class OptimMaxRe<Hoa3d, T> : public Optim<Hoa3d, T>
    {
    private:
        static T* generate(const size_t order) hoa_noexcept
        {
            T* vector = Signal<T>::alloc(order);
            for(size_t i = 1; i <= order; i++)
            {
                vector[i-1] = cos(T(i) *  T(HOA_PI) / (T)(2. * order + 2.));
            }
            return vector;
        }
        const T*  m_weights;
    public:

        //! The optimization constructor.
        /**	The optimization constructor allocates and initialize the member values to computes spherical harmonics weighted coefficients depending on a order of decomposition. The order must be at least 1.
         @param     order	The order.
         */
        OptimMaxRe(const size_t order) hoa_noexcept : Optim<Hoa3d, T>(order),
        m_weights(generate(order)) {}

        //! The optimization destructor.
        /**	The optimization destructor free the memory.
         */
        ~OptimMaxRe() hoa_noexcept
        {
             Signal<T>::free(const_cast<T*>(m_weights));
        }

        //! This method performs the max-re optimization.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        inline void process(T const* inputs, T* outputs) hoa_noexcept hoa_override
        {
            const T* weights = m_weights;
            *outputs    = *inputs;
            *(++outputs)  = *(++inputs) * *weights;
            *(++outputs)  = *(++inputs) * *weights;
            *(++outputs)  = *(++inputs) * *weights;
            for(size_t i = 2; i <= Processor<Hoa3d, T>::Harmonics::getDecompositionOrder(); i++)
            {
                const T weight = *(++weights);
                for(size_t j = 0; j < 2 * i + 1; j++)
                {
                    *(++outputs) = *(++inputs) * weight;    // Hamonic [i, [-i...i]]
                }
            }
        }
    };

    //! The inphase optim.
    /** The inphase optim should be used when the auditory covers the entire channels area and when the channels arrangement is not a perfect circle or a perfect sphere or when the channels are not to equal distance.
     */
    template <typename T> class OptimInPhase<Hoa3d, T> : public Optim<Hoa3d, T>
    {
    private:
        static T* generate(const size_t order) hoa_noexcept
        {
            T* vector = Signal<T>::alloc(order);
            const T facn = Math<T>::factorial(long(order));
            for(size_t i = 1; i <= order; i++)
            {
                vector[i-1] = facn / Math<T>::factorial(long(order - i)) * facn / Math<T>::factorial(long(order + i));
            }
            return vector;
        }

        const T*  m_weights;
    public:

        //! The optimization constructor.
        /**	The optimization constructor allocates and initialize the member values to computes spherical harmonics weighted coefficients depending on a order of decomposition. The order must be at least 1.
         @param     order	The order.
         */
        OptimInPhase(const size_t order) hoa_noexcept : Optim<Hoa3d, T>(order),
        m_weights(generate(order)) {}

        //! The optimization destructor.
        /**	The optimization destructor free the memory.
         */
        ~OptimInPhase() hoa_noexcept
        {
            Signal<T>::free(const_cast<T*>(m_weights));
        }

        //! This method performs the in-phase optimization.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        inline void process(T const* inputs, T* outputs) hoa_noexcept hoa_override
        {
            const T* weights = m_weights;
            *outputs    = *inputs;
            *(++outputs)  = *(++inputs) * *weights;
            *(++outputs)  = *(++inputs) * *weights;
            *(++outputs)  = *(++inputs) * *weights;
            for(size_t i = 2; i <= Processor<Hoa3d, T>::Harmonics::getDecompositionOrder(); i++)
            {
                const T weight = *(++weights);
                for(size_t j = 0; j < 2 * i + 1; j++)
                {
                    *(++outputs) = *(++inputs) * weight;    // Hamonic [i, [-i...i]]
                }
            }
        }
    };
}

#endif
