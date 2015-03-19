/*
// Copyright (c) 2012-2015 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_OPTIM_LIGHT
#define DEF_HOA_OPTIM_LIGHT

#include "Harmonics.hpp"

namespace hoa
{
    //! The optim class optimizes the ambisonic sound field for several restitution systems.
    /** The optim should be used to optimize the ambisonic sound field. There are 3 optimizations, Basic (no optimization), MaxRe (energy vector optimization) and InPhase (energy and velocity vector optimization). Basic has no effect, it should be used (or not) with a perfect ambisonic channels arrengement where all the channels are to equal distance on a circle or a sphere, and for a listener placed at the perfect center of the circle ot the sphere. MaxRe should be used should be used for an auditory confined to the center of the circle ot the sphere. InPhase should be used when the auditory covers the entire channels area and when the channels arragement is not a perfect circle or a perfect sphere or when the channels are not to equal distance. Note that the optimizations decrease the precision of the sound field restitution thus it can be compared to particular cases of the fractional orders.
     */
    template <Dimension D, typename T> class Optim : public Processor< Harmonic<D, T> >
    {
    public:
        
        //! The optim constructor.
        /**	The optim constructor allocates and initialize the member values to computes spherical harmonics weighted coefficients depending of a decomposition order. The order must be at least 1.
         @param     order	The order.
         */
        Optim(const ulong order) noexcept = 0;
        
        //! The optim destructor.
        /**	The optim destructor free the memory.
         */
        virtual ~Optim() noexcept;
        
        //! This method performs the optimization.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        virtual void process(T const* inputs, T* outputs) const noexcept;
        
        //! The basic optim.
        /** The basic optim has no effect, it should be used (or not) with a perfect ambisonic channels arrengement where all the channels are to equal distance on a circle or a sphere, and for a listener placed at the perfect center of the circle ot the sphere.
         */
        class Basic : public Optim
        {
        public:
            
            //! The optim constructor.
            /**	The optim constructor allocates and initialize the member values to computes spherical harmonics weighted coefficients depending of a decomposition order. The order must be at least 1.
             @param     order	The order.
             */
            Basic(const ulong order) noexcept = 0;
            
            //! The optim destructor.
            /**	The optim destructor free the memory.
             */
            virtual ~Basic() noexcept;
            
            //! This method performs the basic optimization.
            /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
             \f[h'_{l,m} = h_{l,m}\f]
             with \f$l\f$ the degree and \f$m\f$ the order.
             @param     inputs   The inputs array.
             @param     outputs  The outputs array.
             */
            virtual void process(T const* inputs, T* outputs) const noexcept;
        };
        
        //! The maxre optim.
        /** The maxre optim should be used for an auditory confined to the center of the circle ot the sphere.
         */
        class MaxRe : public Optim
        {
        public:
            //! The optim constructor.
            /**	The optim constructor allocates and initialize the member values to computes spherical harmonics weighted coefficients depending of a decomposition order. The order must be at least 1.
             @param     order	The order.
             */
            MaxRe(const ulong order) noexcept = 0;
            
            //! The optim destructor.
            /**	The optim destructor free the memory.
             */
            virtual ~MaxRe() noexcept;
            
            //! This method performs the max-re optimization.
            /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
             \f[h'_{l,m} = \cos{(l \times \frac{\pi}{2N + 2})} h_{l,m} \f]
             with \f$N\f$ the order of decomposition, \f$l\f$ the degree and \f$m\f$ the order.
             @param     inputs   The inputs array.
             @param     outputs  The outputs array.
             */
            virtual void process(T const* inputs, T* outputs) const noexcept;
        };
        
        //! The inphase optim.
        /** The inphase optim should be used when the auditory covers the entire channels area and when the channels arragement is not a perfect circle or a perfect sphere or when the channels are not to equal distance.
         */
        class InPhase : public Optim
        {
        public:
            //! The optim constructor.
            /**	The optim constructor allocates and initialize the member values to computes spherical harmonics weighted coefficients depending of a decomposition order. The order must be at least 1.
             \f[h'_{l,m} = \frac{N!^2}{(N + l)!(N -l)!} h_{l,m} \f]
             with \f$N\f$ the order of decomposition, \f$l\f$ the degree and \f$m\f$ the order.
             @param     order	The order.
             */
            InPhase(const ulong order) noexcept = 0;
            
            //! The optim destructor.
            /**	The optim destructor free the memory.
             */
            virtual ~InPhase() noexcept;
            
            //! This method performs the in-phase optimization.
            /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
             @param     inputs   The inputs array.
             @param     outputs  The outputs array.
             */
            virtual void process(T const* inputs, T* outputs) const noexcept;
        };
    };
    
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    
    template <typename T> class Optim<Hoa2d, T> : public Processor< Harmonic<Hoa2d, T> >
    {
    public:
        
        //! The optim constructor.
        /**	The optim constructor allocates and initialize the member values to computes spherical harmonics weighted coefficients depending of a decomposition order. The order must be at least 1.
         @param     order	The order.
         */
        Optim(const ulong order) noexcept : Processor< Harmonic<Hoa2d, T> >(order)
        {
            ;
        }
        
        //! The optim destructor.
        /**	The optim destructor free the memory.
         */
        virtual ~Optim() noexcept
        {
            ;
        }
        
        //! This method performs the optimization.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        virtual void process(T const* inputs, T* outputs) const noexcept = 0;
        
        //! The basic optim.
        /** The basic optim has no effect, it should be used (or not) with a perfect ambisonic channels arrengement where all the channels are to equal distance on a circle or a sphere, and for a listener placed at the perfect center of the circle ot the sphere.
         */
        class Basic;
        
        //! The maxre optim.
        /** The maxre optim should be used for an auditory confined to the center of the circle ot the sphere.
         */
        class MaxRe;
        
        //! The inphase optim.
        /** The inphase optim should be used when the auditory covers the entire channels area and when the channels arragement is not a perfect circle or a perfect sphere or when the channels are not to equal distance.
         */
        class InPhase;
    };
    
    template <typename T> class Optim<Hoa2d, T>::Basic : public  Optim<Hoa2d, T>
    {
    public:
        
        //! The optimization constructor.
        /**	The optimization constructor allocates and initialize the member values to computes spherical harmonics weighted coefficients depending of a decomposition order. The order must be at least 1.
         @param     order	The order.
         */
        Basic(const ulong order) noexcept :  Optim<Hoa2d, T>(order)
        {
            ;
        }
        
        //! The optimization destructor.
        /**	The optimization destructor free the memory.
         */
        ~Basic() noexcept
        {
            ;
        }
        
        //! This method performs the basic optimization.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        inline void process(T const* inputs, T* outputs) const noexcept
        {
            (*outputs++)  = (*inputs++);
            (*outputs++)  = (*inputs++);
            (*outputs++)  = (*inputs++);
            for(ulong i = 2; i <= Processor< Harmonic<Hoa2d, T> >::getDecompositionOrder(); i++)
            {
                (*outputs++)  = (*inputs++);
                (*outputs++)  = (*inputs++);
            }
        }
    };
    
    template <typename T> class Optim<Hoa2d, T>::MaxRe : public  Optim<Hoa2d, T>
    {
    private:
        const T   m_cosmaxRe;
        const T   m_sinmaxRe;
    public:
        
        //! The optimization constructor.
        /**	The optimization constructor allocates and initialize the member values to computes spherical harmonics weighted coefficients depending of a decomposition order. The order must be at least 1.
         @param     order	The order.
         */
        MaxRe(const ulong order) noexcept :  Optim<Hoa2d, T>(order),
        m_cosmaxRe(cos(HOA_PI / (T)(2. * order + 2.))),
        m_sinmaxRe(sin(HOA_PI / (T)(2. * order + 2.)))
        {
            ;
        }
        
        //! The optimization destructor.
        /**	The optimization destructor free the memory.
         */
        ~MaxRe() noexcept
        {
            ;
        }
        
        //! This method performs the max-re optimization.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        inline void process(T const* inputs, T* outputs) const noexcept
        {
            T cos_re = m_cosmaxRe;
            T sin_re = m_sinmaxRe;
            T tcos_re = m_cosmaxRe;
            (*outputs++)  = (*inputs++);
            (*outputs++)  = (*inputs++) * cos_re;
            (*outputs++)  = (*inputs++) * cos_re;
            for(ulong i = 2; i <= Processor< Harmonic<Hoa2d, T> >::getDecompositionOrder(); i++)
            {
                cos_re  = tcos_re * m_cosmaxRe - sin_re * m_sinmaxRe;
                sin_re  = tcos_re * m_sinmaxRe + sin_re * m_cosmaxRe;
                tcos_re  = cos_re;
                
                (*outputs++)  = (*inputs++) * cos_re;
                (*outputs++)  = (*inputs++) * cos_re;
            }
        }
    };
    
    template <typename T> class Optim<Hoa2d, T>::InPhase : public  Optim<Hoa2d, T>
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
         */
        InPhase(const ulong order) noexcept :  Optim<Hoa2d, T>(order),
        m_facorder(Math<T>::factorial(order)),
        m_facinphase(m_facorder * m_facorder * order),
        m_facorder1(m_facorder * (order + 2) * order),
        m_facorder2(m_facorder / order)
        {
            ;
        }
        
        //! The optimization destructor.
        /**	The optimization destructor free the memory.
         */
        ~InPhase() noexcept
        {
            ;
        }
        
        //! This method performs the in-phase optimization.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        inline void process(T const* inputs, T* outputs) const noexcept
        {
            T order1 = Processor< Harmonic<Hoa2d, T> >::getDecompositionOrder() + 3;
            T order2 = Processor< Harmonic<Hoa2d, T> >::getDecompositionOrder() - 1;
            T factor1 = m_facorder1;
            T factor2 = m_facorder2;
            T factor  = m_facinphase / (factor1 * factor2);
            
            (*outputs++)  = (*inputs++);
            (*outputs++)  = (*inputs++) * factor;
            (*outputs++)  = (*inputs++) * factor;
            for(ulong i = 2; i <= Processor< Harmonic<Hoa2d, T> >::getDecompositionOrder(); i++)
            {
                factor1 *= order1++;
                factor2 /= order2--;
                factor = m_facinphase / (factor1 * factor2);
                
                (*outputs++)  = (*inputs++) * factor;
                (*outputs++)  = (*inputs++) * factor;
            }
        }
    };
    
    template <typename T> class Optim<Hoa3d, T> : public Processor< Harmonic<Hoa3d, T> >
    {
    public:
        
        //! The optim constructor.
        /**	The optim constructor allocates and initialize the member values to computes spherical harmonics weighted coefficients depending of a decomposition order. The order must be at least 1.
         @param     order	The order.
         */
        Optim(const ulong order) noexcept : Processor< Harmonic<Hoa3d, T> >(order)
        {
            ;
        }
        
        //! The optim destructor.
        /**	The optim destructor free the memory.
         */
        virtual ~Optim() noexcept
        {
            ;
        }
        
        //! This method performs the optimization.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        virtual void process(T const* inputs, T* outputs) const noexcept = 0;
        
        //! The basic optim.
        /** The basic optim has no effect, it should be used (or not) with a perfect ambisonic channels arrengement where all the channels are to equal distance on a circle or a sphere, and for a listener placed at the perfect center of the circle ot the sphere.
         */
        class Basic;
        
        //! The maxre optim.
        /** The maxre optim should be used for an auditory confined to the center of the circle ot the sphere.
         */
        class MaxRe;
        
        //! The inphase optim.
        /** The inphase optim should be used when the auditory covers the entire channels area and when the channels arragement is not a perfect circle or a perfect sphere or when the channels are not to equal distance.
         */
        class InPhase;
    };
    
    template <typename T> class Optim<Hoa3d, T>::Basic : public Optim<Hoa3d, T>
    {
    public:
        
        //! The optimization constructor.
        /**	The optimization constructor allocates and initialize the member values to computes spherical harmonics weighted coefficients depending of a decomposition order. The order must be at least 1.
         @param     order	The order.
         */
        Basic(const ulong order) noexcept : Optim<Hoa3d, T>(order)
        {
            ;
        }
        
        //! The optimization destructor.
        /**	The optimization destructor free the memory.
         */
        ~Basic() noexcept
        {
            ;
        }
        
        //! This method performs the basic optimization.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        inline void process(T const* inputs, T* outputs) const noexcept
        {
            (*outputs++)  = (*inputs++);
            (*outputs++)  = (*inputs++);
            (*outputs++)  = (*inputs++);
            (*outputs++)  = (*inputs++);
            for(ulong i = 2; i <= Processor< Harmonic<Hoa3d, T> >::getDecompositionOrder(); i++)
            {
                for(ulong j = 0; j < 2 * i + 1; j++)
                {
                    (*outputs++)    = (*inputs++);    // Hamonic [i, ~j]
                }
            }
        }
    };
    
    template <typename T> class Optim<Hoa3d, T>::MaxRe : public Optim<Hoa3d, T>
    {
    private:
        const T   m_cosmaxRe;
        const T   m_sinmaxRe;
    public:
        
        //! The optimization constructor.
        /**	The optimization constructor allocates and initialize the member values to computes spherical harmonics weighted coefficients depending of a decomposition order. The order must be at least 1.
         @param     order	The order.
         */
        MaxRe(const ulong order) noexcept : Optim<Hoa3d, T>(order),
        m_cosmaxRe(cos(HOA_PI / (T)(2. * order + 2.))),
        m_sinmaxRe(sin(HOA_PI / (T)(2. * order + 2.)))
        {
            ;
        }
        
        //! The optimization destructor.
        /**	The optimization destructor free the memory.
         */
        ~MaxRe() noexcept
        {
            ;
        }
        
        //! This method performs the max-re optimization.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        inline void process(T const* inputs, T* outputs) const noexcept
        {
            T cos_re = m_cosmaxRe;
            T sin_re = m_cosmaxRe;
            T tcos_re= cos_re;
            (*outputs++)  = (*inputs++);
            (*outputs++)  = (*inputs++) * cos_re;
            (*outputs++)  = (*inputs++) * cos_re;
            (*outputs++)  = (*inputs++) * cos_re;
            for(ulong i = 2; i <= Processor< Harmonic<Hoa3d, T> >::getDecompositionOrder(); i++)
            {
                cos_re  = tcos_re * m_cosmaxRe - sin_re * m_sinmaxRe;
                sin_re  = tcos_re * m_sinmaxRe + sin_re * m_cosmaxRe;
                tcos_re  = cos_re;
                for(ulong j = 0; j < 2 * i + 1; j++)
                {
                    (*outputs++)    = (*inputs++) * cos_re;    // Hamonic [i, [-i...i]]
                }
            }
        }
    };

    template <typename T> class Optim<Hoa3d, T>::InPhase : public Optim<Hoa3d, T>
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
         */
        InPhase(const ulong order) noexcept : Optim<Hoa3d, T>(order),
        m_facorder(Math<T>::factorial(order)),
        m_facinphase(m_facorder * m_facorder * order),
        m_facorder1(m_facorder * (order + 2) * order),
        m_facorder2(m_facorder / order)
        {
            ;
        }
        
        //! The optimization destructor.
        /**	The optimization destructor free the memory.
         */
        ~InPhase() noexcept
        {
            ;
        }
        
        //! This method performs the in-phase optimization.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        inline void process(T const* inputs, T* outputs) const noexcept
        {
            T order1 = Processor< Harmonic<Hoa3d, T> >::getDecompositionOrder() + 3;
            T order2 = Processor< Harmonic<Hoa3d, T> >::getDecompositionOrder() - 1;
            T factor1 = m_facorder1;
            T factor2 = m_facorder2;
            T factor  = m_facinphase / (factor1 * factor2);
            
            (*outputs++)  = (*inputs++);
            (*outputs++)  = (*inputs++) * factor;
            (*outputs++)  = (*inputs++) * factor;
            (*outputs++)  = (*inputs++) * factor;
            for(ulong i = 2; i <= Processor< Harmonic<Hoa3d, T> >::getDecompositionOrder(); i++)
            {
                factor1 *= order1++;
                factor2 /= order2--;
                factor = m_facinphase / (factor1 * factor2);
                
                for(ulong j = 0; j < 2 * i + 1; j++)
                {
                    (*outputs++)    = (*inputs++) * factor;    // Hamonic [i, ~j]
                }
            }
        }
    };
    
#endif
}

#endif



